#!/usr/bin/env python3
"""严格特异定义的 z 查询（PREREG_v2 补充C §1）：按分子性状逐行记录。

手写极简 tabix 远程客户端：.tbi 线性索引 + HTTP Range 请求 + BGZF 解块。
每个 (other_cond, variant, gene) 查询输出：
  - 该基因在该变异处的每个分子性状一行（trait, beta, se, median_tpm）
  - 位置级标志：n_pos_rows（该位置全部行数）、variant_tested（变异是否被测过）
  - 基因无命中时写一行 trait=""、n_pair_rows=0
  - 查询失败写一行 trait=""、n_pair_rows=-1（续跑时重查）
节流 >=1.4s/请求（EBI 防火墙会把频繁请求当 DoS），断点续跑，单写者。
"""
import gzip
import os
import struct
import sys
import time
import zlib

import pandas as pd
import requests

BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/{qts}/{qtd}/{qtd}"
TBI_CACHE = "eqtl_cache/sumstats"
THROTTLE = 1.4
BLOCK_FETCH = 131072
MAX_ROUNDS = 24

STUDIES = {
    "fairfax": dict(
        qts="QTS000010",
        cond={"naive": "QTD000081", "IFNg_24h": "QTD000082",
              "LPS_2h": "QTD000083", "LPS_24h": "QTD000084"},
        labels="b_data/fairfax_labels_fixed.csv",
        out="b_data/strict_z_fairfax.csv"),
    "quach": dict(
        qts="QTS000024",
        cond={"naive": "QTD000409", "LPS_6h": "QTD000414",
              "Pam3CSK4_6h": "QTD000419", "R848_6h": "QTD000424",
              "Influenza_6h": "QTD000429"},
        labels="b_data/Quach_2016_labels_fixed.csv",
        out="b_data/strict_z_quach.csv"),
}

session = requests.Session()
last_req = [0.0]


def throttle():
    dt = time.time() - last_req[0]
    if dt < THROTTLE:
        time.sleep(THROTTLE - dt)
    last_req[0] = time.time()


def http_range(url, start, end):
    """end 含端点。失败重试 3 次；带 cache-buster 绕过本地代理的坏缓存。"""
    for att in range(3):
        throttle()
        u = f"{url}?cb={time.time_ns()}"
        try:
            r = session.get(u, headers={"Range": f"bytes={start}-{end}",
                                        "Cache-Control": "no-cache"},
                            timeout=(20, 60))
            if r.status_code in (200, 206) and r.content:
                return r.content
            print(f"  HTTP {r.status_code} 空/异常，重试", flush=True)
        except Exception as e:
            print(f"  网络异常 {type(e).__name__}: {e}，重试", flush=True)
        time.sleep(10)
    raise RuntimeError(f"重试用尽: {url} {start}-{end}")


def fetch_whole(url, path):
    for att in range(5):
        throttle()
        try:
            r = session.get(f"{url}?cb={time.time_ns()}",
                            headers={"Cache-Control": "no-cache"},
                            timeout=(20, 180))
            if r.status_code == 200 and len(r.content) > 100:
                with open(path, "wb") as f:
                    f.write(r.content)
                return
        except Exception:
            pass
        time.sleep(15 * (att + 1))
    raise RuntimeError(f"下载失败: {url}")


class Tbi:
    """最小 tabix 索引：参考名列表 + 每条染色体的线性索引。"""

    def __init__(self, path):
        raw = gzip.decompress(open(path, "rb").read())
        assert raw[:4] == b"TBI\x01", "非 tabix 索引"
        (self.n_ref, self.fmt, self.col_seq, self.col_beg, self.col_end,
         self.meta, self.skip) = struct.unpack_from("<7i", raw, 4)
        (l_nm,) = struct.unpack_from("<i", raw, 32)
        off = 36
        self.names = raw[off:off + l_nm].decode().split("\0")[:-1]
        off += l_nm
        self.lin = []
        for _ in range(self.n_ref):
            (n_bin,) = struct.unpack_from("<i", raw, off)
            off += 4
            for _ in range(n_bin):
                (n_chunk,) = struct.unpack_from("<i", raw, off + 4)
                off += 8 + 16 * n_chunk
            (n_intv,) = struct.unpack_from("<i", raw, off)
            off += 4
            lin = struct.unpack_from(f"<{n_intv}Q", raw, off) if n_intv else ()
            self.lin.append(lin)
            off += 8 * n_intv


def bgzf_blocks(buf, url="", start=0):
    """切出 buf 里所有完整 BGZF 块，返回 (消耗字节数, 解压文本列表)。"""
    out, i, n = [], 0, len(buf)
    while i + 18 <= n:
        if buf[i:i + 2] != b"\x1f\x8b":
            raise RuntimeError(
                f"非BGZF头 @{start+i} 前16字节={buf[i:i+16].hex()} len={n}")
        xlen = struct.unpack_from("<H", buf, i + 10)[0]
        j, bsize = i + 12, None
        while j < i + 12 + xlen:
            si1, si2 = buf[j], buf[j + 1]
            slen = struct.unpack_from("<H", buf, j + 2)[0]
            if si1 == 66 and si2 == 67:
                bsize = struct.unpack_from("<H", buf, j + 4)[0] + 1
            j += 4 + slen
        if bsize is None or i + bsize > n:
            break
        block = buf[i:i + bsize]
        d = zlib.decompressobj(16 + zlib.MAX_WBITS)
        try:
            out.append(d.decompress(block).decode("utf-8", "replace"))
        except zlib.error as e:
            raise RuntimeError(
                f"zlib解压失败 块@{start + i} bsize={bsize} msg={e} "
                f"头16字节={block[:16].hex()} resp_len={n}") from e
        i += bsize
    return i, out


class RemoteTabix:
    def __init__(self, qts, qtd):
        self.qtd = qtd
        self.url = BASE.format(qts=qts, qtd=qtd) + ".all.tsv.gz"
        os.makedirs(TBI_CACHE, exist_ok=True)
        tbi = f"{TBI_CACHE}/{qtd}.all.tsv.gz.tbi"
        if not os.path.exists(tbi):
            fetch_whole(self.url + ".tbi", tbi)
        self.idx = Tbi(tbi)
        first = http_range(self.url, 0, BLOCK_FETCH - 1)
        _, texts = bgzf_blocks(first, self.url, 0)
        head = texts[0].split("\n")[0].replace("\r", "").split("\t")
        self.cols = {c: k for k, c in enumerate(head)}
        for need in ("variant", "gene_id", "beta", "se",
                     "molecular_trait_id", "median_tpm"):
            assert need in self.cols, f"{qtd} 缺列 {need}"

    def query_point(self, chrom, pos):
        """取 (chrom, pos) 的全部行（所有基因），返回 list[dict]。"""
        pos = int(pos)
        if chrom not in self.idx.names:
            return []
        ri = self.idx.names.index(chrom)
        lin = self.idx.lin[ri]
        tile = (pos - 1) >> 14
        if tile >= len(lin) or lin[tile] == 0:
            return []
        cursor = lin[tile] >> 16
        pos_col = self.idx.col_beg - 1
        chr_col = self.idx.col_seq - 1
        pending = ""
        tail = b""
        for _ in range(MAX_ROUNDS):
            resp = http_range(self.url, cursor, cursor + BLOCK_FETCH - 1)
            chunk = tail + resp
            used, blocks = bgzf_blocks(chunk, self.url, cursor - len(tail))
            tail = chunk[used:]
            cursor += len(resp)
            if not blocks:
                return []
            rows = []
            for text in blocks:
                pending += text
                lines = pending.split("\n")
                pending = lines[-1]
                for line in lines[:-1]:
                    f = line.rstrip("\r").split("\t")
                    if len(f) <= pos_col:
                        continue
                    try:
                        p = int(f[pos_col])
                    except ValueError:
                        continue
                    if p > pos:
                        return self._to_dicts(rows)
                    if p == pos and f[chr_col] == chrom:
                        rows.append(f)
        raise RuntimeError(f"6×128KB 内未越过目标位置 {chrom}:{pos}，数据可疑")

    def _to_dicts(self, rows):
        c = self.cols
        return [dict(variant=r[c["variant"]], gene_id=r[c["gene_id"]],
                     trait=r[c["molecular_trait_id"]],
                     beta=float(r[c["beta"]]), se=float(r[c["se"]]),
                     tpm=r[c["median_tpm"]])
                for r in rows]


def query_window(tab, chrom, pos, variant, tiles=3):
    """扩窗查询：±tiles 个 16kb tile 内是否存在该变异（供未检出复核，精确匹配变异串）。"""
    pos = int(pos)
    lo = max(1, pos - tiles * 16384)
    hi = pos + tiles * 16384
    ri = tab.idx.names.index(chrom)
    lin = tab.idx.lin[ri]
    start_tile = (lo - 1) >> 14
    while start_tile < len(lin) and lin[start_tile] == 0:
        start_tile += 1
    if start_tile >= len(lin):
        return False
    cursor = lin[start_tile] >> 16
    pos_col = tab.idx.col_beg - 1
    chr_col = tab.idx.col_seq - 1
    vcol = tab.cols["variant"]
    pending, tail = "", b""
    for _ in range(24):
        resp = http_range(tab.url, cursor, cursor + BLOCK_FETCH - 1)
        chunk = tail + resp
        used, blocks = bgzf_blocks(chunk, tab.url, cursor - len(tail))
        tail = chunk[used:]
        cursor += len(resp)
        if not blocks:
            break
        stop = False
        for text in blocks:
            pending += text
            lines = pending.split("\n")
            pending = lines[-1]
            for line in lines[:-1]:
                f = line.rstrip("\r").split("\t")
                if len(f) <= pos_col:
                    continue
                try:
                    p = int(f[pos_col])
                except ValueError:
                    continue
                if p > hi:
                    stop = True
                    break
                if lo <= p <= hi and f[chr_col] == chrom and f[vcol] == variant:
                    return True
            if stop:
                break
    return False


def clean_resume(path):
    """失败行(n_pair_rows=-1)重查后旧失败行残留 → 清理：有成功行的对丢弃失败行。"""
    prev = pd.read_csv(path)
    ok_pairs = set(map(tuple, prev[prev.n_pair_rows >= 0]
                       [["other_cond", "variant", "gene"]].values))
    prev = prev[[t in ok_pairs for t in
                 zip(prev.other_cond, prev.variant, prev.gene)]].copy()
    assert not prev.duplicated(
        ["other_cond", "variant", "gene", "trait", "row_id"]).any()
    prev.to_csv(path, index=False)
    return prev


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    todo = list(STUDIES) if which == "both" else [which]
    for key in todo:
        cfg = STUDIES[key]
        lab = pd.read_csv(cfg["labels"])
        need = {}
        for r in lab.itertuples(index=False):
            if r.group == "shared":
                continue
            other = "naive" if r.group == "act_specific" else r.stim
            need[(other, r.variant, r.gene)] = (other, r.variant, r.gene,
                                                r.chrom, int(r.pos))
        done = set()
        if os.path.exists(cfg["out"]):
            prev = clean_resume(cfg["out"])
            done = set(map(tuple, prev[["other_cond", "variant", "gene"]]
                           .drop_duplicates().values))
        tasks = [v for k, v in need.items() if k not in done]
        print(f"[{key}] 需查 {len(need)}，已完成 {len(done)}，待查 {len(tasks)}",
              flush=True)
        if not tasks:
            continue
        by_cond = {}
        for t in tasks:
            by_cond.setdefault(t[0], []).append(t)
        header = not os.path.exists(cfg["out"])
        t0, n = time.time(), 0
        for cond, items in by_cond.items():
            items.sort(key=lambda x: (x[3], x[4]))
            tab = RemoteTabix(cfg["qts"], cfg["cond"][cond])
            for other_cond, variant, gene, chrom, pos in items:
                cname = chrom[3:] if chrom.startswith("chr") else chrom
                pos_rows = None
                for att in range(2):
                    try:
                        pos_rows = tab.query_point(cname, pos)
                        break
                    except Exception as e:
                        print(f"  {variant}/{cond} 第{att+1}次失败: {e}",
                              flush=True)
                        time.sleep(15)
                if pos_rows is None:                     # 查询失败，续跑重查
                    rec = [dict(other_cond=other_cond, variant=variant, gene=gene,
                                chrom=chrom, pos=pos, trait="", beta="", se="",
                                tpm="", n_pair_rows=-1, n_pos_rows=-1,
                                variant_tested="", row_id=0)]
                else:
                    tested = any(r["variant"] == variant for r in pos_rows)
                    pair = [r for r in pos_rows if r["variant"] == variant
                            and r["gene_id"] == gene]
                    seen, uniq = set(), []
                    for r in pair:                       # rsid 重复行去重
                        k2 = (r["trait"], r["beta"], r["se"])
                        if k2 not in seen:
                            seen.add(k2)
                            uniq.append(r)
                    if uniq:
                        rec = [dict(other_cond=other_cond, variant=variant,
                                    gene=gene, chrom=chrom, pos=pos,
                                    trait=r["trait"], beta=r["beta"], se=r["se"],
                                    tpm=r["tpm"], n_pair_rows=len(uniq),
                                    n_pos_rows=len(pos_rows),
                                    variant_tested=tested, row_id=i)
                               for i, r in enumerate(sorted(
                                   uniq, key=lambda x: (x["trait"], x["beta"])))]
                    else:
                        rec = [dict(other_cond=other_cond, variant=variant,
                                    gene=gene, chrom=chrom, pos=pos, trait="",
                                    beta="", se="", tpm="", n_pair_rows=0,
                                    n_pos_rows=len(pos_rows),
                                    variant_tested=tested, row_id=0)]
                pd.DataFrame(rec).to_csv(cfg["out"], mode="a", index=False,
                                         header=header)
                header = False
                n += 1
                if n % 50 == 0:
                    el = time.time() - t0
                    print(f"  [{key}] {n}/{len(tasks)}｜{n/max(el,1)*60:.0f} 查/分",
                          flush=True)
        print(f"[{key}] 完成 → {cfg['out']}", flush=True)


if __name__ == "__main__":
    main()
