#!/usr/bin/env python3
"""Verify the plan hash chains without any external timestamp.

The two plan files were hashed repeatedly as they were edited. This script
checks that the text we hold is the text that was hashed, and that the later
edits of PREREG_v2.md were appends rather than rewrites:

  * hash the file as it stands and look it up in the chain;
  * for PREREG_v2.md, cut the file at the start of each appended supplement and
    hash the prefix, then look that up too.

A prefix of the final text hashing to an earlier chain entry proves that
everything before that point was already there, byte for byte, when the hash was
recorded. It does not prove when the hash was computed: the timestamps come from
local files and have no external anchor.

Usage: python build/verify_hash_chain.py
"""
import hashlib
import io

def chain(path):
    out = []
    for line in io.open(path, encoding="utf-8"):
        head = line.split()[0]
        if len(head) == 64:
            out.append(head)
    return out


def sha(b):
    if isinstance(b, str):
        b = b.encode("utf-8")
    return hashlib.sha256(b).hexdigest()


def audit(md, hashfile, cuts):
    text = io.open(md, "rb").read().decode("utf-8")   # keep CRLF, truncate on characters
    ch = chain(hashfile)
    print(f"\n{md}  ({len(ch)} entries in {hashfile})")
    h = sha(text.rstrip("\n") + "\n")
    print(f"  whole file            {h[:12]}  {'= chain[%d]' % ch.index(h) if h in ch else 'NOT IN CHAIN'}")
    for label, anchor in cuts:
        i = text.find(anchor)
        if i < 0:
            print(f"  prefix before {label:<9} anchor not found")
            continue
        h = sha(text[:i].rstrip("\n") + "\n")
        state = f"= chain[{ch.index(h)}]" if h in ch else "not in chain"
        print(f"  prefix before {label:<9} {h[:12]}  {state}")


audit("PREREG_v2.md", "prereg_v2_hash.txt",
      [("§7", "## 7. 补充 A"), ("§8", "## 8. 补充 B"), ("§9", "## 9. 补充 C"),
       ("§10", "## 10. 补充 D"), ("§11", "## 11. 补充 E")])
audit("PREREG_fairfax_main.md", "prereg_hash.txt", [("appendix", "## 补充条款")])
print("\nA prefix match shows the earlier text was not rewritten after hashing.")
print("It does not show when the hash was made: the times are local file records.")
