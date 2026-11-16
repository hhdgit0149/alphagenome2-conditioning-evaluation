# Input data

None of the raw data is redistributed here. Download from the sources below; the scripts expect
the files inside `b_data/` (and `mpra/`) at the repository root, see the top-level README.

| What | Source | Accession / path |
|---|---|---|
| Single-cell PBMC, control and IFN-beta-stimulated, 8 donors | Kang et al. 2018, Nat Biotechnol | GEO **GSE96583**, batch 2 (`GSM2560248`, `GSM2560249`), plus the tSNE/cell-type table |
| Nominal summary statistics and credible sets, monocytes | eQTL Catalogue r8 | `https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/` — QTS000010 (Fairfax) and QTS000024 (Quach), studies QTD000081–84 and QTD000409/414/419/424/429 |
| Fairfax naive / IFN-gamma / LPS expression (for the match score) | ArrayExpress | **E-MTAB-2232**, processed files per condition |
| U2OS-GR dexamethasone time course (for the scale check) | GEO | **GSE46448** |
| AR enhancer activity regions (element level, extension only) | GEO | **GSE151064** |
| Glucocorticoid MPRA (extension only) | PNAS 2023 supplement | `pnas.2305773120.sd01.xlsx`, sheets S1–S3 |

## Notes

- The eQTL files are large (tens of GB in total). `35_strict_z_lookup.py` avoids downloading them by
  issuing throttled HTTP range requests against the `.tbi` index; it needs only the index files,
  which go in `eqtl_cache/`.
- Expression values must be on a linear scale before they are used as a conditioning vector.
  Log-scale microarray values produce predictions that look inert (see the manuscript's Methods).
- Respect the rate limits of the EBI FTP service: the client throttles to one request every 1.4 s
  and was run serially.
