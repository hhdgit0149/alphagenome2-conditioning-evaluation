# Evaluation code for "Cell-state conditioning in AlphaGenome 2 reproduces an unseen interferon response but shows no reliable stimulus-specific eQTL effects"

Archived at **https://doi.org/10.5281/zenodo.24494903 (v1.1.0)**, together with the derived tables
and the manuscript sources. This repository holds the code only.

Licences: the code is MIT. Files in the archive that contain AlphaGenome outputs, or numbers
computed from them, remain subject to the AlphaGenome Output Terms of Use
(http://deepmind.google.com/science/alphagenome/output-terms); see the archive README for the
file-by-file list.

This repository contains the analysis scripts behind the manuscript. It does **not** contain
the input data (several GB, third-party, different redistribution terms); the sources and the
exact accessions are listed in `data/README.md`.

## What the code does

| Script | Purpose |
|---|---|
| `17_build_B.py` | Builds the conditioning AnnData objects (256 cells per condition) and the pseudobulk profile |
| `18_B_predict.py` | Scores 40 one-megabase windows per condition with the denoised single-cell expression head |
| `20_B_leaveout.py` | Gene hold-out: removes the 119 target genes from the conditioning input |
| `54_disjoint_obs.py` | Recomputes the measured response on cells disjoint from the conditioning input |
| `57_calibration.py`, `59_calib_middle.py`, `60_calib_middle_score.py`, `61_reliability.py`, `65_rebuild15.py` | The 15-point calibration, its nulls and the reliability correction |
| `35_strict_z_lookup.py` | Minimal remote tabix client: parses `.tbi`, issues HTTP range requests, decodes BGZF blocks |
| `36_strict_labels.py`, `38_validate_client.py` | Strict condition-specific labels; validation of the tabix client against the credible sets |
| `42_head_sensitivity.py`, `47_flu_sensitivity.py`, `48_strict19.py` | Head sensitivity, leave-one-out, visibility gate, threshold sensitivity |
| `44_arm4_mismatch.py`, `45_build_mismatch_cond.py`, `46_arm3_v1.py` | Mismatched-condition arm and the first-generation reference-track arm |
| `52_null_tb.py`, `53_noise_splits.py`, `56_null_near.py` | Split-half null and magnitude-matched nulls |
| `62_condition_marker_check.py` | Three-layer check that condition labels are not swapped |
| `63_export_for_audit.py`, `66_export_inventory.py`, `67_marker_genes.py`, `68_keyword_table.py` | Export the per-gene tables, inventories and keyword counts |
| `69_figures.py` | All six figures, from archived CSVs, without calling the model |
| `build_pdf.py`, `make_supp_tables.py`, `make_supp_pandoc.py`, `make_references.py` | Typesetting and table generation |

## Layout expected at run time

Every script uses relative paths and expects to run from a working directory containing:

```
b_data/        per-gene prediction tables, credible sets, count matrices, cached AnnData
eqtl_cache/    downloaded .tbi indexes
mpra/          MPRA supplementary files (extension only, not used in v1)
```

Set the API key before any scoring script:

```bash
export ALPHAGENOME_API_KEY=...      # no key is stored in this repository
```

## Environment

Python 3.13, numpy 2.5, pandas 3.0, scipy 1.18, anndata 0.13, matplotlib 3.11, `alphagenome` 2.0.0.
Versions actually used for the reported numbers are in `requirements.txt`.

## Reproducing the reported numbers without re-scoring

The model calls are the expensive part and are not needed to check the statistics. The archived
per-gene tables in `tables/` and the supplementary CSVs are enough to recompute:

- the 15-point calibration (`calibration_final.csv`)
- the split-half and magnitude-matched nulls
- the attrition counts in Table S1
- the keyword counts in Table S3

`69_figures.py` reads only these files.

## Data

See `data/README.md`. Nothing in this repository is derived from data with redistribution
restrictions beyond the terms of the original sources.

## Licence

MIT — see `LICENSE`. The archived tables under `tables/` and the manuscript text are
additionally available under CC-BY-4.0.
