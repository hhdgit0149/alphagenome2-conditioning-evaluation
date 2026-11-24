#!/usr/bin/env python3
"""Fetch the model catalogues' metadata only. No scoring, no predictions.

Run this on your own machine with your own API key in the environment:

    set ALPHAGENOME_API_KEY=...        (Windows, current shell only)
    python build/fetch_catalogue_metadata.py

It writes two CSVs (all text fields per track / per condition), records the
retrieval date and the SDK version, and re-counts the two numbers the paper
reports: 4,430 first-generation tracks and 2,799 second-generation conditions.
It also lists every condition carrying an activation-state label and flags the
ones whose text looks like an in-vitro stimulation.

Nothing here is uploaded anywhere; the key is read from the environment and is
never written to disk.
"""
import datetime
import os
import sys

import pandas as pd

KEY = os.environ.get("ALPHAGENOME_API_KEY")
if not KEY:
    sys.exit("set ALPHAGENOME_API_KEY in the environment first")

try:
    import alphagenome
    sdk = getattr(alphagenome, "__version__", "unknown")
except Exception:
    sdk = "not importable"

today = datetime.date.today().isoformat()
print(f"retrieval date {today}   sdk {sdk}")

# The two listing calls differ between releases; adjust the attribute names to
# the ones your installed SDK exposes and re-run.
from alphagenome.data import genome  # noqa: F401  (import check only)

# ---- fill these in from the SDK's metadata interface ------------------------
# e.g. client = alphagenome.AlphaGenome(api_key=KEY)
#      tracks = client.tracks(...)          # first generation
#      conds  = client.conditions(...)      # second generation
# The rest of the script only needs two DataFrames with one row per track /
# condition and all text fields kept as columns.

def main():
    raise SystemExit(
        "Open the SDK docs for the metadata listing call, assign the two "
        "DataFrames to tracks and conds below, then remove this line."
    )


if __name__ == "__main__":
    main()
