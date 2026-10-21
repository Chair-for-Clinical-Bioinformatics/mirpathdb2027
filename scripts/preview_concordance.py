#!/usr/bin/env python
"""Quick cell-for-cell concordance preview: our matrices vs the frozen v2 ones.

This is a lightweight PREVIEW for sanity-checking the controlled reproduction
output (it is NOT the full validation report — that is the validator's job, with
figures and BH-FDR agreement). For each produced cell it prints, on the cells
tested in BOTH matrices (p < 1.0), the Spearman rho of ``-log10 p``, the median
absolute ``-log10 p`` difference, and the raw-p<=0.05 directional agreement.

By construction our matrices share v2's exact miRNA-row and pathway-column order,
so this asserts that alignment first (a HARD concordance check on the comparison
frame) and fails loudly if it ever drifts.

Usage::

    python scripts/preview_concordance.py --organisms hsa mmu
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
from scipy.stats import spearmanr

from mirpathdb_prod.config import load_config
from mirpathdb_prod.reproduction.runner import VARIANTS
from mirpathdb_prod.reproduction.v2_matrix import read_v2_matrix, v2_matrix_path

FAMILIES = {
    "hsa": ("GO_BP", "GO_CC", "GO_MF", "REACTOME", "WIKIPATHWAYS"),
    "mmu": ("GO_BP", "GO_CC", "GO_MF", "WIKIPATHWAYS"),
}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--organisms", nargs="+", default=["mmu", "hsa"])
    args = parser.parse_args(argv)

    cfg = load_config()
    v2_root = cfg.path("dragon_backup_v2_dir")
    out_root = cfg.path("data_processed_dir") / "enrichment"

    hdr = f"{'cell':48s} {'both':>8s} {'rho':>7s} {'medΔlog':>8s} {'agree%':>7s}"
    print(hdr)
    print("-" * len(hdr))
    for organism in args.organisms:
        for variant in VARIANTS:
            for family in FAMILIES[organism]:
                new_path = out_root / organism / f"{family}_{variant}.csv"
                if not new_path.exists():
                    continue
                v2p = v2_matrix_path(v2_root, organism, family, variant)
                nm = read_v2_matrix(new_path)
                vm = read_v2_matrix(v2p)
                # HARD concordance check on the comparison frame.
                assert nm.pathways == vm.pathways, f"{organism}/{family}/{variant}: pathway drift"
                assert nm.mirnas == vm.mirnas, f"{organism}/{family}/{variant}: miRNA drift"
                new, v2 = nm.values, vm.values
                both = (new < 1.0) & (v2 < 1.0)
                label = f"{organism}/{family}/{variant}"
                if both.sum() < 5:
                    print(f"{label:48s} {both.sum():8d}   (too few both-tested)")
                    continue
                a = -np.log10(new[both])
                b = -np.log10(v2[both])
                rho, _ = spearmanr(a, b)
                med = float(np.median(np.abs(a - b)))
                tested = (new < 1.0) | (v2 < 1.0)
                agree = float(np.mean(((new <= 0.05) == (v2 <= 0.05))[tested]))
                print(f"{label:48s} {both.sum():8d} {rho:7.3f} {med:8.3f} {100 * agree:6.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
