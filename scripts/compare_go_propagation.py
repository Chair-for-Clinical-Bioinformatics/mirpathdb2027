#!/usr/bin/env python3
"""GO ``all`` vs ``direct`` propagation sensitivity comparison.

Thin entrypoint delegating to the ``compare-go-propagation`` CLI command, which
itself delegates to :mod:`mirpathdb_prod.validation.go_propagation_driver`. The
production ``all`` (GOALL) GO enrichment is reused when its build fingerprint
matches; the ``direct`` (GO) side is computed in an isolated directory with the
production loaders + clusterProfiler kernel. Nothing under ``data/processed/`` or
``config.yaml`` is modified.

Usage (from the repo root)::

    MPD_RSCRIPT=.../envs/mpd_cp/bin/Rscript PYTHONPATH=src \\
      python scripts/compare_go_propagation.py [-c config.yaml] \\
        [--organisms hsa mmu rno mml] [--tiers HIGH MEDIUM LOW] \\
        [--families GO_BP GO_CC GO_MF] [--out-dir DIR] [--fdr 0.05] \\
        [--secondary-fdr 0.005] [--structural-only] [--recompute-all] \\
        [--reuse-only] [--max-cores N]

``-c/--config`` is a top-level option of ``mirpathdb-prod``; this wrapper hoists
it in front of the subcommand so it can be given anywhere on the command line.
"""

from __future__ import annotations

import sys

from mirpathdb_prod.cli import main as cli_main


def _split_config(argv: list[str]) -> tuple[list[str], list[str]]:
    """Return ``(top_level, rest)``: ``-c/--config <path>`` hoisted before the subcommand."""
    top: list[str] = []
    rest: list[str] = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a in ("-c", "--config") and i + 1 < len(argv):
            top.extend([a, argv[i + 1]])
            i += 2
            continue
        if a.startswith("--config="):
            top.append(a)
            i += 1
            continue
        rest.append(a)
        i += 1
    return top, rest


def main() -> int:
    top, rest = _split_config(sys.argv[1:])
    return cli_main([*top, "compare-go-propagation", *rest])


if __name__ == "__main__":
    sys.exit(main())
