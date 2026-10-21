#!/usr/bin/env python3
"""Validate organism concordance and emit an explicit pass/fail result.

Thin entrypoint delegating to the ``validate-organism-parity`` CLI command,
which itself delegates to :mod:`mirpathdb_prod.validation.organism_parity`.
Compares by canonical taxid under the "no-loss FLOOR + expand" model: the
miRPathDB2 floor (hsa, mmu) must be preserved, and the production set may expand
with the approved extensions rat (10116, validated+predicted) and macaque (9544,
predicted-dominant — lone validated MTI marked).
"""

from __future__ import annotations

import sys

from mirpathdb_prod.cli import main as cli_main


def main() -> int:
    return cli_main(["validate-organism-parity", *sys.argv[1:]])


if __name__ == "__main__":
    sys.exit(main())
