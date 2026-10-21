#!/usr/bin/env python3
"""Extract the prototype/legacy organism baseline.

Thin entrypoint delegating to the ``extract-legacy-organisms`` CLI command,
which reads ``validation/organisms_old.tsv`` via
:mod:`mirpathdb_prod.sources.legacy` (read-only; never imports prototype code).
"""

from __future__ import annotations

import sys

from mirpathdb_prod.cli import main as cli_main


def main() -> int:
    return cli_main(["extract-legacy-organisms", *sys.argv[1:]])


if __name__ == "__main__":
    sys.exit(main())
