#!/usr/bin/env python3
"""Fetch/load the miRPathDB2 reference (primary external validation source).

Thin entrypoint delegating to the ``build-mirpathdb2-reference`` CLI command,
which reads ``validation/organisms_mirpathdb2.tsv`` via
:mod:`mirpathdb_prod.sources.mirpathdb2`.
"""

from __future__ import annotations

import sys

from mirpathdb_prod.cli import main as cli_main


def main() -> int:
    return cli_main(["build-mirpathdb2-reference", *sys.argv[1:]])


if __name__ == "__main__":
    sys.exit(main())
