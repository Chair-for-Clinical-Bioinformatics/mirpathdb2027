#!/usr/bin/env python3
"""Build the new-production organism table.

Thin entrypoint delegating to the ``build-new-organisms`` CLI command, which
itself uses :mod:`mirpathdb_prod.normalize.organisms` and
:mod:`mirpathdb_prod.reporting.tables` to write a deterministic
``validation/organisms_new.tsv`` from the config organism list.
"""

from __future__ import annotations

import sys

from mirpathdb_prod.cli import main as cli_main


def main() -> int:
    return cli_main(["build-new-organisms", *sys.argv[1:]])


if __name__ == "__main__":
    sys.exit(main())
