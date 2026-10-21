"""Legacy/prototype source adapter (INTERNAL reference only).

Fetch/load only: reads the prototype/legacy organism baseline that was extracted
from the deprecated v3 prototype snapshot into ``validation/organisms_old.tsv``
(the canonical-schema TSV derived from the prototype ``ORGANISMS`` constant and
cross-checked against the live Postgres schemas hsa/mmu/mml).

This adapter READS reference data only; it must never import prototype Python
code. Normalization and validation happen elsewhere.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from mirpathdb_prod.sources import read_organism_tsv

# Filename of the extracted legacy baseline within the validation dir.
LEGACY_ORGANISMS_FILENAME = "organisms_old.tsv"

# Source identifier stamped on the loaded records (matches the TSV `source` col).
SOURCE = "old_prototype"


def load_organisms(tsv_path: Path) -> list[dict[str, Any]]:
    """Load raw organism records from the extracted legacy baseline TSV.

    Args:
        tsv_path: Path to ``organisms_old.tsv`` (canonical-schema TSV).

    Returns:
        Raw, un-normalized organism row dicts (hsa / mmu / mml).
    """
    return read_organism_tsv(Path(tsv_path))


def legacy_tsv_path(validation_dir: Path) -> Path:
    """Return the conventional path to the legacy baseline within a validation dir."""
    return Path(validation_dir) / LEGACY_ORGANISMS_FILENAME
