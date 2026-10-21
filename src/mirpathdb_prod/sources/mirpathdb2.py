"""miRPathDB2 source adapter (PRIMARY external validation source).

Fetch/load only: reads the miRPathDB2 organism reference that was derived from
the offline ``version_2/extracted/`` p-value matrices (hsa_raw/ and mmu_raw/)
into ``validation/organisms_mirpathdb2.tsv`` (canonical-schema TSV).

This is the dedicated adapter the briefing requires; it keeps the miRPathDB2
reference explicit and separate from the legacy/prototype loader.

miRPathDB2 (Backes et al. 2020, NAR, doi:10.1093/nar/gkz1022) covers
Human + Mouse only — there is NO rhesus (mml) in this source.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from mirpathdb_prod.sources import read_organism_tsv

# Filename of the miRPathDB2 reference within the validation dir.
MIRPATHDB2_ORGANISMS_FILENAME = "organisms_mirpathdb2.tsv"

# Source identifier stamped on the loaded records (matches the TSV `source` col).
SOURCE = "mirpathdb2"

# miRPathDB2 is Human + Mouse only — this is an invariant of the reference and is
# asserted at load time so a reshaped/extended file is caught immediately.
EXPECTED_TAXIDS: frozenset[int] = frozenset({9606, 10090})


def load_organisms(tsv_path: Path) -> list[dict[str, Any]]:
    """Load raw organism records from the miRPathDB2 reference TSV.

    Args:
        tsv_path: Path to ``organisms_mirpathdb2.tsv`` (canonical-schema TSV).

    Returns:
        Raw, un-normalized organism row dicts (hsa / mmu only).

    Raises:
        ValueError: if the reference contains a taxid outside the documented
            miRPathDB2 set (Human + Mouse), which would indicate a corrupted or
            wrong reference file.
    """
    rows = read_organism_tsv(Path(tsv_path))
    taxids = {int(r["canonical_tax_id"]) for r in rows if str(r["canonical_tax_id"]).strip()}
    unexpected = taxids - EXPECTED_TAXIDS
    if unexpected:
        raise ValueError(
            f"{tsv_path}: miRPathDB2 reference contains unexpected taxids {sorted(unexpected)}; "
            f"miRPathDB2 covers Human (9606) + Mouse (10090) only"
        )
    return rows


def mirpathdb2_tsv_path(validation_dir: Path) -> Path:
    """Return the conventional path to the miRPathDB2 reference within a validation dir."""
    return Path(validation_dir) / MIRPATHDB2_ORGANISMS_FILENAME
