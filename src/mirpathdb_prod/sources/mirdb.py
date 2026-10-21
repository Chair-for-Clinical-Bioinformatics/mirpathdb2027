"""miRDB target-source adapter (fetch/load, real implementation).

miRDB v6.0 is a **predicted** MTI source and one of the FOUR official target
sources. Its **MirTarget** prediction model is fit on miRNA-overexpression
studies and is heavily used in the field — its rationale for inclusion. It is a
sanctioned, free/redistributable predictor that feeds the **LOW** confidence
band (``pred_union`` with TargetScan; their intersection is a stronger "low+";
see ``docs/target_evidence.md``). It is NOT an organism-selection anchor: miRDB
v6.0 ships exactly 5 species (hsa, mmu, rat, dog, chicken — NO macaque); the owner
EXCLUDED dog + chicken (2026-06-24), so among the production organisms miRDB covers
hsa/mmu/rno (NOT macaque).

Per-organism coverage (organism availability matrix — verified):
    hsa ✓ · mmu ✓ · rno ✓ (production)
    cfa / gga shipped by miRDB but EVALUATED + EXCLUDED (owner-omitted)
    mml ✗ (macaque NOT among miRDB's 5 species; its predicted source is TargetScan
    — its substantive band; its only validated evidence is a single DIANA-TarBase v9
    MTI, marked)

The cached download is the all-species prediction-result table (gzip TSV, no
header) with three columns::

    miRNA (species-prefixed display name) | target (RefSeq mRNA accession) |
    target_score (MirTarget score, 50–100)

There is no per-species file: rows are sliced by the miRNA name's species
prefix (``hsa-``/``mmu-``/``rno-``). The target is a RefSeq accession (NOT a
symbol); mapping accession->symbol is a normalization concern handled elsewhere.
Fetch/load ONLY; no normalization/judgment; no prototype code imported.
"""

from __future__ import annotations

import csv
import gzip
from collections.abc import Iterator
from pathlib import Path
from typing import Any

SOURCE = "mirdb"
SOURCE_VERSION = "v6.0"

# Predicted evidence only -> pred_union / (with TargetScan) pred_intersection -> LOW.
PROVENANCE_TIERS: tuple[str, ...] = ("pred_union", "pred_intersection")

# Production organisms with native miRDB coverage (canonical taxids). miRDB
# v6.0 also ships cfa (9615) + gga (9031), but those were evaluated + EXCLUDED
# from the production set; mml (9544) is not shipped by miRDB.
SUPPORTED_TAXIDS: frozenset[int] = frozenset({9606, 10090, 10116})

# The all-species miRDB table has no header; these are the canonical column
# names this adapter assigns positionally.
MIRDB_COLUMNS: tuple[str, ...] = ("mirna", "target", "target_score")

# Canonical taxid -> miRNA-name species prefix used to slice the all-species
# table (miRDB names are miRBase-style ``<prefix>-...``).
_TAXID_TO_PREFIX: dict[int, str] = {
    9606: "hsa-",
    10090: "mmu-",
    10116: "rno-",
}


def prefix_for_taxid(taxid: int) -> str:
    """Return the miRNA-name species prefix used to slice miRDB for an organism."""
    try:
        return _TAXID_TO_PREFIX[taxid]
    except KeyError as exc:
        raise ValueError(
            f"no miRDB species prefix for taxid {taxid}; "
            f"supported: {sorted(SUPPORTED_TAXIDS)}"
        ) from exc


def iter_predictions(path: Path, *, taxid: int) -> Iterator[dict[str, Any]]:
    """Stream predicted MTI rows for one organism (slice the all-species table).

    Yields one dict per row whose miRNA name starts with the organism's species
    prefix (keyed by :data:`MIRDB_COLUMNS`). Streaming keeps the large gzip TSV
    off the heap.

    Raises:
        ValueError: if ``taxid`` is unsupported.
        FileNotFoundError: if ``path`` does not exist.
    """
    if taxid not in SUPPORTED_TAXIDS:
        raise ValueError(
            f"miRDB has no production slice for taxid {taxid}; "
            f"supported: {sorted(SUPPORTED_TAXIDS)}"
        )
    prefix = prefix_for_taxid(taxid)

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"miRDB prediction file not found: {path}")

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or not row[0].startswith(prefix):
                continue
            # pad/truncate defensively to the 3 canonical columns
            values = (row + ["", "", ""])[: len(MIRDB_COLUMNS)]
            yield dict(zip(MIRDB_COLUMNS, values, strict=True))


def load_predictions(path: Path, *, taxid: int) -> list[dict[str, Any]]:
    """Load raw predicted MTIs for one organism (eager wrapper of stream).

    Fail-loud (M2): raises ``ValueError`` if the species-prefix slice yields 0
    rows for a requested production taxid -- miRDB covers every production organism
    it supports (hsa/mmu/rno), so a 0-row slice means the prefix drifted or the
    file was reshaped, not a real absence. The streaming :func:`iter_predictions`
    stays a generator (no eager guard) for callers that legitimately stream/probe.

    Args:
        path: Path to the miRDB v6.0 all-species prediction file (gzip TSV).
        taxid: Canonical NCBI taxid to load (must be in
            :data:`SUPPORTED_TAXIDS`).

    Returns:
        Raw, un-normalized predicted-MTI row dicts (verbatim; no judgment).
        ``target`` is a RefSeq accession (mapping to symbol happens in
        normalization).

    Raises:
        ValueError: if 0 rows match ``taxid``'s species prefix (a requested slice).
    """
    rows = list(iter_predictions(Path(path), taxid=taxid))
    if not rows:
        raise ValueError(
            f"miRDB species-prefix slice for taxid {taxid} "
            f"({prefix_for_taxid(taxid)!r}) yielded 0 rows from {path}: this is a "
            f"requested production slice, so a 0-row result means the prefix "
            f"drifted or the file was reshaped, not a real absence."
        )
    return rows
