"""TargetScan target-source adapter (fetch/load, real implementation).

TargetScan v8 is a **predicted** MTI source: human-anchored, conserved-site,
broad-vertebrate. It covers human/mouse/rat/dog/chicken/macaque (and more). For
macaque (mml) it supplies the **predicted (LOW)** evidence (293 miRNAs) —
macaque's **substantive** evidence band. Macaque's only validated evidence is a
**single** DIANA-TarBase v9 MTI (``mml-miR-125b-5p -> ANGPT2``), marked into HIGH;
so macaque is a **predicted-dominant** production organism (its volume is
predicted, with the lone validated MTI flagged). miRTarBase and miRDB do NOT cover
macaque. TargetScan ships one
multi-species file keyed by *Gene Tax ID*, so a single download covers every
production organism. It feeds the **LOW** confidence band (``pred_union`` with
miRDB; their intersection is a stronger "low+"). It is also the authority for the
**conserved target SITE** conservation axis (conserved-vs-all sites; see
``docs/conservation.md``).

Per-organism coverage (organism availability matrix — verified):
    hsa ✓ anchor · mmu ✓ · rno ✓ · mml ✓ (macaque — predicted/LOW, 293 miRNAs,
    its substantive band; lone DIANA validated MTI marked)
    cfa / gga covered by TargetScan too but evaluated + EXCLUDED (owner-omitted)

The cached download is the default-predictions context-scores file (one big TSV
inside a ``.zip``); columns (verified)::

    Gene ID | Gene Symbol | Transcript ID | Gene Tax ID | miRNA | Site Type |
    UTR_start | UTR end | context++ score | context++ score percentile |
    weighted context++ score | weighted context++ score percentile |
    Predicted relative KD

``Gene Symbol`` is the symbol (v2 ORA space); ``Gene ID`` is an Ensembl id;
``miRNA`` is the miRBase-family display name. Fetch/load ONLY; no
normalization/judgment; no prototype code imported.
"""

from __future__ import annotations

import csv
import io
import zipfile
from collections.abc import Iterator
from pathlib import Path
from typing import Any

SOURCE = "targetscan"
SOURCE_VERSION = "v8"  # vert80 default predictions (see provenance/manifest.yaml)

# Predicted evidence only -> pred_union / (with miRDB) pred_intersection -> LOW.
PROVENANCE_TIERS: tuple[str, ...] = ("pred_union", "pred_intersection")

# TargetScan ships one multi-species file; these are the production organisms
# with usable (non-zero, conserved) rows. It is macaque's (mml 9544) PREDICTED
# (LOW) source — its substantive band; macaque's only validated evidence is a
# single DIANA-TarBase v9 MTI (marked).
SUPPORTED_TAXIDS: frozenset[int] = frozenset({9606, 10090, 10116, 9544})

# Column header of the default-predictions context-scores TSV (verified).
TARGETSCAN_COLUMNS: tuple[str, ...] = (
    "Gene ID",
    "Gene Symbol",
    "Transcript ID",
    "Gene Tax ID",
    "miRNA",
    "Site Type",
    "UTR_start",
    "UTR end",
    "context++ score",
    "context++ score percentile",
    "weighted context++ score",
    "weighted context++ score percentile",
    "Predicted relative KD",
)

_TAXID_COL = "Gene Tax ID"


def _open_text(path: Path) -> tuple[io.TextIOBase, zipfile.ZipFile | None]:
    """Open the predictions file as text, transparently entering a ``.zip``.

    Returns the text stream and the owning ``ZipFile`` (if any) so the caller can
    close both. For a plain ``.txt`` the second element is ``None``.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"TargetScan predictions file not found: {path}")

    if path.suffix == ".zip":
        zf = zipfile.ZipFile(path)
        members = [n for n in zf.namelist() if n.lower().endswith(".txt")]
        if not members:
            zf.close()
            raise ValueError(f"{path}: no .txt member inside the TargetScan zip")
        raw = zf.open(members[0])
        return io.TextIOWrapper(raw, encoding="utf-8", newline=""), zf

    return path.open("r", encoding="utf-8", newline=""), None


def iter_predictions(path: Path, *, taxid: int) -> Iterator[dict[str, Any]]:
    """Stream predicted MTI rows for one organism (multi-species file slice).

    Yields one dict per row whose ``Gene Tax ID`` matches ``taxid`` (keyed by
    :data:`TARGETSCAN_COLUMNS`). Streaming keeps the 137 MB file off the heap.

    Raises:
        ValueError: if ``taxid`` is unsupported or the header is unexpected.
        FileNotFoundError: if ``path`` does not exist.
    """
    if taxid not in SUPPORTED_TAXIDS:
        raise ValueError(
            f"TargetScan production slice not defined for taxid {taxid}; "
            f"supported: {sorted(SUPPORTED_TAXIDS)}"
        )
    wanted = str(taxid)
    stream, zf = _open_text(Path(path))
    try:
        reader = csv.DictReader(stream, delimiter="\t")
        header = tuple(reader.fieldnames or ())
        if header != TARGETSCAN_COLUMNS:
            raise ValueError(
                f"{path}: unexpected TargetScan header\n"
                f"  expected: {TARGETSCAN_COLUMNS}\n"
                f"  got:      {header}"
            )
        for row in reader:
            if row.get(_TAXID_COL) == wanted:
                yield dict(row)
    finally:
        stream.close()
        if zf is not None:
            zf.close()


def load_predictions(path: Path, *, taxid: int) -> list[dict[str, Any]]:
    """Load raw predicted MTIs for one organism (eager wrapper of stream).

    Fail-loud (M2): raises ``ValueError`` if the Gene-Tax-ID slice yields 0 rows
    for a requested production taxid -- TargetScan covers every production organism
    (hsa/mmu/rno/mml), so a 0-row slice means the taxid column drifted or the file
    was reshaped, not a real absence. The streaming :func:`iter_predictions` stays
    a generator (no eager guard) for callers that legitimately stream/probe.

    Args:
        path: Path to the multi-species TargetScan predictions file
            (``.zip`` or plain ``.txt``), keyed by Gene Tax ID.
        taxid: Canonical NCBI taxid to slice out
            (must be in :data:`SUPPORTED_TAXIDS`).

    Returns:
        Raw, un-normalized predicted-MTI row dicts (verbatim; no judgment).

    Raises:
        ValueError: if 0 rows match ``taxid`` (a requested production slice).
    """
    rows = list(iter_predictions(Path(path), taxid=taxid))
    if not rows:
        raise ValueError(
            f"TargetScan Gene-Tax-ID slice for taxid {taxid} yielded 0 rows from "
            f"{path}: every production organism (hsa/mmu/rno/mml) is covered, so a "
            f"0-row slice means the taxid column drifted or the file was reshaped."
        )
    return rows
