"""Reader for the frozen miRPathDB2 p-value matrices.

The v2 ``{FAMILY}_{variant}.csv`` matrices have a quirky shape (verified against
the real files): the **first line is the quoted pathway-name header with NO
leading label cell**, then every data row is ``miRNA_display_name`` followed by
one p-value per pathway. So the header has exactly one *fewer* field than each
data row, and header column ``i`` lines up with data field ``i + 1``.

Example (``hsa/GO_MF_validated_miRTarBase_strong.csv``)::

    "1-phosphatidylinositol-3-kinase activity",...,"ATP binding",...   <- header (702 cols)
    hsa-let-7a-5p,1.0,1.0,...,0.0293925,...                            <- data (703 fields)

``ATP binding`` is header index 7; ``hsa-let-7a-5p``'s value for it is data
field 8 (== 7 + 1) == ``0.0293925``.

A cell of ``1.0`` means "not tested" in v2 (no overlap / term not scored), which
the comparison treats distinctly from a genuinely computed ``p == 1.0`` only when
documented; for rank/agreement we keep the literal value.

No network, no R; just CSV parsing with :mod:`csv` (handles the quoted commas in
pathway names).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# v2's "not tested" sentinel (a literal 1.0 in the matrix).
NOT_TESTED = 1.0


@dataclass
class V2Matrix:
    """A parsed v2 p-value matrix.

    Attributes:
        pathways: Pathway display names in column order (the header).
        mirnas: miRNA display names in row order.
        values: ``float64`` array shaped ``(n_mirnas, n_pathways)`` of raw
            ORA p-values (``1.0`` == not tested in v2).
    """

    pathways: list[str]
    mirnas: list[str]
    values: np.ndarray

    def pathway_index(self) -> dict[str, int]:
        """Return ``pathway_name -> column index``."""
        return {p: i for i, p in enumerate(self.pathways)}

    def mirna_index(self) -> dict[str, int]:
        """Return ``mirna_name -> row index``."""
        return {m: i for i, m in enumerate(self.mirnas)}


def read_v2_matrix(path: str | Path) -> V2Matrix:
    """Parse a v2 ``{FAMILY}_{variant}.csv`` matrix.

    Args:
        path: Path to the CSV.

    Returns:
        A :class:`V2Matrix`. Empty (0-byte) files yield an all-empty matrix.

    Raises:
        FileNotFoundError: if ``path`` does not exist.
        ValueError: if a data row's field count does not match the header
            (mis-shaped file — fail loudly rather than mis-align cells).
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"v2 matrix not found: {p}")
    if p.stat().st_size == 0:
        return V2Matrix(pathways=[], mirnas=[], values=np.empty((0, 0), dtype=np.float64))

    with p.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            return V2Matrix(pathways=[], mirnas=[], values=np.empty((0, 0), dtype=np.float64))
        pathways = [h.strip() for h in header]
        n_path = len(pathways)

        mirnas: list[str] = []
        rows: list[list[float]] = []
        for fields in reader:
            if not fields:
                continue
            mirna = fields[0].strip()
            if not mirna:
                continue
            pvals = fields[1:]
            if len(pvals) != n_path:
                raise ValueError(
                    f"{p}: row {mirna!r} has {len(pvals)} p-values but header has "
                    f"{n_path} pathways (mis-shaped matrix)"
                )
            mirnas.append(mirna)
            rows.append([float(v) for v in pvals])

    values = (
        np.asarray(rows, dtype=np.float64)
        if rows
        else np.empty((0, n_path), dtype=np.float64)
    )
    return V2Matrix(pathways=pathways, mirnas=mirnas, values=values)


# v2 p-value-matrix directory layout under the dragon backup. The deep nested
# path is fixed by the original GeneTrail2 workspace tree.
_HSA_MATRIX_SUBPATH = (
    "extracted/hsa_raw/home/HPL/tkehl/workspace/graviton/data/"
    "mirna_pathway_dictionary/p_value_matrices/hsa"
)
_MMU_MATRIX_SUBPATH = (
    "extracted/mmu_raw/home/HPL/tkehl/workspace/graviton/data/"
    "mirna_pathway_dictionary/p_value_matrices/mmu"
)


def v2_matrix_dir(v2_root: str | Path, organism: str) -> Path:
    """Return the directory holding ``{FAMILY}_{variant}.csv`` for an organism.

    Args:
        v2_root: The ``version_2`` root (``paths.dragon_backup_v2_dir``).
        organism: ``"hsa"`` or ``"mmu"``.

    Raises:
        ValueError: for an unknown organism.
    """
    root = Path(v2_root)
    if organism == "hsa":
        return root / _HSA_MATRIX_SUBPATH
    if organism == "mmu":
        return root / _MMU_MATRIX_SUBPATH
    raise ValueError(f"unknown organism {organism!r}; expected 'hsa' or 'mmu'")


def v2_matrix_path(v2_root: str | Path, organism: str, family: str, variant: str) -> Path:
    """Return the path of one ``{family}_{variant}.csv`` v2 matrix.

    ``variant`` here is the **matrix** variant token (``validated_miRTarBase_strong``
    / ``validated_miRTarBase_all``), NOT the GeneTrail2 subtree token.
    """
    return v2_matrix_dir(v2_root, organism) / f"{family}_{variant}.csv"
