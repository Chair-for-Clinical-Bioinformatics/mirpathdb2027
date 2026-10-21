"""GMT pathway-set parser + family classification.

A GMT (Gene Matrix Transposed) file is one term per line, tab-separated:

    <term_name>\t<description>\t<member1>\t<member2>\t...

This module parses that shape into a ``term -> set(members)`` mapping
(``TERM2GENE`` style) and classifies a GMT *file* into its pathway family
(``GO_BP`` / ``GO_CC`` / ``GO_MF`` / ``REACTOME`` / ``WIKIPATHWAYS`` / ``KEGG``)
from its filename, so a directory of v2 GMTs can be split family-by-family.

IMPORTANT — what the v2 ``*.gmt`` members are
---------------------------------------------
The miRPathDB2 ``*.gmt`` files ship **miRNA** members per pathway (the miRNAs
enriched for that pathway), NOT gene symbols. They are therefore the wrong
source for a ``TERM2GENE`` (pathway -> gene) mapping used to drive ORA. The
gene-level ``TERM2GENE`` for the controlled v2 reproduction is reconstructed
from the GeneTrail2 result trees instead (see
:mod:`mirpathdb_prod.pathways.v2_genesets`). This parser is the generic GMT
reader (it does not assume members are genes vs miRNAs); callers label the
member space.

No network, no R, pure parsing.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

# Canonical pathway families (the comparison/build axis). KEGG is recognised so
# v2 KEGG GMTs are not silently mis-bucketed, but it is EXCLUDED from the
# published DB / controlled comparison by license (see config/pathway_sources.yaml).
FAMILY_GO_BP = "GO_BP"
FAMILY_GO_CC = "GO_CC"
FAMILY_GO_MF = "GO_MF"
FAMILY_REACTOME = "REACTOME"
FAMILY_WIKIPATHWAYS = "WIKIPATHWAYS"
FAMILY_KEGG = "KEGG"

#: All recognised families (longest tokens first so "GO_BP" wins over a bare
#: "GO" and "WIKIPATHWAYS" is matched whole).
KNOWN_FAMILIES: tuple[str, ...] = (
    FAMILY_GO_BP,
    FAMILY_GO_CC,
    FAMILY_GO_MF,
    FAMILY_REACTOME,
    FAMILY_WIKIPATHWAYS,
    FAMILY_KEGG,
)

#: Families that ship in the published DB / take part in the controlled v2
#: comparison (KEGG excluded by license).
PUBLISHABLE_FAMILIES: tuple[str, ...] = (
    FAMILY_GO_BP,
    FAMILY_GO_CC,
    FAMILY_GO_MF,
    FAMILY_REACTOME,
    FAMILY_WIKIPATHWAYS,
)


@dataclass(frozen=True)
class GmtTerm:
    """One parsed GMT line.

    Attributes:
        name: The term name (column 1).
        description: The term description (column 2; often a URL in v2 GMTs).
        members: The member identifiers (columns 3..N) as a frozenset
            (duplicates collapsed, blank fields dropped).
    """

    name: str
    description: str
    members: frozenset[str]


def classify_family(path: str | Path) -> str:
    """Classify a GMT file into its pathway family from its filename.

    v2 GMTs are named ``<FAMILY>_<variant>.gmt`` (e.g.
    ``GO_MF_validated_miRTarBase_strong.gmt``). The family is matched
    case-sensitively against :data:`KNOWN_FAMILIES`, longest token first so
    ``GO_BP`` is not shortened to ``GO`` and ``WIKIPATHWAYS`` matches whole.

    Args:
        path: GMT file path (only its basename is inspected).

    Returns:
        One of :data:`KNOWN_FAMILIES`.

    Raises:
        ValueError: if no known family prefixes the filename.
    """
    stem = Path(path).name
    for family in KNOWN_FAMILIES:
        # Family must be a prefix followed by a separator (``_`` or ``.``) so
        # "GO_BPX..." cannot match and "GO_BP_..." / "GO_BP.gmt" both do.
        if stem.startswith(family) and (
            len(stem) == len(family) or stem[len(family)] in ("_", ".")
        ):
            return family
    raise ValueError(
        f"cannot classify GMT {stem!r} into a known family; known: {KNOWN_FAMILIES}"
    )


def iter_gmt(path: str | Path) -> Iterator[GmtTerm]:
    """Yield :class:`GmtTerm` records from a GMT file (streaming).

    Blank lines and lines without at least a name + description are skipped.
    Member fields that are empty (e.g. a trailing tab) are dropped.

    Args:
        path: GMT file path.

    Yields:
        One :class:`GmtTerm` per non-empty term line.

    Raises:
        FileNotFoundError: if ``path`` does not exist.
    """
    gmt_path = Path(path)
    if not gmt_path.exists():
        raise FileNotFoundError(f"GMT file not found: {gmt_path}")
    with gmt_path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                # Not a well-formed GMT line (need at least name + description).
                continue
            name = cols[0].strip()
            description = cols[1].strip()
            members = frozenset(c.strip() for c in cols[2:] if c.strip())
            if not name:
                continue
            yield GmtTerm(name=name, description=description, members=members)


def parse_gmt(path: str | Path) -> dict[str, set[str]]:
    """Parse a GMT file into a ``term -> set(members)`` mapping (TERM2GENE-style).

    Args:
        path: GMT file path.

    Returns:
        Mapping from term name to the set of its members. If a term name
        repeats, the member sets are unioned.

    Raises:
        FileNotFoundError: if ``path`` does not exist.
    """
    out: dict[str, set[str]] = {}
    for term in iter_gmt(path):
        out.setdefault(term.name, set()).update(term.members)
    return out


def parse_gmt_dir(directory: str | Path) -> dict[str, dict[str, set[str]]]:
    """Parse every ``*.gmt`` under ``directory``, split by pathway family.

    Files that cannot be classified into a known family are skipped (not an
    error) so a directory may carry unrelated files. Within a family, multiple
    GMTs (e.g. different variants) are merged term-by-term.

    Args:
        directory: Directory holding ``*.gmt`` files.

    Returns:
        Mapping ``family -> {term -> set(members)}`` for each family present.

    Raises:
        FileNotFoundError: if ``directory`` does not exist or is not a directory.
    """
    root = Path(directory)
    if not root.is_dir():
        raise FileNotFoundError(f"GMT directory not found: {root}")
    by_family: dict[str, dict[str, set[str]]] = {}
    for gmt_file in sorted(root.glob("*.gmt")):
        try:
            family = classify_family(gmt_file)
        except ValueError:
            continue
        family_map = by_family.setdefault(family, {})
        for term in iter_gmt(gmt_file):
            family_map.setdefault(term.name, set()).update(term.members)
    return by_family


def term2gene_rows(term2gene: dict[str, set[str]]) -> list[tuple[str, str]]:
    """Flatten a ``term -> set(genes)`` mapping into ``(term, gene)`` row pairs.

    This is the two-column long form clusterProfiler's ``enricher`` consumes as
    its ``TERM2GENE`` data frame. Rows are emitted in a deterministic order
    (terms sorted, genes sorted within a term) so a written temp file is stable.

    Args:
        term2gene: Mapping from term name to its gene set.

    Returns:
        List of ``(term, gene)`` tuples.
    """
    rows: list[tuple[str, str]] = []
    for term in sorted(term2gene):
        for gene in sorted(term2gene[term]):
            rows.append((term, gene))
    return rows
