"""Frozen-release input loaders and exact production query reconstruction.

The study consumes EXACTLY the artefacts the production ORA consumed, as shipped
in the release payload (``mirpathdb-data-<date>/``):

* ``enrichment_inputs/<org>/<FAMILY>_universe.tsv``  (column ``gene``): the
  background universe of one organism x pathway-family cell (the union of that
  family's gene sets, in the family's own symbol space);
* ``enrichment_inputs/<org>/<FAMILY>_term2gene.tsv`` (``pathway_id``, ``gene``):
  the complete pathway memberships of that cell;
* ``target_pairs/<org>_target_pairs.tsv`` (``mirna``, ``target_gene``,
  ``atomic_evidence`` in ``High``/``Medium``/``Low``): the released miRNA-target
  pairs, one row per pair at its strongest evidence band, already filtered by
  the release species policy (native miRNA prefixes).

Per-miRNA queries are rebuilt the way every production runner built them
(``scripts/run_production_pathways.py``, ``run_production_enrichment.py``,
``run_macaque_inferred_reactome.py``):

1. cumulative union of the evidence bands
   (``HIGH`` = High; ``MEDIUM`` = High + Medium; ``LOW`` = High + Medium + Low),
2. :func:`mirpathdb_prod.sources.production_target_sets.normalize_to_go_symbols`
   with the casemap built from the family universe (case-normalised symbol
   matching; symbols absent from the family space are dropped),
3. miRNAs whose normalised set is empty are not scored (they are not eligible).

Nothing here is reconstructed from the FDR-filtered ``enrichment_significant``
exports: those carry only retained rows and cannot define a query or a family.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from mirpathdb_prod.sources import production_target_sets as pts

TIERS: tuple[str, ...] = ("HIGH", "MEDIUM", "LOW")

# Cumulative tier -> atomic_evidence labels it unions (the export's disjoint bands).
CUMULATIVE_ATOMIC: dict[str, tuple[str, ...]] = {
    "HIGH": ("High",),
    "MEDIUM": ("High", "Medium"),
    "LOW": ("High", "Medium", "Low"),
}

ORGANISMS: tuple[str, ...] = ("hsa", "mmu", "rno", "mml")

# Production size limits of the enricher call (minGSSize = 1, maxGSSize = 1e8):
# an eligible term is a non-empty membership within the background.
PRODUCTION_MIN_GS_SIZE = 1
PRODUCTION_MAX_GS_SIZE = 100_000_000


def sha256_file(path: Path) -> str:
    """SHA-256 hex digest of a file (streamed)."""
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _clean(value: object) -> str:
    return "" if value is None else str(value).strip()


@dataclass
class FamilySpace:
    """One organism x pathway-family background + complete memberships.

    ``matrix`` is a ``terms x genes`` CSR matrix of 0/1 memberships restricted to
    ``genes`` (the background). ``terms`` keeps only memberships that are
    non-empty within the background and inside the production size limits, i.e.
    the *eligible* terms of the fixed-family comparator; eligibility is defined
    from the background alone, before any query is seen.
    """

    organism: str
    family: str
    genes: tuple[str, ...]
    terms: tuple[str, ...]
    matrix: sp.csr_matrix
    background_label: str = "annotation"
    provenance: dict[str, object] = field(default_factory=dict)
    dropped_terms: tuple[str, ...] = ()
    _gene_index: dict[str, int] | None = field(default=None, repr=False)
    _term_sizes: np.ndarray | None = field(default=None, repr=False)

    # ---- derived -----------------------------------------------------------
    @property
    def N(self) -> int:  # noqa: N802 - statistical notation
        return len(self.genes)

    @property
    def n_terms(self) -> int:
        return len(self.terms)

    @property
    def term_sizes(self) -> np.ndarray:
        """K per eligible term (membership size within the background)."""
        if self._term_sizes is None:
            self._term_sizes = np.asarray(self.matrix.sum(axis=1)).ravel().astype(np.int64)
        return self._term_sizes

    def gene_index(self) -> dict[str, int]:
        if self._gene_index is None:
            self._gene_index = {g: i for i, g in enumerate(self.genes)}
        return self._gene_index

    def term_index(self) -> dict[str, int]:
        return {t: i for i, t in enumerate(self.terms)}

    def term_genes(self, term: str) -> list[str]:
        j = self.term_index()[term]
        row = self.matrix.getrow(j)
        return [self.genes[i] for i in row.indices]

    def term2gene(self) -> dict[str, set[str]]:
        """``{term: set(gene)}`` of the eligible memberships (for the R kernel)."""
        out: dict[str, set[str]] = {}
        indptr, indices = self.matrix.indptr, self.matrix.indices
        for j, term in enumerate(self.terms):
            out[term] = {self.genes[i] for i in indices[indptr[j] : indptr[j + 1]]}
        return out

    # ---- restriction -------------------------------------------------------
    def restrict_background(self, keep_genes: Iterable[str], label: str) -> FamilySpace:
        """Restrict background AND memberships to ``keep_genes`` (subset of genes).

        Terms whose membership becomes empty inside the new background are
        dropped from the eligible set and recorded in ``dropped_terms`` (these
        are the "terms that become unavailable" under the alternative
        background). Genes not in the current background are ignored.
        """
        index = self.gene_index()
        keep_idx = np.array(sorted({index[g] for g in keep_genes if g in index}), dtype=np.int64)
        sub = self.matrix[:, keep_idx].tocsr()
        sizes = np.asarray(sub.sum(axis=1)).ravel()
        ok = (sizes >= PRODUCTION_MIN_GS_SIZE) & (sizes <= PRODUCTION_MAX_GS_SIZE)
        dropped = tuple(t for t, flag in zip(self.terms, ok, strict=True) if not flag)
        sub = sub[np.flatnonzero(ok), :].tocsr()
        sub.eliminate_zeros()
        # Same enricher semantics as build_family_space: a kept gene that no longer
        # belongs to any eligible term is outside the effective universe.
        annotated = np.asarray(sub.sum(axis=0)).ravel() > 0
        n_unannotated = int((~annotated).sum())
        if n_unannotated:
            keep_cols = np.flatnonzero(annotated)
            sub = sub[:, keep_cols].tocsr()
            keep_idx = keep_idx[keep_cols]
        return FamilySpace(
            organism=self.organism,
            family=self.family,
            genes=tuple(self.genes[i] for i in keep_idx),
            terms=tuple(t for t, flag in zip(self.terms, ok, strict=True) if flag),
            matrix=sub,
            background_label=label,
            provenance=dict(
                self.provenance,
                restricted_from=self.background_label,
                n_restricted_unannotated_dropped=n_unannotated,
            ),
            dropped_terms=dropped,
        )


def load_universe(path: Path) -> list[str]:
    """Read a ``<FAMILY>_universe.tsv`` (header ``gene``) into a sorted unique list."""
    genes: set[str] = set()
    with Path(path).open("r", encoding="utf-8", newline="") as fh:
        header = fh.readline().strip()
        if header != "gene":
            raise ValueError(f"{path}: expected header 'gene', got {header!r}")
        for line in fh:
            g = line.strip()
            if g:
                genes.add(g)
    return sorted(genes)


def load_term2gene(path: Path) -> pd.DataFrame:
    """Read a ``<FAMILY>_term2gene.tsv`` (``pathway_id``, ``gene``) as a DataFrame."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_filter=False, quoting=3)
    if list(df.columns) != ["pathway_id", "gene"]:
        raise ValueError(f"{path}: expected columns pathway_id, gene; got {list(df.columns)}")
    df["pathway_id"] = df["pathway_id"].str.strip()
    df["gene"] = df["gene"].str.strip()
    return df[(df["pathway_id"] != "") & (df["gene"] != "")].drop_duplicates()


def build_family_space(
    organism: str,
    family: str,
    universe_path: Path,
    term2gene_path: Path,
    *,
    with_hashes: bool = True,
) -> FamilySpace:
    """Materialise a :class:`FamilySpace` from the shipped universe + term2gene."""
    genes = load_universe(universe_path)
    index = {g: i for i, g in enumerate(genes)}
    t2g = load_term2gene(term2gene_path)
    # Every term in the file is a candidate; a term whose members all fall outside
    # the universe ends up with K = 0 and is recorded as ineligible (dropped).
    terms_all = sorted(t2g["pathway_id"].unique())
    in_universe = t2g["gene"].map(index)
    n_outside = int(in_universe.isna().sum())
    t2g = t2g[in_universe.notna()]
    gene_codes = in_universe[in_universe.notna()].astype(np.int64).to_numpy()
    tindex = {t: j for j, t in enumerate(terms_all)}
    term_codes = t2g["pathway_id"].map(tindex).astype(np.int64).to_numpy()
    mat = sp.coo_matrix(
        (np.ones(len(term_codes), dtype=np.int32), (term_codes, gene_codes)),
        shape=(len(terms_all), len(genes)),
    ).tocsr()
    mat.data[:] = 1  # duplicates summed -> clip to membership
    sizes = np.asarray(mat.sum(axis=1)).ravel()
    ok = (sizes >= PRODUCTION_MIN_GS_SIZE) & (sizes <= PRODUCTION_MAX_GS_SIZE)
    dropped = tuple(t for t, flag in zip(terms_all, ok, strict=True) if not flag)
    mat = mat[np.flatnonzero(ok), :].tocsr()
    # enricher's effective universe is (genes annotated to >= 1 term) INTERSECT the
    # supplied universe (DOSE:::enricher_internal: extID <- intersect(ALLEXTID, universe)),
    # and n counts annotated query genes only. A universe gene in no eligible term
    # therefore contributes to neither N nor n; drop it here so the independent kernel
    # matches the production kernel exactly, and record how many were dropped
    # (0 in production, where the universe IS the union of the memberships).
    annotated = np.asarray(mat.sum(axis=0)).ravel() > 0
    n_unannotated = int((~annotated).sum())
    if n_unannotated:
        keep = np.flatnonzero(annotated)
        mat = mat[:, keep].tocsr()
        genes = [genes[i] for i in keep]
    prov: dict[str, object] = {
        "universe_path": str(universe_path),
        "term2gene_path": str(term2gene_path),
        "n_universe_file": int(annotated.size),
        "n_universe_unannotated_dropped": n_unannotated,
        "universe_equals_annotated_union": n_unannotated == 0,
        "n_universe": len(genes),
        "n_terms_in_file": len(terms_all),
        "n_terms_eligible": int(ok.sum()),
        "n_membership_rows": int(len(t2g)),
        "n_membership_rows_outside_universe": n_outside,
        "min_gs_size": PRODUCTION_MIN_GS_SIZE,
        "max_gs_size": PRODUCTION_MAX_GS_SIZE,
    }
    if with_hashes:
        prov["universe_sha256"] = sha256_file(universe_path)
        prov["term2gene_sha256"] = sha256_file(term2gene_path)
    return FamilySpace(
        organism=organism,
        family=family,
        genes=tuple(genes),
        terms=tuple(t for t, flag in zip(terms_all, ok, strict=True) if flag),
        matrix=mat,
        provenance=prov,
        dropped_terms=dropped,
    )


# --------------------------------------------------------------------------- #
# Target pairs -> production queries
# --------------------------------------------------------------------------- #
TARGET_PAIR_COLUMNS = ("organism", "mirna", "target_gene", "atomic_evidence")


def load_target_pairs(path: Path) -> pd.DataFrame:
    """Read ``<org>_target_pairs.tsv`` (the released atomic pairs)."""
    df = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        na_filter=False,
        quoting=3,
        usecols=list(TARGET_PAIR_COLUMNS),
    )
    missing = [c for c in TARGET_PAIR_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    bad = sorted(set(df["atomic_evidence"].unique()) - {"High", "Medium", "Low"})
    if bad:
        raise ValueError(f"{path}: unexpected atomic_evidence labels {bad}")
    return df


def raw_cumulative_sets(pairs: pd.DataFrame, tier: str) -> dict[str, set[str]]:
    """``{mirna: set(raw symbol)}`` for the cumulative tier (High..tier)."""
    if tier not in CUMULATIVE_ATOMIC:
        raise ValueError(f"unknown tier {tier!r}; expected one of {TIERS}")
    sub = pairs[pairs["atomic_evidence"].isin(CUMULATIVE_ATOMIC[tier])]
    out: dict[str, set[str]] = {}
    for mirna, genes in sub.groupby("mirna", sort=True)["target_gene"]:
        s = {g for g in genes if g}
        if s:
            out[mirna] = s
    return out


@dataclass
class Stratum:
    """One organism x family x cumulative-tier cell with its eligible queries."""

    organism: str
    family: str
    tier: str
    space: FamilySpace
    queries: dict[str, np.ndarray]  # eligible miRNA -> sorted gene indices in space
    n_raw: dict[str, int]  # raw cumulative set size per miRNA (before mapping)
    empty_mirnas: tuple[str, ...]  # had raw targets, none mapped into the family space
    normalization: dict[str, int]  # n_kept_pairs / n_dropped_pairs (production stats)

    @property
    def key(self) -> str:
        return f"{self.organism}_{self.family}_{self.tier}"

    def query_sizes(self) -> dict[str, int]:
        return {m: int(len(idx)) for m, idx in self.queries.items()}

    def query_genes(self, mirna: str) -> list[str]:
        return [self.space.genes[i] for i in self.queries[mirna]]

    def target_union(self) -> np.ndarray:
        """Sorted gene indices targeted by ANY eligible miRNA of this tier."""
        if not self.queries:
            return np.zeros(0, dtype=np.int64)
        return np.unique(np.concatenate(list(self.queries.values())))


def normalize_queries(
    raw_sets: Mapping[str, set[str]], space: FamilySpace
) -> tuple[dict[str, np.ndarray], dict[str, int], list[str]]:
    """Production normalisation of raw target sets onto ``space`` (exact reuse).

    Returns ``(queries, stats, empty_mirnas)`` where ``queries`` holds sorted
    gene indices for every miRNA with >= 1 mapped gene, ``stats`` the production
    ``n_kept_pairs``/``n_dropped_pairs`` counters, and ``empty_mirnas`` the
    miRNAs that had raw targets but none in the family space.
    """
    casemap = pts.build_symbol_casemap(space.genes)  # genes are sorted -> deterministic
    norm, kept, dropped = pts.normalize_to_go_symbols(dict(raw_sets), casemap)
    index = space.gene_index()
    queries: dict[str, np.ndarray] = {}
    for mirna in sorted(norm):
        genes = norm[mirna]
        if genes:
            queries[mirna] = np.array(sorted(index[g] for g in genes), dtype=np.int64)
    empty = sorted(set(raw_sets) - set(queries))
    return queries, {"n_kept_pairs": kept, "n_dropped_pairs": dropped}, empty


class ReleaseInputs:
    """Lazy, cached access to one release payload directory."""

    def __init__(self, inputs_dir: str | Path, *, with_hashes: bool = True) -> None:
        self.root = Path(inputs_dir).resolve()
        self.enrichment_inputs = self.root / "enrichment_inputs"
        self.target_pairs_dir = self.root / "target_pairs"
        if not self.enrichment_inputs.is_dir() or not self.target_pairs_dir.is_dir():
            raise FileNotFoundError(
                f"{self.root}: expected enrichment_inputs/ and target_pairs/ sub-directories"
            )
        self.with_hashes = with_hashes
        self._spaces: dict[tuple[str, str], FamilySpace] = {}
        self._pairs: dict[str, pd.DataFrame] = {}
        self._strata: dict[tuple[str, str, str], Stratum] = {}

    # ---- discovery ---------------------------------------------------------
    def organisms(self) -> list[str]:
        return sorted(p.name for p in self.enrichment_inputs.iterdir() if p.is_dir())

    def available_families(self, organism: str) -> list[str]:
        d = self.enrichment_inputs / organism
        fams = []
        for p in sorted(d.glob("*_universe.tsv")):
            fam = p.name[: -len("_universe.tsv")]
            if (d / f"{fam}_term2gene.tsv").exists():
                fams.append(fam)
        return fams

    def release_info(self) -> dict[str, object]:
        p = self.root / "RELEASE_INFO.json"
        if not p.exists():
            return {}
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            return {}

    # ---- loaders -----------------------------------------------------------
    def family_space(self, organism: str, family: str) -> FamilySpace:
        key = (organism, family)
        if key not in self._spaces:
            d = self.enrichment_inputs / organism
            self._spaces[key] = build_family_space(
                organism,
                family,
                d / f"{family}_universe.tsv",
                d / f"{family}_term2gene.tsv",
                with_hashes=self.with_hashes,
            )
        return self._spaces[key]

    def target_pairs(self, organism: str) -> pd.DataFrame:
        if organism not in self._pairs:
            self._pairs[organism] = load_target_pairs(
                self.target_pairs_dir / f"{organism}_target_pairs.tsv"
            )
        return self._pairs[organism]

    def target_pairs_provenance(self, organism: str) -> dict[str, object]:
        p = self.target_pairs_dir / f"{organism}_target_pairs.tsv"
        out: dict[str, object] = {"path": str(p), "bytes": p.stat().st_size}
        if self.with_hashes:
            out["sha256"] = sha256_file(p)
        return out

    def stratum(self, organism: str, family: str, tier: str) -> Stratum:
        key = (organism, family, tier)
        if key not in self._strata:
            space = self.family_space(organism, family)
            raw = raw_cumulative_sets(self.target_pairs(organism), tier)
            queries, stats, empty = normalize_queries(raw, space)
            self._strata[key] = Stratum(
                organism=organism,
                family=family,
                tier=tier,
                space=space,
                queries=queries,
                n_raw={m: len(s) for m, s in raw.items()},
                empty_mirnas=tuple(empty),
                normalization=stats,
            )
        return self._strata[key]

    def drop_cached_stratum(self, organism: str, family: str, tier: str) -> None:
        self._strata.pop((organism, family, tier), None)


# --------------------------------------------------------------------------- #
# Provenance cross-check: expected per-cell miRNA counts from the run manifests
# --------------------------------------------------------------------------- #
def load_manifest_cell_counts(provenance_dir: Path) -> dict[tuple[str, str, str], int]:
    """Collect ``n_mirnas`` per (organism, family, tier) from the production manifests.

    Sources (committed records of the production runs):
    ``production_pathways_manifest.json`` (hsa/mmu all families, rno Reactome +
    WikiPathways), ``production_enrichment_manifest.json`` (rno/mml GO) and
    ``macaque_inferred_reactome_manifest.json`` (mml inferred Reactome). The
    counts are the miRNAs the production runner SCORED, before the release
    species policy removed non-native miRNA names from the shipped contract, so
    a reconstructed cohort must be <= the manifest count.
    """
    out: dict[tuple[str, str, str], int] = {}

    def walk(obj: object, org_hint: str | None) -> None:
        if isinstance(obj, dict):
            org = (
                obj.get("organism", org_hint) if isinstance(obj.get("organism"), str) else org_hint
            )
            if "n_mirnas" in obj and "tier" in obj and "family" in obj and isinstance(org, str):
                out[(org, str(obj["family"]), str(obj["tier"]))] = int(obj["n_mirnas"])
            for v in obj.values():
                walk(v, org)
        elif isinstance(obj, list):
            for v in obj:
                walk(v, org_hint)

    for name in (
        "production_pathways_manifest.json",
        "production_enrichment_manifest.json",
        "macaque_inferred_reactome_manifest.json",
    ):
        p = Path(provenance_dir) / name
        if p.exists():
            try:
                walk(json.loads(p.read_text(encoding="utf-8")), None)
            except json.JSONDecodeError:
                continue
    return out


def percentile_sizes(sizes: Sequence[int], percentiles: Sequence[float]) -> dict[str, int]:
    """Nearest-rank percentiles of actual query sizes (always an observed size)."""
    if not sizes:
        return {}
    arr = np.asarray(sorted(int(s) for s in sizes), dtype=np.int64)
    out: dict[str, int] = {}
    for q in percentiles:
        out[f"p{int(q):02d}"] = int(np.percentile(arr, q, method="nearest"))
    return out
