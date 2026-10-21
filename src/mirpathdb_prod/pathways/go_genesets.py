"""Production GO gene-set source (``org.*.eg.db`` -> TERM2GENE), Bioconductor bridge.

This is the **production** GO gene-set source for the rat/macaque (and any
``org.*.eg.db``-backed) enrichment build. Unlike the controlled v2 reproduction
-- which reconstructs pathway -> gene membership from frozen GeneTrail2 outputs --
production GO gene sets come straight from the Bioconductor organism annotation
package (``org.Rn.eg.db`` for rat, ``org.Mmu.eg.db`` for macaque, ``org.Hs/Mm``
for human/mouse). For each organism we materialise, **split by GO ontology**
(GO_BP / GO_CC / GO_MF), a ``TERM2GENE`` map ``GO_id -> set(gene SYMBOL)`` plus a
``TERM2NAME`` map ``GO_id -> GO term name`` for clusterProfiler ``enricher``.

ID space (the #1 silent-failure risk): the TERM2GENE is built in **gene-SYMBOL**
space (``org.*.eg.db`` ``SYMBOL`` column). Every target-set builder
(:mod:`mirpathdb_prod.sources.production_target_sets`) emits the SAME symbol space
(DIANA ``gene_name`` / miRTarBase ``Target Gene`` / TargetScan ``Gene Symbol`` are
already symbols; miRDB RefSeq accessions are mapped to symbol **through the same
``org.*.eg.db``** so the spaces cannot drift). Because clusterProfiler ``enricher``
is ID-agnostic and the TERM2GENE drives the universe, symbol-vs-symbol is the only
contract that matters.

Reproducibility (REQUIRED): the GO annotation in ``org.*.eg.db`` drifts between
Bioconductor releases, so a real build MUST pin the ``org_db_package`` version and
the Bioconductor/GO release. The materialiser records the package version it read
alongside the TERM2GENE so the produced gene sets are auditable.

This module shells out to R (like :mod:`mirpathdb_prod.enrichment.clusterprofiler`);
it does not import ``rpy2``. The R binary is ``$MPD_RSCRIPT`` if set, else the
first ``Rscript`` on ``PATH``.
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

# Canonical taxid -> Bioconductor organism annotation package (the GO source).
ORG_DB_PACKAGES: dict[int, str] = {
    9606: "org.Hs.eg.db",
    10090: "org.Mm.eg.db",
    10116: "org.Rn.eg.db",
    9544: "org.Mmu.eg.db",
}

# The GO families (ontologies) we split into. These are the column tokens used in
# data/processed/enrichment/<org>/ filenames and in the build manifest.
GO_FAMILIES: tuple[str, ...] = ("GO_BP", "GO_CC", "GO_MF")

# Ontology code (org.db ONTOLOGY column) -> family token.
_ONTOLOGY_TO_FAMILY: dict[str, str] = {"BP": "GO_BP", "CC": "GO_CC", "MF": "GO_MF"}

# Environment override for the R interpreter (mirrors the clusterProfiler bridge).
RSCRIPT_ENV = "MPD_RSCRIPT"

# GO annotation propagation mode -> (gene->term column, ontology column) in
# org.*.eg.db via AnnotationDbi::select. ``all`` uses the PROPAGATED columns
# (GOALL / ONTOLOGYALL): a gene annotated to a child term is also counted toward
# every ancestor term (the GO true-path rule), matching canonical GO ORA
# (clusterProfiler::enrichGO, GeneTrail2, topGO). ``direct`` uses the DIRECT
# columns (GO / ONTOLOGY): only the gene's own annotations, which systematically
# under-populates higher-level terms. Default is ``all`` (canonical); overridable
# via ``config.yaml``'s ``go_propagation:`` key.
GO_PROPAGATION_COLUMNS: dict[str, tuple[str, str]] = {
    "all": ("GOALL", "ONTOLOGYALL"),
    "direct": ("GO", "ONTOLOGY"),
}
DEFAULT_GO_PROPAGATION = "all"


def resolve_go_propagation(raw_config: dict, *, default: str = DEFAULT_GO_PROPAGATION) -> str:
    """Read + validate the ``go_propagation`` mode from a raw config dict.

    Returns ``"all"`` (canonical, propagated GOALL) unless ``config.yaml`` sets
    ``go_propagation: direct``. Raises ``ValueError`` on an unknown mode so a typo
    fails loud rather than silently falling back.
    """
    mode = (raw_config or {}).get("go_propagation", default)
    if mode not in GO_PROPAGATION_COLUMNS:
        raise ValueError(
            f"unknown go_propagation {mode!r}; choose from "
            f"{sorted(GO_PROPAGATION_COLUMNS)}"
        )
    return mode

# The R driver. Loads the organism annotation package, pulls every
# ENTREZID -> (GO, ONTOLOGY, SYMBOL) row, drops missing values, and writes a
# long-form TSV (family<TAB>go_id<TAB>go_name<TAB>symbol) plus a provenance JSON
# (org.db package + version, GO.db version, R version). GO term NAMES come from
# GO.db (a dependency of org.*.eg.db). Kept inline so there is no extra file to
# ship; written to a temp dir per call.
_R_GO_DRIVER = r"""
suppressMessages({
  library(AnnotationDbi)
  library(GO.db)
})
args    <- commandArgs(trailingOnly = TRUE)
pkg     <- args[[1]]
out_file  <- args[[2]]
prov_file <- args[[3]]
go_col  <- args[[4]]   # "GOALL" (propagated) or "GO" (direct)
ont_col <- args[[5]]   # "ONTOLOGYALL" or "ONTOLOGY" (matches go_col)

suppressMessages(library(pkg, character.only = TRUE))
db <- get(pkg)

# GOALL/ONTOLOGYALL = the propagated (true-path) mapping so a gene counts toward
# every ancestor term; GO/ONTOLOGY = direct annotations only. Column indexed by
# name so the same driver serves both modes. (A GOALL select also returns an
# EVIDENCEALL column, discarded by the projection below.)
go <- suppressMessages(AnnotationDbi::select(
  db,
  keys     = keys(db, keytype = "ENTREZID"),
  columns  = c(go_col, ont_col, "SYMBOL"),
  keytype  = "ENTREZID"
))
go <- go[!is.na(go[[go_col]]) & !is.na(go[[ont_col]]) & !is.na(go$SYMBOL), ]
go <- unique(go[, c(go_col, ont_col, "SYMBOL")])

# GO term names from GO.db (one lookup over the distinct GO ids present).
ids   <- unique(go[[go_col]])
terms <- suppressMessages(AnnotationDbi::select(
  GO.db, keys = ids, columns = "TERM", keytype = "GOID"))
name_of <- setNames(terms$TERM, terms$GOID)

fam <- c(BP = "GO_BP", CC = "GO_CC", MF = "GO_MF")[go[[ont_col]]]
out <- data.frame(
  family = fam,
  go_id  = go[[go_col]],
  go_name = ifelse(is.na(name_of[go[[go_col]]]), go[[go_col]], name_of[go[[go_col]]]),
  symbol = go$SYMBOL,
  stringsAsFactors = FALSE
)
out <- out[!is.na(out$family), ]
write.table(out, file = out_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

manual_json <- function(x) {
  parts <- vapply(names(x), function(k)
    sprintf('"%s": "%s"', k, gsub('"', '\\"', as.character(x[[k]]))), character(1))
  paste0("{", paste(parts, collapse = ", "), "}")
}
prov <- list(
  org_db_package = pkg,
  org_db_version = as.character(packageVersion(pkg)),
  go_db_version  = as.character(packageVersion("GO.db")),
  r_version      = paste(R.version$major, R.version$minor, sep = "."),
  go_columns     = paste(go_col, ont_col, sep = "/"),
  n_pairs        = as.character(nrow(out))
)
writeLines(manual_json(prov), prov_file)
"""


@dataclass
class GoGeneSets:
    """Materialised GO gene sets for one organism, split by ontology.

    Attributes:
        taxid: Canonical NCBI taxid the gene sets were built for.
        org_db_package: The ``org.*.eg.db`` package the GO annotation came from.
        term2gene: ``{family: {go_id: set(symbol)}}`` for GO_BP / GO_CC / GO_MF.
        term2name: ``{go_id: go_term_name}`` (shared across families).
        provenance: Version pins (org.db package + version, GO.db version, R).
    """

    taxid: int
    org_db_package: str
    term2gene: dict[str, dict[str, set[str]]]
    term2name: dict[str, str]
    provenance: dict[str, str] = field(default_factory=dict)

    def family(self, family: str) -> dict[str, set[str]]:
        """Return the ``go_id -> set(symbol)`` map for one GO family."""
        if family not in self.term2gene:
            raise KeyError(f"no GO family {family!r}; have {sorted(self.term2gene)}")
        return self.term2gene[family]

    def universe(self, family: str) -> set[str]:
        """Return the background = union of all symbols in a family's terms."""
        genes: set[str] = set()
        for s in self.term2gene.get(family, {}).values():
            genes |= s
        return genes

    def summary(self) -> dict[str, dict[str, int]]:
        """Per-family ``{n_terms, n_genes, n_pairs}`` count summary (for manifests)."""
        out: dict[str, dict[str, int]] = {}
        for fam, t2g in self.term2gene.items():
            genes: set[str] = set()
            pairs = 0
            for g in t2g.values():
                genes |= g
                pairs += len(g)
            out[fam] = {"n_terms": len(t2g), "n_genes": len(genes), "n_pairs": pairs}
        return out


def _resolve_rscript(rscript: str | None = None) -> str:
    """Return the Rscript binary ($MPD_RSCRIPT, explicit arg, or PATH lookup)."""
    if rscript:
        return rscript
    override = os.environ.get(RSCRIPT_ENV)
    if override:
        return override
    found = shutil.which("Rscript")
    if found:
        return found
    raise RuntimeError(
        "Rscript not found: set $MPD_RSCRIPT to the clusterProfiler conda env's "
        "Rscript (e.g. .../envs/mpd_cp/bin/Rscript) or put it on PATH"
    )


def org_db_for_taxid(taxid: int) -> str:
    """Return the ``org.*.eg.db`` package name for a canonical taxid."""
    try:
        return ORG_DB_PACKAGES[taxid]
    except KeyError as exc:
        raise ValueError(
            f"no org.*.eg.db package for taxid {taxid}; known: {sorted(ORG_DB_PACKAGES)}"
        ) from exc


def _parse_go_tsv(path: Path) -> tuple[dict[str, dict[str, set[str]]], dict[str, str]]:
    """Parse the R-written long-form GO TSV into (term2gene, term2name)."""
    term2gene: dict[str, dict[str, set[str]]] = {f: {} for f in GO_FAMILIES}
    term2name: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            family = (row.get("family") or "").strip()
            go_id = (row.get("go_id") or "").strip()
            symbol = (row.get("symbol") or "").strip()
            if family not in term2gene or not go_id or not symbol:
                continue
            term2gene[family].setdefault(go_id, set()).add(symbol)
            name = (row.get("go_name") or "").strip()
            if name:
                term2name.setdefault(go_id, name)
    return term2gene, term2name


def build_go_genesets(
    taxid: int, *, rscript: str | None = None, propagation: str = DEFAULT_GO_PROPAGATION
) -> GoGeneSets:
    """Materialise GO gene sets (TERM2GENE/TERM2NAME) for ``taxid`` from org.*.eg.db.

    Shells out to R once, reads ``ENTREZID -> (go_col, ont_col, SYMBOL)`` from the
    organism's annotation package, names terms from GO.db, and returns the
    symbol-space gene sets split by GO_BP / GO_CC / GO_MF.

    Args:
        taxid: Canonical NCBI taxid (must be in :data:`ORG_DB_PACKAGES`).
        rscript: Override the Rscript binary (else $MPD_RSCRIPT / PATH).
        propagation: ``"all"`` (default; propagated GOALL/ONTOLOGYALL = canonical
            GO ORA, ancestors included) or ``"direct"`` (GO/ONTOLOGY, direct
            annotations only). See :data:`GO_PROPAGATION_COLUMNS`.

    Returns:
        A :class:`GoGeneSets` with the per-family TERM2GENE, shared TERM2NAME, and
        the org.db/GO.db/R version pins + the ``go_propagation`` mode in
        ``provenance``.

    Raises:
        ValueError: if ``taxid`` has no org.*.eg.db mapping, or ``propagation`` is
            unknown.
        RuntimeError: if the R subprocess fails (stderr is surfaced).
    """
    if propagation not in GO_PROPAGATION_COLUMNS:
        raise ValueError(
            f"unknown propagation {propagation!r}; choose from "
            f"{sorted(GO_PROPAGATION_COLUMNS)}"
        )
    go_col, ont_col = GO_PROPAGATION_COLUMNS[propagation]
    pkg = org_db_for_taxid(taxid)
    rscript_bin = _resolve_rscript(rscript)

    with tempfile.TemporaryDirectory(prefix="mpd_go_t2g_") as tmp:
        tmp_path = Path(tmp)
        out_file = tmp_path / "go_term2gene.tsv"
        prov_file = tmp_path / "provenance.json"
        driver_file = tmp_path / "driver_go.R"
        driver_file.write_text(_R_GO_DRIVER, encoding="utf-8")

        proc = subprocess.run(
            [
                rscript_bin,
                str(driver_file),
                pkg,
                str(out_file),
                str(prov_file),
                go_col,
                ont_col,
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"GO TERM2GENE Rscript failed for {pkg} (exit {proc.returncode}).\n"
                f"--- stderr ---\n{proc.stderr}\n--- stdout ---\n{proc.stdout}"
            )

        term2gene, term2name = _parse_go_tsv(out_file)
        provenance: dict[str, str] = {}
        if prov_file.exists():
            try:
                provenance = json.loads(prov_file.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                provenance = {}
    # Stamp the propagation mode so the manifest self-documents the GO semantics.
    provenance["go_propagation"] = propagation

    return GoGeneSets(
        taxid=taxid,
        org_db_package=pkg,
        term2gene=term2gene,
        term2name=term2name,
        provenance=provenance,
    )


def write_go_genesets_tsv(gene_sets: GoGeneSets, out_path: str | Path) -> Path:
    """Persist a materialised GO TERM2GENE as a long-form TSV (cache/audit).

    Columns: ``family<TAB>go_id<TAB>go_name<TAB>symbol`` (one row per pair),
    deterministic order. The companion :func:`read_go_genesets_tsv` reloads it so
    a build can reuse a cached materialisation without re-shelling to R.
    """
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", quoting=csv.QUOTE_NONE, escapechar="\\")
        writer.writerow(["family", "go_id", "go_name", "symbol"])
        for family in GO_FAMILIES:
            t2g = gene_sets.term2gene.get(family, {})
            for go_id in sorted(t2g):
                name = gene_sets.term2name.get(go_id, go_id)
                for symbol in sorted(t2g[go_id]):
                    writer.writerow([family, go_id, name, symbol])
    return out


def read_go_genesets_tsv(path: str | Path, *, taxid: int, org_db_package: str) -> GoGeneSets:
    """Reload a cached GO TERM2GENE TSV written by :func:`write_go_genesets_tsv`."""
    term2gene, term2name = _parse_go_tsv(Path(path))
    return GoGeneSets(
        taxid=taxid,
        org_db_package=org_db_package,
        term2gene=term2gene,
        term2name=term2name,
    )
