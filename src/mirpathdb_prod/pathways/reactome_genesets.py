"""Production Reactome gene-set source (``reactome.db`` x ``org.*.eg.db`` -> TERM2GENE).

This is the **production** Reactome gene-set source for the uniform GO + Reactome
+ WikiPathways enrichment build. It mirrors :mod:`mirpathdb_prod.pathways.go_genesets`:
for each organism we materialise a single-family ``TERM2GENE`` map
``reactome_pathway_id -> set(gene SYMBOL)`` plus a ``TERM2NAME`` map
``reactome_pathway_id -> pathway name`` for clusterProfiler ``enricher``.

Source + ID space
-----------------
The gene sets come from the Bioconductor **``reactome.db``** annotation package
(``reactomePATHID2EXTID`` = Reactome pathway -> Entrez gene id;
``reactomePATHID2NAME`` = pathway id -> ``"<Species>: <name>"``), joined to the
organism's **``org.*.eg.db``** to turn Entrez ids into gene SYMBOLs. Using the
SAME annotation package the GO sets came from is what keeps Reactome in the GO /
target-set gene-SYMBOL space (the #1 silent-failure risk). Entrez ids with no
symbol are dropped; pathways left empty after the join are dropped.

Per-organism split
-------------------
Reactome encodes the organism in the pathway id stable prefix
``R-<SPECIES>-<number>`` (``R-HSA-`` human, ``R-MMU-`` mouse, ``R-RNO-`` rat).
We slice ``PATHID2EXTID`` by that prefix so each organism gets ONLY its own
pathways.

Organism coverage (the #1 thing to *check*, not assume)
-------------------------------------------------------
Reactome ships **NO curated Macaca mulatta** (``R-MML-``) pathways -- 0 in
``reactome.db`` v1.95.0 -- so **macaque gets no native / curated Reactome** here
and :data:`REACTOME_SPECIES_CODE` has no 9544 entry; :func:`species_code_for_taxid`
raises for it. This map is the **native / curated** source only. Reactome *does*
list ``Macaca mulatta`` as an **orthology-inferred** species; that INFERRED macaque
Reactome (human ``R-HSA-`` projected to ``R-MML-`` via NCBI orthology) is built
separately by :mod:`mirpathdb_prod.pathways.inferred_reactome_genesets` and is NOT
part of this native map. Native coverage in v1.95.0: HSA 2848 / MMU 1822 / RNO 1809
pathways.

This module shells out to R (like the GO source); it does not import ``rpy2``.
The R binary is ``$MPD_RSCRIPT`` if set, else the first ``Rscript`` on ``PATH``.
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

from mirpathdb_prod.pathways.go_genesets import ORG_DB_PACKAGES, org_db_for_taxid

# The single pathway family this source emits (matches the gmt.py token).
REACTOME_FAMILY = "REACTOME"

# Canonical taxid -> Reactome stable-id species code (the ``R-<CODE>-`` prefix).
# ONLY species Reactome ships NATIVE / CURATED pathways for are listed; macaque
# (9544) is deliberately ABSENT -- reactome.db has 0 curated ``R-MML-`` pathways,
# so macaque gets no native / curated Reactome. (Inferred macaque Reactome --
# human ``R-HSA-`` orthology-projected to ``R-MML-`` -- is built separately by
# ``inferred_reactome_genesets`` and is intentionally NOT in this native map.)
REACTOME_SPECIES_CODE: dict[int, str] = {
    9606: "HSA",
    10090: "MMU",
    10116: "RNO",
}

RSCRIPT_ENV = "MPD_RSCRIPT"

# R driver: for one organism (species code + org.db package), pull every Reactome
# pathway whose stable id starts with R-<CODE>-, expand it to Entrez members via
# reactomePATHID2EXTID, join Entrez->SYMBOL through org.*.eg.db, and write a
# long-form TSV (pathway_id<TAB>pathway_name<TAB>symbol) + a provenance JSON. The
# pathway NAME comes from reactomePATHID2NAME with the leading "<Species>: " prefix
# stripped so names match the human-readable Reactome label.
_R_REACTOME_DRIVER = r"""
suppressMessages({
  library(AnnotationDbi)
  library(reactome.db)
})
args      <- commandArgs(trailingOnly = TRUE)
species   <- args[[1]]              # e.g. "HSA"
pkg       <- args[[2]]              # e.g. "org.Hs.eg.db"
out_file  <- args[[3]]
prov_file <- args[[4]]

suppressMessages(library(pkg, character.only = TRUE))
db <- get(pkg)

# All Reactome pathway ids for this organism (R-<CODE>-...) that actually carry
# Entrez members. Not every PATHID has an EXTID mapping (container/event nodes
# have none), so restrict to mappedkeys() before subsetting or as.list() errors.
prefix    <- paste0("R-", species, "-")
mapped    <- mappedkeys(reactomePATHID2EXTID)
ids       <- mapped[startsWith(mapped, prefix)]

# pathway -> Entrez members (list keyed by pathway id).
p2e <- as.list(reactomePATHID2EXTID[ids])
p2e <- p2e[!vapply(p2e, is.null, logical(1))]

# pathway -> name; strip the leading "<Species>: " label. Restrict to ids that
# actually have a name mapping (mappedkeys) so the subset cannot error.
named_ids <- intersect(names(p2e), mappedkeys(reactomePATHID2NAME))
nm_list   <- if (length(named_ids)) as.list(reactomePATHID2NAME[named_ids]) else list()
name_of <- vapply(names(p2e), function(id) {
  v <- nm_list[[id]]
  if (is.null(v) || length(v) == 0 || is.na(v[[1]])) return(id)
  sub("^[^:]+:\\s*", "", v[[1]])
}, character(1))

# Build the full Entrez universe across this organism's pathways, then one
# org.db join Entrez -> SYMBOL (de-duplicated by Entrez).
all_entrez <- unique(unlist(p2e, use.names = FALSE))
all_entrez <- all_entrez[!is.na(all_entrez) & nzchar(all_entrez)]
map <- suppressMessages(AnnotationDbi::select(
  db, keys = all_entrez, columns = "SYMBOL", keytype = "ENTREZID"))
map <- map[!is.na(map$SYMBOL), ]
map <- map[!duplicated(map$ENTREZID), ]
sym_of <- setNames(map$SYMBOL, map$ENTREZID)

# Emit long-form pathway_id / pathway_name / symbol rows.
con <- file(out_file, "w")
writeLines("pathway_id\tpathway_name\tsymbol", con)
n_pairs <- 0L
for (id in names(p2e)) {
  ent  <- unique(p2e[[id]])
  syms <- unique(sym_of[ent])
  syms <- syms[!is.na(syms)]
  if (length(syms) == 0) next
  nm <- name_of[[id]]
  for (s in syms) {
    writeLines(paste(id, nm, s, sep = "\t"), con)
    n_pairs <- n_pairs + 1L
  }
}
close(con)

manual_json <- function(x) {
  parts <- vapply(names(x), function(k)
    sprintf('"%s": "%s"', k, gsub('"', '\\"', as.character(x[[k]]))), character(1))
  paste0("{", paste(parts, collapse = ", "), "}")
}
prov <- list(
  reactome_db_version = as.character(packageVersion("reactome.db")),
  org_db_package      = pkg,
  org_db_version      = as.character(packageVersion(pkg)),
  r_version           = paste(R.version$major, R.version$minor, sep = "."),
  species_code        = species,
  n_pathways          = as.character(length(p2e)),
  n_pairs             = as.character(n_pairs)
)
writeLines(manual_json(prov), prov_file)
"""


@dataclass
class ReactomeGeneSets:
    """Materialised Reactome gene sets for one organism (single family).

    Attributes:
        taxid: Canonical NCBI taxid the gene sets were built for.
        species_code: The Reactome stable-id species code (e.g. ``HSA``).
        org_db_package: The ``org.*.eg.db`` used for the Entrez->SYMBOL join.
        term2gene: ``{pathway_id: set(symbol)}`` (the single REACTOME family).
        term2name: ``{pathway_id: pathway_name}``.
        provenance: Version pins (reactome.db, org.db package + version, R).
    """

    taxid: int
    species_code: str
    org_db_package: str
    term2gene: dict[str, set[str]]
    term2name: dict[str, str]
    provenance: dict[str, str] = field(default_factory=dict)

    def family(self, family: str = REACTOME_FAMILY) -> dict[str, set[str]]:
        """Return the ``pathway_id -> set(symbol)`` map (only REACTOME exists)."""
        if family != REACTOME_FAMILY:
            raise KeyError(f"no Reactome family {family!r}; only {REACTOME_FAMILY!r}")
        return self.term2gene

    def universe(self) -> set[str]:
        """Return the background = union of all symbols across all pathways."""
        genes: set[str] = set()
        for s in self.term2gene.values():
            genes |= s
        return genes

    def summary(self) -> dict[str, int]:
        """``{n_terms, n_genes, n_pairs}`` count summary (for manifests)."""
        genes: set[str] = set()
        pairs = 0
        for g in self.term2gene.values():
            genes |= g
            pairs += len(g)
        return {"n_terms": len(self.term2gene), "n_genes": len(genes), "n_pairs": pairs}


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


def species_code_for_taxid(taxid: int) -> str:
    """Return the NATIVE / curated Reactome stable-id species code for a taxid.

    Raises:
        ValueError: if Reactome ships no *native / curated* pathways for ``taxid``
            (e.g. macaque 9544) -- the caller must treat that organism as native
            Reactome = none. Inferred macaque Reactome is built separately by
            :mod:`mirpathdb_prod.pathways.inferred_reactome_genesets`.
    """
    try:
        return REACTOME_SPECIES_CODE[taxid]
    except KeyError as exc:
        raise ValueError(
            f"Reactome ships no native / curated pathways for taxid {taxid}; "
            f"available: {sorted(REACTOME_SPECIES_CODE)} "
            f"(no curated R-MML- -> macaque has no native Reactome; "
            f"inferred macaque Reactome is built by inferred_reactome_genesets)"
        ) from exc


def _parse_reactome_tsv(path: Path) -> tuple[dict[str, set[str]], dict[str, str]]:
    """Parse the R-written long-form Reactome TSV into (term2gene, term2name)."""
    term2gene: dict[str, set[str]] = {}
    term2name: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pid = (row.get("pathway_id") or "").strip()
            symbol = (row.get("symbol") or "").strip()
            if not pid or not symbol:
                continue
            term2gene.setdefault(pid, set()).add(symbol)
            name = (row.get("pathway_name") or "").strip()
            if name:
                term2name.setdefault(pid, name)
    return term2gene, term2name


def build_reactome_genesets(taxid: int, *, rscript: str | None = None) -> ReactomeGeneSets:
    """Materialise Reactome gene sets (SYMBOL space) for ``taxid`` from reactome.db.

    Shells out to R once, slices ``reactomePATHID2EXTID`` to the organism's
    ``R-<CODE>-`` pathways, joins Entrez members to SYMBOL via the organism's
    ``org.*.eg.db``, and returns the symbol-space single-family TERM2GENE plus a
    TERM2NAME with the ``"<Species>: "`` label stripped.

    Args:
        taxid: Canonical NCBI taxid (must be in :data:`REACTOME_SPECIES_CODE`).
        rscript: Override the Rscript binary (else $MPD_RSCRIPT / PATH).

    Returns:
        A :class:`ReactomeGeneSets` with the symbol-space TERM2GENE/TERM2NAME and
        the reactome.db/org.db/R version pins in ``provenance``.

    Raises:
        ValueError: if ``taxid`` has no native/curated Reactome coverage (e.g.
            macaque -- whose Reactome is built INFERRED by
            :mod:`mirpathdb_prod.pathways.inferred_reactome_genesets`).
        RuntimeError: if the R subprocess fails (stderr is surfaced).
    """
    species = species_code_for_taxid(taxid)
    pkg = org_db_for_taxid(taxid)
    rscript_bin = _resolve_rscript(rscript)

    with tempfile.TemporaryDirectory(prefix="mpd_reactome_t2g_") as tmp:
        tmp_path = Path(tmp)
        out_file = tmp_path / "reactome_term2gene.tsv"
        prov_file = tmp_path / "provenance.json"
        driver_file = tmp_path / "driver_reactome.R"
        driver_file.write_text(_R_REACTOME_DRIVER, encoding="utf-8")

        proc = subprocess.run(
            [rscript_bin, str(driver_file), species, pkg, str(out_file), str(prov_file)],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"Reactome TERM2GENE Rscript failed for {species}/{pkg} "
                f"(exit {proc.returncode}).\n--- stderr ---\n{proc.stderr}\n"
                f"--- stdout ---\n{proc.stdout}"
            )

        term2gene, term2name = _parse_reactome_tsv(out_file)
        provenance: dict[str, str] = {}
        if prov_file.exists():
            try:
                provenance = json.loads(prov_file.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                provenance = {}

    return ReactomeGeneSets(
        taxid=taxid,
        species_code=species,
        org_db_package=pkg,
        term2gene=term2gene,
        term2name=term2name,
        provenance=provenance,
    )


def write_reactome_genesets_tsv(gene_sets: ReactomeGeneSets, out_path: str | Path) -> Path:
    """Persist a materialised Reactome TERM2GENE as a long-form TSV (cache/audit).

    Columns: ``family<TAB>pathway_id<TAB>pathway_name<TAB>symbol`` (one row per
    pair), deterministic order. :func:`read_reactome_genesets_tsv` reloads it.
    """
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", quoting=csv.QUOTE_NONE, escapechar="\\")
        writer.writerow(["family", "pathway_id", "pathway_name", "symbol"])
        for pid in sorted(gene_sets.term2gene):
            name = gene_sets.term2name.get(pid, pid)
            for symbol in sorted(gene_sets.term2gene[pid]):
                writer.writerow([REACTOME_FAMILY, pid, name, symbol])
    return out


def read_reactome_genesets_tsv(
    path: str | Path, *, taxid: int, species_code: str, org_db_package: str
) -> ReactomeGeneSets:
    """Reload a cached Reactome TERM2GENE TSV (written by the writer above)."""
    term2gene: dict[str, set[str]] = {}
    term2name: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pid = (row.get("pathway_id") or "").strip()
            symbol = (row.get("symbol") or "").strip()
            if not pid or not symbol:
                continue
            term2gene.setdefault(pid, set()).add(symbol)
            name = (row.get("pathway_name") or "").strip()
            if name:
                term2name.setdefault(pid, name)
    return ReactomeGeneSets(
        taxid=taxid,
        species_code=species_code,
        org_db_package=org_db_package,
        term2gene=term2gene,
        term2name=term2name,
    )


__all__ = [
    "ORG_DB_PACKAGES",
    "REACTOME_FAMILY",
    "REACTOME_SPECIES_CODE",
    "ReactomeGeneSets",
    "build_reactome_genesets",
    "read_reactome_genesets_tsv",
    "species_code_for_taxid",
    "write_reactome_genesets_tsv",
]
