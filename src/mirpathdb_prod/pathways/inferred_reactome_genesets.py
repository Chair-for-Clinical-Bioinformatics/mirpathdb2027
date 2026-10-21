"""INFERRED macaque Reactome gene-sets (human pathways projected via orthology).

Reactome ships **no curated / authored Macaca mulatta pathways** -- the
Bioconductor ``reactome.db`` v1.95.0 and the flat ``Ensembl2Reactome_All_Levels``
file both carry **0** ``R-MML-`` rows -- which is why
:mod:`mirpathdb_prod.pathways.reactome_genesets` deliberately excludes macaque
(9544) from its native species map. Reactome itself, however, *does* list
``Macaca mulatta`` (``R-MML-`` stable ids) as an **ORTHOLOGY-INFERRED** species:
its macaque pathways are not human-curated, they are *projected* from the human
(``R-HSA-``) pathways through human<->macaque gene orthology.

This module reproduces **Reactome's own inference method** for ``R-MML-`` so the
production build can include macaque Reactome as an **INFERRED, orthology-
projected, LIMITED** family (distinct from the native, curated Reactome families
hsa/mmu/rno get). The recipe::

    human Reactome pathway gene set (Entrez, from reactome.db, R-HSA-)
      -> hsa Entrez  --[NCBI gene_orthologs, 9606<->9544, strict 1:1]-->  mml Entrez
      -> mml SYMBOL  (org.Mmu.eg.db)
      = the INFERRED macaque Reactome gene set for that pathway

Every emitted pathway/output is **marked INFERRED** (orthology-projected from
human, limited): the inferred pathway id is ``R-MML-<n>`` (the macaque analogue of
the human ``R-HSA-<n>``; ``REACTOME_INFERRED_PREFIX``), the gene sets are in the
macaque ``org.Mmu.eg.db`` SYMBOL space (so they live in the SAME symbol space as
the macaque GO + TargetScan target sets), and the provenance records
``inferred = True`` plus the pinned NCBI ``gene_orthologs`` source (size + sha256
+ Last-Modified) alongside the reactome.db / org.db / R version pins.

Reproducibility: the ortholog source is the **cached** NCBI
``gene_orthologs.gz`` (no live API call); its version is pinned by file size,
sha256, and the recorded NCBI ``Last-Modified`` date. The pathway-membership
source is the versioned Bioconductor ``reactome.db`` (R-HSA- pathways), and both
the human and macaque Entrez->SYMBOL joins go through the SAME ``org.*.eg.db``
release the GO sets came from -- so the inferred macaque Reactome cannot drift
from the rest of the macaque symbol space.

This module shells out to R (like the native Reactome + GO sources); it does not
import ``rpy2``. The R binary is ``$MPD_RSCRIPT`` if set, else the first
``Rscript`` on ``PATH``.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import warnings
from dataclasses import dataclass, field
from pathlib import Path

from mirpathdb_prod.pathways.go_genesets import org_db_for_taxid

# The single pathway family this source emits (same token as native Reactome, so
# the downstream ORA / matrix naming is identical -- the INFERRED marking lives in
# the provenance + manifest, NOT the family token). Re-exported from the native
# Reactome module so the token cannot drift between native and inferred Reactome.
from mirpathdb_prod.pathways.reactome_genesets import REACTOME_FAMILY

# Canonical taxids for the projection (human source -> macaque target).
HUMAN_TAXID = 9606
MACAQUE_TAXID = 9544

# Human Reactome stable-id prefix we project FROM, and the inferred macaque
# stable-id prefix we relabel TO. ``R-MML-`` is Reactome's own code for Macaca
# mulatta (the orthology-inferred species), so an inferred pathway keeps the
# human pathway number under the macaque species code: R-HSA-<n> -> R-MML-<n>.
HUMAN_PREFIX = "R-HSA-"
REACTOME_INFERRED_PREFIX = "R-MML-"

# NCBI taxids used to slice gene_orthologs (human <-> macaque).
HUMAN_NCBI_TAXID = "9606"
MACAQUE_NCBI_TAXID = "9544"

# Default cached NCBI ortholog source (no live API call -- see module docstring).
DEFAULT_ORTHOLOG_PATH = "/tmp/mpd-data/gene_orthologs.gz"

RSCRIPT_ENV = "MPD_RSCRIPT"

# The R driver. Given (human org.db, macaque org.db, gene_orthologs path), it:
#   1. pulls every R-HSA- Reactome pathway -> Entrez members from reactome.db;
#   2. reads gene_orthologs, slices the strict 9606<->9544 ortholog pairs into an
#      hsa-Entrez -> mml-Entrez map;
#   3. joins mml Entrez -> mml SYMBOL via org.Mmu.eg.db;
#   4. for each human pathway, projects its hsa Entrez members through the
#      ortholog map to mml Entrez, then to mml SYMBOL = the INFERRED macaque set;
#   5. relabels R-HSA-<n> -> R-MML-<n> and writes a long-form TSV
#      (pathway_id<TAB>pathway_name<TAB>symbol) + a provenance JSON.
# The pathway NAME comes from the human reactomePATHID2NAME with the leading
# "Homo sapiens: " label stripped (the macaque inferred pathway carries the human
# pathway's name -- it IS the human pathway, orthology-projected).
_R_INFERRED_DRIVER = r"""
suppressMessages({
  library(AnnotationDbi)
  library(reactome.db)
})
args        <- commandArgs(trailingOnly = TRUE)
hsa_pkg     <- args[[1]]              # "org.Hs.eg.db"
mml_pkg     <- args[[2]]              # "org.Mmu.eg.db"
orth_file   <- args[[3]]              # gene_orthologs.gz
human_prefix<- args[[4]]              # "R-HSA-"
infer_prefix<- args[[5]]              # "R-MML-"
out_file    <- args[[6]]
prov_file   <- args[[7]]

suppressMessages(library(hsa_pkg, character.only = TRUE))
suppressMessages(library(mml_pkg, character.only = TRUE))
mml_db <- get(mml_pkg)

# 1. human R-HSA- Reactome pathway -> Entrez members (restrict to mappedkeys so
#    the subset cannot error on container/event nodes with no EXTID).
mapped <- mappedkeys(reactomePATHID2EXTID)
ids    <- mapped[startsWith(mapped, human_prefix)]
p2e    <- as.list(reactomePATHID2EXTID[ids])
p2e    <- p2e[!vapply(p2e, is.null, logical(1))]

# human pathway -> name; strip the leading "<Species>: " label.
named_ids <- intersect(names(p2e), mappedkeys(reactomePATHID2NAME))
nm_list   <- if (length(named_ids)) as.list(reactomePATHID2NAME[named_ids]) else list()
name_of <- vapply(names(p2e), function(id) {
  v <- nm_list[[id]]
  if (is.null(v) || length(v) == 0 || is.na(v[[1]])) return(id)
  sub("^[^:]+:\\s*", "", v[[1]])
}, character(1))

# 2. gene_orthologs -> strict 9606<->9544 hsa-Entrez -> mml-Entrez map. The NCBI
#    gene_orthologs file lists ONLY one-to-one orthologs; we still de-dup defensively.
orth <- read.delim(gzfile(orth_file), header = TRUE, sep = "\t",
                   colClasses = "character", quote = "", comment.char = "")
colnames(orth)[1] <- "tax_id"
sel <- orth$tax_id == "9606" & orth$Other_tax_id == "9544"
hsa_e <- orth$GeneID[sel]
mml_e <- orth$Other_GeneID[sel]
keep  <- !duplicated(hsa_e) & !duplicated(mml_e)   # strict 1:1
hsa2mml <- setNames(mml_e[keep], hsa_e[keep])
n_orth  <- length(hsa2mml)

# 3. mml Entrez -> mml SYMBOL (org.Mmu.eg.db), de-duplicated by Entrez.
mml_entrez <- unique(unname(hsa2mml))
mml_map <- suppressMessages(AnnotationDbi::select(
  mml_db, keys = mml_entrez, columns = "SYMBOL", keytype = "ENTREZID"))
mml_map <- mml_map[!is.na(mml_map$SYMBOL), ]
mml_map <- mml_map[!duplicated(mml_map$ENTREZID), ]
mml_sym <- setNames(mml_map$SYMBOL, mml_map$ENTREZID)

# 4 + 5. project each human pathway and relabel R-HSA-<n> -> R-MML-<n>.
con <- file(out_file, "w")
writeLines("pathway_id\tpathway_name\tsymbol", con)
n_terms <- 0L; n_pairs <- 0L
for (id in names(p2e)) {
  hsa_ent <- unique(p2e[[id]])
  m_ent   <- unique(hsa2mml[hsa_ent]); m_ent <- m_ent[!is.na(m_ent)]
  syms    <- unique(mml_sym[m_ent]);   syms  <- syms[!is.na(syms)]
  if (length(syms) == 0) next
  infer_id <- sub(paste0("^", human_prefix), infer_prefix, id)
  nm <- name_of[[id]]
  for (s in syms) {
    writeLines(paste(infer_id, nm, s, sep = "\t"), con)
    n_pairs <- n_pairs + 1L
  }
  n_terms <- n_terms + 1L
}
close(con)

manual_json <- function(x) {
  parts <- vapply(names(x), function(k)
    sprintf('"%s": "%s"', k, gsub('"', '\\"', as.character(x[[k]]))), character(1))
  paste0("{", paste(parts, collapse = ", "), "}")
}
prov <- list(
  inferred              = "true",
  inference_method      = paste0(
    "human R-HSA- Reactome projected to macaque via NCBI 1:1 orthology ",
    "(replicates Reactome's R-MML- inference)"),
  source_prefix         = human_prefix,
  inferred_prefix       = infer_prefix,
  reactome_db_version   = as.character(packageVersion("reactome.db")),
  human_org_db_package  = hsa_pkg,
  human_org_db_version  = as.character(packageVersion(hsa_pkg)),
  macaque_org_db_package= mml_pkg,
  macaque_org_db_version= as.character(packageVersion(mml_pkg)),
  r_version             = paste(R.version$major, R.version$minor, sep = "."),
  n_orthologs           = as.character(n_orth),
  n_human_pathways      = as.character(length(p2e)),
  n_inferred_pathways   = as.character(n_terms),
  n_pairs               = as.character(n_pairs)
)
writeLines(manual_json(prov), prov_file)
"""


@dataclass
class InferredReactomeGeneSets:
    """INFERRED macaque Reactome gene sets (orthology-projected from human).

    Mirrors :class:`mirpathdb_prod.pathways.reactome_genesets.ReactomeGeneSets`
    but is explicitly **inferred**: ``inferred`` is always ``True`` and the
    provenance carries the ortholog-source pin + the inference method.

    Attributes:
        taxid: The macaque taxid (always 9544 -- the projection target).
        species_code: The inferred Reactome species code (``MML``).
        org_db_package: The macaque ``org.Mmu.eg.db`` used for Entrez->SYMBOL.
        term2gene: ``{R-MML-<n>: set(macaque SYMBOL)}`` (the single REACTOME family).
        term2name: ``{R-MML-<n>: human pathway name}`` (carried from the human pathway).
        provenance: Version pins (reactome.db, both org.db pkgs + versions, R) plus
            the inference metadata and the pinned NCBI ortholog source.
        inferred: Always ``True`` -- this family is orthology-projected, not curated.
    """

    taxid: int
    species_code: str
    org_db_package: str
    term2gene: dict[str, set[str]]
    term2name: dict[str, str]
    provenance: dict[str, str] = field(default_factory=dict)
    inferred: bool = True

    def family(self, family: str = REACTOME_FAMILY) -> dict[str, set[str]]:
        """Return the ``R-MML-<n> -> set(symbol)`` map (only REACTOME exists)."""
        if family != REACTOME_FAMILY:
            raise KeyError(f"no Reactome family {family!r}; only {REACTOME_FAMILY!r}")
        return self.term2gene

    def universe(self) -> set[str]:
        """Return the background = union of all inferred macaque symbols."""
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


def _file_sha256(path: Path, *, chunk: int = 1 << 20) -> str:
    """Stream a sha256 of a (possibly large) file for ortholog-source pinning."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def ortholog_source_provenance(
    ortholog_path: str | Path, *, last_modified: str | None = None
) -> dict[str, str]:
    """Pin the NCBI ``gene_orthologs`` source (path, size, sha256, Last-Modified).

    The macaque Reactome is INFERRED, so the ortholog source is the load-bearing
    input that must be reproducible; this records its identity so a rebuild from
    the same file is verifiable.

    Args:
        ortholog_path: Path to the cached ``gene_orthologs.gz``.
        last_modified: The NCBI ``Last-Modified`` date (HTTP header) if known;
            recorded verbatim so the upstream snapshot is pinned.

    Returns:
        ``{source, path, size_bytes, sha256, ncbi_last_modified, taxids}``.
    """
    p = Path(ortholog_path)
    if not p.exists():
        raise FileNotFoundError(f"NCBI gene_orthologs not found: {p}")
    prov: dict[str, str] = {
        "source": "NCBI Gene gene_orthologs.gz (ftp.ncbi.nlm.nih.gov/gene/DATA/)",
        "path": str(p),
        "size_bytes": str(p.stat().st_size),
        "sha256": _file_sha256(p),
        "taxids": f"{HUMAN_NCBI_TAXID}<->{MACAQUE_NCBI_TAXID} (one-to-one)",
    }
    if last_modified:
        prov["ncbi_last_modified"] = last_modified
    return prov


def _parse_inferred_tsv(path: Path) -> tuple[dict[str, set[str]], dict[str, str]]:
    """Parse the R-written long-form inferred-Reactome TSV -> (term2gene, term2name)."""
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


def build_inferred_macaque_reactome(
    *,
    ortholog_path: str | Path = DEFAULT_ORTHOLOG_PATH,
    ortholog_last_modified: str | None = None,
    rscript: str | None = None,
) -> InferredReactomeGeneSets:
    """Materialise the INFERRED macaque Reactome gene sets (orthology projection).

    Projects every human ``R-HSA-`` Reactome pathway onto macaque through the NCBI
    ``gene_orthologs`` 9606<->9544 strict 1:1 ortholog map (mapped to macaque
    SYMBOL via ``org.Mmu.eg.db``), relabelling each pathway ``R-HSA-<n>`` ->
    ``R-MML-<n>``. This replicates Reactome's OWN inference method for ``R-MML-``.

    Args:
        ortholog_path: Cached NCBI ``gene_orthologs.gz`` (default the build cache).
        ortholog_last_modified: NCBI ``Last-Modified`` date for the file (pinned).
        rscript: Override the Rscript binary (else $MPD_RSCRIPT / PATH).

    Returns:
        An :class:`InferredReactomeGeneSets` (``inferred=True``) with the
        macaque-SYMBOL-space TERM2GENE/TERM2NAME and a provenance dict carrying
        the inference metadata + the pinned reactome.db / org.db / R versions +
        the pinned NCBI ortholog source.

    Raises:
        FileNotFoundError: if the ortholog file is missing.
        RuntimeError: if the R subprocess fails (stderr is surfaced).
    """
    orth = Path(ortholog_path)
    if not orth.exists():
        raise FileNotFoundError(f"NCBI gene_orthologs not found: {orth}")
    hsa_pkg = org_db_for_taxid(HUMAN_TAXID)
    mml_pkg = org_db_for_taxid(MACAQUE_TAXID)
    rscript_bin = _resolve_rscript(rscript)

    with tempfile.TemporaryDirectory(prefix="mpd_inferred_reactome_") as tmp:
        tmp_path = Path(tmp)
        out_file = tmp_path / "inferred_reactome_term2gene.tsv"
        prov_file = tmp_path / "provenance.json"
        driver_file = tmp_path / "driver_inferred_reactome.R"
        driver_file.write_text(_R_INFERRED_DRIVER, encoding="utf-8")

        proc = subprocess.run(
            [
                rscript_bin,
                str(driver_file),
                hsa_pkg,
                mml_pkg,
                str(orth),
                HUMAN_PREFIX,
                REACTOME_INFERRED_PREFIX,
                str(out_file),
                str(prov_file),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"inferred macaque Reactome Rscript failed (exit {proc.returncode}).\n"
                f"--- stderr ---\n{proc.stderr}\n--- stdout ---\n{proc.stdout}"
            )

        term2gene, term2name = _parse_inferred_tsv(out_file)
        provenance: dict[str, str] = {}
        if prov_file.exists():
            try:
                provenance = json.loads(prov_file.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                # L5: surface a lost reactome.db/org.db/R version pin, don't drop {}.
                warnings.warn(
                    f"inferred-Reactome provenance JSON unreadable ({prov_file}): "
                    f"{exc} -- recording empty provenance (version pins lost).",
                    stacklevel=2,
                )
                provenance = {}
        else:
            warnings.warn(
                f"inferred-Reactome provenance JSON not written ({prov_file}); "
                f"reactome.db/org.db/R version pins unavailable.",
                stacklevel=2,
            )

    # Always stamp the ortholog-source pin + the inferred flag (independent of R).
    provenance["inferred"] = "true"
    provenance["ortholog_source"] = json.dumps(
        ortholog_source_provenance(orth, last_modified=ortholog_last_modified)
    )

    return InferredReactomeGeneSets(
        taxid=MACAQUE_TAXID,
        species_code="MML",
        org_db_package=mml_pkg,
        term2gene=term2gene,
        term2name=term2name,
        provenance=provenance,
        inferred=True,
    )


def write_inferred_reactome_genesets_tsv(
    gene_sets: InferredReactomeGeneSets, out_path: str | Path
) -> Path:
    """Persist an inferred macaque Reactome TERM2GENE as a long-form TSV (cache/audit).

    Columns: ``family<TAB>pathway_id<TAB>pathway_name<TAB>symbol`` (one row per
    pair), deterministic order. :func:`read_inferred_reactome_genesets_tsv` reloads it.
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


def read_inferred_reactome_genesets_tsv(
    path: str | Path, *, org_db_package: str | None = None
) -> InferredReactomeGeneSets:
    """Reload a cached inferred macaque Reactome TERM2GENE TSV (written by the writer)."""
    term2gene, term2name = _parse_inferred_tsv(Path(path))
    return InferredReactomeGeneSets(
        taxid=MACAQUE_TAXID,
        species_code="MML",
        org_db_package=org_db_package or org_db_for_taxid(MACAQUE_TAXID),
        term2gene=term2gene,
        term2name=term2name,
        inferred=True,
    )


__all__ = [
    "DEFAULT_ORTHOLOG_PATH",
    "HUMAN_PREFIX",
    "MACAQUE_TAXID",
    "REACTOME_FAMILY",
    "REACTOME_INFERRED_PREFIX",
    "InferredReactomeGeneSets",
    "build_inferred_macaque_reactome",
    "ortholog_source_provenance",
    "read_inferred_reactome_genesets_tsv",
    "write_inferred_reactome_genesets_tsv",
]
