"""clusterProfiler over-representation engine (Bioconductor / R bridge).

clusterProfiler is the production GO ORA engine, replacing the retired custom
ORA/BH kernel. It runs in R via a small ``Rscript`` subprocess that calls
``clusterProfiler::enricher`` with a caller-supplied ``TERM2GENE`` — so the same
bridge serves both the **controlled v2 reproduction** (TERM2GENE = reconstructed
v2 gene sets, no annotation DB needed; the space is whatever the TERM2GENE uses,
typically gene symbols) and the **production GO build** (TERM2GENE materialised
from ``org.*.eg.db``).

Numerical contract (matches v2's GeneTrail2 ORA and the recon spec)::

    enricher(gene, TERM2GENE, TERM2NAME, universe,
             pvalueCutoff = 1, qvalueCutoff = 1,        # return ALL terms' p
             pAdjustMethod = "BH",
             minGSSize = 1, maxGSSize = <very large>)

``pvalueCutoff = qvalueCutoff = 1`` is essential: it makes every testable term
come back with a p-value (not just the significant ones), which is what a
cell-for-cell comparison against a dense v2 matrix needs.

Performance — BATCHED + CORE-PARALLEL (the workhorse path)
----------------------------------------------------------
A v2 reproduction cell tests hundreds of miRNAs against ONE TERM2GENE. Shelling
out one ``Rscript`` per miRNA (the :func:`run_ora` single-set path) pays the R
startup + ``clusterProfiler`` attach + TERM2GENE-parse cost once *per miRNA* —
hours of serial subprocess churn. :func:`run_ora_batch` is the production path:
it writes every gene set + the shared TERM2GENE + universe ONCE, then a SINGLE
R session loads ``clusterProfiler`` once, builds the TERM2GENE/TERM2NAME frames
once, and runs ``parallel::mclapply`` over all gene sets across
``min(detectCores() - 4, $MPD_MAX_CORES)`` cores (default cap 96). The full
hsa+mmu controlled reproduction then completes in minutes, not hours.

:func:`run_ora` is a thin one-entry wrapper over :func:`run_ora_batch` so the
recon single-set interface and the batch path share one numerical driver.

Reproducibility (REQUIRED, see ``docs/enrichment.md``): provenance pins the
clusterProfiler + R versions on every run, and for the production GO build the
``org.*.eg.db`` annotation-package version (GO annotations drift between
Bioconductor releases). The controlled reproduction does not consult an
annotation DB, so its provenance records ``org_db_package = None``.

This module shells out to R; it does not import ``rpy2``. The R binary is taken
from ``$MPD_RSCRIPT`` if set, else the first ``Rscript`` on ``PATH``.
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import tempfile
import warnings
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

from mirpathdb_prod.enrichment.base import (
    ENGINE_CLUSTERPROFILER,
    EnrichmentRequest,
    EnrichmentResult,
)

# Canonical taxid -> Bioconductor GO-annotation package (production GO build).
# Its version MUST be pinned in provenance for a reproducible production run.
ORG_DB_PACKAGES: dict[int, str] = {
    9606: "org.Hs.eg.db",  # Homo sapiens
    10090: "org.Mm.eg.db",  # Mus musculus
    10116: "org.Rn.eg.db",  # Rattus norvegicus
    9544: "org.Mmu.eg.db",  # Macaca mulatta (macaque — predicted-dominant; lone validated MTI)
}

# Provenance pin slots a real run must populate.
PROVENANCE_PINS: tuple[str, ...] = (
    "clusterprofiler_version",
    "r_version",
    "org_db_package",
    "org_db_version",
    "go_release",
    "bioconductor_release",
)

# Environment override for the R interpreter; falls back to PATH lookup.
RSCRIPT_ENV = "MPD_RSCRIPT"

# Environment cap on the mclapply core count for the batched driver (the host
# may have far more cores than we want one cell to grab). Default 96 per the
# performance mandate; the driver further caps at detectCores() - 4.
MAX_CORES_ENV = "MPD_MAX_CORES"
_DEFAULT_MAX_CORES = 96
_CORE_HEADROOM = 4  # leave a few cores for the parent / OS under mclapply forks

# A "no upper bound" gene-set size for enricher (R has no clean +Inf for an
# integer arg; a value far larger than any real term suffices and is explicit).
_MAX_GS_SIZE = 100_000_000

# The R driver. Reads three TSVs (gene set, TERM2GENE, universe) + an optional
# TERM2NAME TSV, runs enricher with the locked parameters, and writes the result
# table + a small JSON provenance blob. Kept inline so there is no extra file to
# ship; written to a temp dir per call.
_R_DRIVER = r"""
suppressMessages({
  library(clusterProfiler)
})

args <- commandArgs(trailingOnly = TRUE)
gene_file      <- args[[1]]
t2g_file       <- args[[2]]
universe_file  <- args[[3]]
t2n_file       <- args[[4]]   # may be "NA" -> no TERM2NAME
out_file       <- args[[5]]
prov_file      <- args[[6]]

read_one_col <- function(path) {
  v <- scan(path, what = character(), sep = "\n", quiet = TRUE)
  v <- trimws(v)
  v[nzchar(v)]
}

gene     <- read_one_col(gene_file)
universe <- read_one_col(universe_file)

t2g <- read.delim(t2g_file, header = FALSE, sep = "\t",
                  colClasses = "character", quote = "",
                  stringsAsFactors = FALSE)
colnames(t2g) <- c("term", "gene")

t2n <- NULL
if (t2n_file != "NA") {
  t2n <- read.delim(t2n_file, header = FALSE, sep = "\t",
                    colClasses = "character", quote = "",
                    stringsAsFactors = FALSE)
  colnames(t2n) <- c("term", "name")
}

res <- enricher(
  gene          = gene,
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  pAdjustMethod = "BH",
  universe      = universe,
  minGSSize     = 1,
  maxGSSize     = %MAX_GS_SIZE%,
  TERM2GENE     = t2g,
  TERM2NAME     = t2n
)

df <- as.data.frame(res)
if (nrow(df) == 0L) {
  out <- data.frame(ID = character(0), Description = character(0),
                    pvalue = numeric(0), p.adjust = numeric(0),
                    Count = integer(0))
} else {
  # Description == ID when no TERM2NAME was supplied (enricher's default).
  out <- df[, c("ID", "Description", "pvalue", "p.adjust", "Count")]
}
write.table(out, file = out_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

prov <- list(
  clusterprofiler_version = as.character(packageVersion("clusterProfiler")),
  r_version = paste(R.version$major, R.version$minor, sep = ".")
)
writeLines(jsonlite_or_manual(prov), prov_file)
"""

# clusterProfiler depends on jsonlite transitively, but to avoid assuming it is
# attached we emit the tiny provenance JSON by hand in R.
_R_JSON_HELPER = r"""
jsonlite_or_manual <- function(x) {
  parts <- vapply(names(x), function(k) {
    sprintf('"%s": "%s"', k, gsub('"', '\\"', as.character(x[[k]])))
  }, character(1))
  paste0("{", paste(parts, collapse = ", "), "}")
}
"""

# The BATCHED + PARALLEL R driver (the workhorse). It loads clusterProfiler ONCE,
# parses TERM2GENE / TERM2NAME / universe ONCE, splits the long-form gene-set
# table into a per-name list, then mclapply()s enricher over every gene set with
# the locked parameters, and writes one combined result table (name, ID,
# Description, pvalue, p.adjust, Count). mc.cores comes from the parent (already
# capped at min($MPD_MAX_CORES, detectCores() - 4)).
_R_BATCH_DRIVER = r"""
suppressMessages({
  library(clusterProfiler)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
genesets_file <- args[[1]]   # 2-col TSV: name<TAB>gene  (long form, all sets)
t2g_file      <- args[[2]]
universe_file <- args[[3]]
t2n_file      <- args[[4]]   # may be "NA" -> no TERM2NAME
out_file      <- args[[5]]
prov_file     <- args[[6]]
mc_cores      <- as.integer(args[[7]])
core_headroom <- as.integer(args[[8]])

# Enforce the documented cap min(requested, detectCores() - headroom), never < 1.
# The clamp runs HERE (in R): detectCores() is the enrichment host's true core
# count, which the Python side cannot know. This leaves `core_headroom` cores for
# the parent/OS under mclapply forks so a host with fewer cores than the default
# request (96) is not oversubscribed.
.avail <- parallel::detectCores()
if (!is.na(.avail) && .avail > 0L) {
  mc_cores <- max(1L, min(mc_cores, .avail - core_headroom))
} else {
  mc_cores <- max(1L, mc_cores)
}

read_one_col <- function(path) {
  if (file.info(path)$size == 0) return(character(0))
  v <- scan(path, what = character(), sep = "\n", quiet = TRUE)
  v <- trimws(v)
  v[nzchar(v)]
}

universe <- read_one_col(universe_file)

t2g <- read.delim(t2g_file, header = FALSE, sep = "\t",
                  colClasses = "character", quote = "",
                  stringsAsFactors = FALSE)
colnames(t2g) <- c("term", "gene")

t2n <- NULL
if (t2n_file != "NA") {
  t2n <- read.delim(t2n_file, header = FALSE, sep = "\t",
                    colClasses = "character", quote = "",
                    stringsAsFactors = FALSE)
  colnames(t2n) <- c("term", "name")
}

# Long-form (name, gene) -> named list of gene-symbol character vectors. Done
# ONCE; reused across every enricher call.
gs <- read.delim(genesets_file, header = FALSE, sep = "\t",
                 colClasses = "character", quote = "",
                 stringsAsFactors = FALSE)
colnames(gs) <- c("name", "gene")
gene_sets <- split(gs$gene, gs$name)
set_names <- names(gene_sets)

# A worker returns one of three things, which we MUST keep distinct (B2):
#   * a data.frame of hits           -> a scored set;
#   * the sentinel "__NO_OVERLAP__"  -> a LEGITIMATE no-overlap set (enricher ran,
#                                       found nothing testable / nrow(df)==0);
#   * the sentinel "__ERROR__"       -> enricher actually THREW (a real failure).
# Distinguishing the last two is the whole point: the old code mapped BOTH to NULL
# and could exit 0 with an empty table even when every set errored.
NO_OVERLAP <- "__NO_OVERLAP__"
run_one <- function(nm) {
  res <- tryCatch(
    enricher(
      gene          = gene_sets[[nm]],
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      pAdjustMethod = "BH",
      universe      = universe,
      minGSSize     = 1,
      maxGSSize     = %MAX_GS_SIZE%,
      TERM2GENE     = t2g,
      TERM2NAME     = t2n
    ),
    error = function(e) structure(conditionMessage(e), class = "mpd_worker_error")
  )
  if (inherits(res, "mpd_worker_error")) return(res)        # real enricher failure
  if (is.null(res)) return(NO_OVERLAP)                       # enricher -> NULL: no overlap
  df <- as.data.frame(res)
  if (nrow(df) == 0L) return(NO_OVERLAP)                     # legitimate no-overlap
  data.frame(
    name        = nm,
    ID          = df$ID,
    Description = df$Description,
    pvalue      = df$pvalue,
    p.adjust    = df$p.adjust,
    Count       = df$Count,
    stringsAsFactors = FALSE
  )
}

parts <- mclapply(set_names, run_one, mc.cores = mc_cores)

# A try-error slot means mclapply's fork itself crashed (e.g. the worker was
# killed) -- surface it immediately, as before.
crashed <- vapply(parts, function(p) inherits(p, "try-error"), logical(1))
# An mpd_worker_error means enricher() THREW for that set (now captured, not
# swallowed to NULL). Count these so a cell where EVERY set failed cannot exit 0
# with an empty table (the degenerate all-p=1.0 trap).
worker_err <- vapply(parts, function(p) inherits(p, "mpd_worker_error"), logical(1))
n_failed <- sum(crashed) + sum(worker_err)
n_total  <- length(parts)
if (n_total > 0L && n_failed == n_total) {
  msgs <- unique(as.character(unlist(parts[worker_err])))
  stop(sprintf(
    paste0("clusterProfiler batch: ALL %d gene set(s) failed in workers (none ",
           "scored) -- refusing an empty/degenerate result. First error(s): %s"),
    n_total, paste(head(msgs, 3), collapse = " | ")
  ))
}
if (any(crashed)) {
  stop(sprintf("mclapply: %d gene set(s) crashed in workers", sum(crashed)))
}

# Partial worker errors (some sets threw but NOT all) were previously dropped
# silently by the is_df filter below, so the affected miRNAs vanished from the
# result with no trace. Count + example them (recorded in provenance) and emit an
# R warning so a partial enricher failure is visible, not a silent hole. Computed
# BEFORE the is_df filter reassigns `parts`, while the error entries still exist.
n_worker_err <- sum(worker_err)
worker_err_examples <- ""
if (n_worker_err > 0L) {
  worker_err_examples <- paste(
    head(unique(as.character(unlist(parts[worker_err]))), 3), collapse = " | ")
  warning(sprintf(
    paste0("clusterProfiler batch: %d of %d gene set(s) threw in workers and ",
           "were dropped (kept the rest). First: %s"),
    n_worker_err, n_total, worker_err_examples))
}

# Keep only the real data.frame parts; both sentinels (no-overlap + the captured
# worker errors, which were NOT all-failures here) are dropped from the table.
is_df <- vapply(parts, is.data.frame, logical(1))
parts <- parts[is_df]

if (length(parts) == 0L) {
  out <- data.frame(name = character(0), ID = character(0),
                    Description = character(0), pvalue = numeric(0),
                    p.adjust = numeric(0), Count = integer(0))
} else {
  out <- do.call(rbind, parts)
}
write.table(out, file = out_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

prov <- list(
  clusterprofiler_version = as.character(packageVersion("clusterProfiler")),
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  mc_cores = as.character(mc_cores),
  detected_cores = as.character(.avail),
  n_gene_sets = as.character(length(set_names)),
  n_worker_errors = as.character(n_worker_err),
  worker_error_examples = worker_err_examples
)
writeLines(jsonlite_or_manual(prov), prov_file)
"""


def _resolve_rscript() -> str:
    """Return the Rscript binary path ($MPD_RSCRIPT or PATH lookup)."""
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


def _resolve_max_cores() -> int:
    """Return the mclapply core cap from ``$MPD_MAX_CORES`` (default 96).

    The R driver further clamps this to ``detectCores() - _CORE_HEADROOM`` at
    runtime, so this is only the Python-side ceiling; a malformed/zero env value
    falls back to the default.
    """
    raw = os.environ.get(MAX_CORES_ENV)
    if raw:
        try:
            val = int(raw)
            if val >= 1:
                return val
        except ValueError:
            pass
    return _DEFAULT_MAX_CORES


def _write_lines(path: Path, items: Iterable[str]) -> None:
    """Write one item per line (UTF-8), skipping blanks."""
    with path.open("w", encoding="utf-8") as fh:
        for item in items:
            s = str(item).strip()
            if s:
                fh.write(s + "\n")


def _write_pairs(path: Path, pairs: Iterable[tuple[str, str]]) -> None:
    """Write tab-separated (a, b) pairs (no header, csv-quoting disabled)."""
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", quoting=csv.QUOTE_NONE, escapechar="\\")
        for a, b in pairs:
            writer.writerow([a, b])


def _term2gene_pairs(term2gene: Mapping[str, Iterable[str]]) -> list[tuple[str, str]]:
    """Flatten ``term -> genes`` into deterministic ``(term, gene)`` pairs."""
    pairs: list[tuple[str, str]] = []
    for term in sorted(term2gene):
        for gene in sorted(set(term2gene[term])):
            pairs.append((term, gene))
    return pairs


def run_ora(
    gene_set: Iterable[str],
    term2gene: Mapping[str, Iterable[str]],
    universe: Iterable[str],
    *,
    term2name: Mapping[str, str] | None = None,
    rscript: str | None = None,
) -> dict[str, tuple[float, float, int]]:
    """Run clusterProfiler over-representation for ONE gene set.

    Honors the recon kernel interface. A thin one-entry wrapper over
    :func:`run_ora_batch` (single shared numerical driver). Runs ``enricher``
    with the locked parameters (``pvalueCutoff = qvalueCutoff = 1``,
    ``minGSSize = 1``, ``maxGSSize`` effectively unbounded, ``pAdjustMethod =
    "BH"``) and reads back ``ID / pvalue / p.adjust / Count``.

    Args:
        gene_set: The query gene identifiers (same ID space as ``term2gene``).
        term2gene: ``term -> iterable(gene)`` membership (drives the ID space;
            clusterProfiler's effective universe is the genes here ∩ ``universe``).
        universe: Background gene identifiers.
        term2name: Optional ``term -> display name`` (passed as TERM2NAME so the
            returned ``ID`` carries the name; defaults to identity term names).
        rscript: Override the Rscript binary (else $MPD_RSCRIPT / PATH).

    Returns:
        ``{term: (pvalue, p_adjust, overlap_count)}`` for every term enricher
        reports. Terms with no testable overlap are absent (consistent with
        clusterProfiler returning only terms it scored).

    Raises:
        RuntimeError: if the R subprocess fails (stderr is surfaced).
    """
    result, _provenance = _run_ora_full(
        gene_set, term2gene, universe, term2name=term2name, rscript=rscript
    )
    return result


# Sentinel gene-set name used by the single-set wrapper through the batch path.
_SINGLE_SET_NAME = "__single__"


def _run_ora_full(
    gene_set: Iterable[str],
    term2gene: Mapping[str, Iterable[str]],
    universe: Iterable[str],
    *,
    term2name: Mapping[str, str] | None = None,
    rscript: str | None = None,
) -> tuple[dict[str, tuple[float, float, int]], dict[str, str]]:
    """Run the enricher for one gene set, returning ``(result, r_provenance)``.

    Internal: delegates to :func:`run_ora_batch` (one entry) so the single-set
    and batch paths share one numerical driver. The batch driver writes a
    provenance JSON (clusterProfiler + R versions probed inside the *same*
    interpreter that ran the analysis — the authoritative pin). The public
    :func:`run_ora` discards the provenance to keep the recon return contract;
    :class:`ClusterProfilerEngine` keeps it.
    """
    batch, provenance = _run_ora_batch_full(
        {_SINGLE_SET_NAME: gene_set},
        term2gene,
        universe,
        term2name=term2name,
        rscript=rscript,
    )
    return batch.get(_SINGLE_SET_NAME, {}), provenance


def run_ora_batch(
    gene_sets: Mapping[str, Iterable[str]],
    term2gene: Mapping[str, Iterable[str]],
    universe: Iterable[str],
    *,
    term2name: Mapping[str, str] | None = None,
    rscript: str | None = None,
) -> dict[str, dict[str, tuple[float, float, int]]]:
    """Run clusterProfiler ORA for MANY gene sets in ONE parallel R session.

    The performance workhorse. Writes every gene set (long-form ``name<TAB>gene``
    TSV) + the shared TERM2GENE + universe + optional TERM2NAME ONCE, then a
    single ``Rscript`` loads ``clusterProfiler`` once, parses the shared frames
    once, and ``parallel::mclapply``s ``enricher`` over all gene sets across
    ``min($MPD_MAX_CORES, detectCores() - 4)`` cores. This amortises R startup +
    library attach + TERM2GENE parse across the whole cell (a v2 reproduction
    cell = hundreds of miRNAs against one TERM2GENE), so a cell runs in seconds
    rather than the hours a per-miRNA subprocess would take.

    All gene sets share the SAME ``term2gene``, ``universe``, and (effective)
    background — exactly the v2-cell contract (one TERM2GENE per (org, variant,
    family); many miRNAs).

    Args:
        gene_sets: ``set_name -> iterable(gene)`` (e.g. ``miRNA_name -> targets``).
        term2gene: Shared ``term -> iterable(gene)`` membership.
        universe: Shared background gene identifiers.
        term2name: Optional shared ``term -> display name`` (keys the per-set
            results by the display name, matching v2's name-keyed matrix columns).
        rscript: Override the Rscript binary (else $MPD_RSCRIPT / PATH).

    Returns:
        ``{set_name: {term: (pvalue, p_adjust, overlap_count)}}``. A gene set
        for which enricher scored no term maps to an empty dict.

    Raises:
        RuntimeError: if the R subprocess fails (stderr is surfaced).
    """
    result, _provenance = _run_ora_batch_full(
        gene_sets, term2gene, universe, term2name=term2name, rscript=rscript
    )
    return result


def run_ora_batch_with_provenance(
    gene_sets: Mapping[str, Iterable[str]],
    term2gene: Mapping[str, Iterable[str]],
    universe: Iterable[str],
    *,
    term2name: Mapping[str, str] | None = None,
    rscript: str | None = None,
) -> tuple[dict[str, dict[str, tuple[float, float, int]]], dict[str, str]]:
    """Run the batched ORA and return BOTH results AND the run provenance.

    The PUBLIC provenance-carrying entrypoint (H6-Q6b): identical numerics to
    :func:`run_ora_batch` but also returns the R-written provenance JSON
    (clusterProfiler + R versions, mc_cores, n_gene_sets) probed inside the same
    interpreter that ran the analysis -- the authoritative version pin every
    production runner records in its manifest. Production runners and the
    controlled-reproduction runner call THIS instead of the private
    :func:`_run_ora_batch_full` (which remains the internal implementation).

    Args:
        gene_sets: ``set_name -> iterable(gene)`` (e.g. ``miRNA_name -> targets``).
        term2gene: Shared ``term -> iterable(gene)`` membership.
        universe: Shared background gene identifiers.
        term2name: Optional shared ``term -> display name``.
        rscript: Override the Rscript binary (else $MPD_RSCRIPT / PATH).

    Returns:
        ``({set_name: {term: (pvalue, p_adjust, overlap)}}, provenance)``.

    Raises:
        RuntimeError: if the R subprocess fails (stderr is surfaced).
    """
    return _run_ora_batch_full(
        gene_sets, term2gene, universe, term2name=term2name, rscript=rscript
    )


def _run_ora_batch_full(
    gene_sets: Mapping[str, Iterable[str]],
    term2gene: Mapping[str, Iterable[str]],
    universe: Iterable[str],
    *,
    term2name: Mapping[str, str] | None = None,
    rscript: str | None = None,
) -> tuple[dict[str, dict[str, tuple[float, float, int]]], dict[str, str]]:
    """Batched enricher run, returning ``(per_set_results, r_provenance)``.

    The single R subprocess writes the combined result table AND a provenance
    JSON (versions probed inside the same interpreter that ran the analysis).
    """
    rscript_bin = rscript or _resolve_rscript()
    mc_cores = _resolve_max_cores()
    universe_list = [str(g).strip() for g in universe if str(g).strip()]
    t2g_pairs = _term2gene_pairs(term2gene)

    # Long-form (name, gene) pairs, deterministic; blank genes/names dropped.
    geneset_pairs: list[tuple[str, str]] = []
    set_names: list[str] = []
    for name in gene_sets:
        sname = str(name).strip()
        if not sname:
            continue
        set_names.append(sname)
        for gene in sorted({str(g).strip() for g in gene_sets[name] if str(g).strip()}):
            geneset_pairs.append((sname, gene))

    with tempfile.TemporaryDirectory(prefix="mpd_cp_batch_") as tmp:
        tmp_path = Path(tmp)
        genesets_file = tmp_path / "gene_sets.tsv"
        t2g_file = tmp_path / "term2gene.tsv"
        universe_file = tmp_path / "universe.txt"
        out_file = tmp_path / "result.tsv"
        prov_file = tmp_path / "provenance.json"
        driver_file = tmp_path / "driver_batch.R"

        _write_pairs(genesets_file, geneset_pairs)
        _write_lines(universe_file, universe_list)
        _write_pairs(t2g_file, t2g_pairs)

        if term2name:
            t2n_file = tmp_path / "term2name.tsv"
            _write_pairs(t2n_file, sorted((str(t), str(n)) for t, n in term2name.items()))
            t2n_arg = str(t2n_file)
        else:
            t2n_arg = "NA"

        driver = _R_JSON_HELPER + _R_BATCH_DRIVER.replace("%MAX_GS_SIZE%", str(_MAX_GS_SIZE))
        driver_file.write_text(driver, encoding="utf-8")

        proc = subprocess.run(
            [
                rscript_bin,
                str(driver_file),
                str(genesets_file),
                str(t2g_file),
                str(universe_file),
                t2n_arg,
                str(out_file),
                str(prov_file),
                str(mc_cores),
                str(_CORE_HEADROOM),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"clusterProfiler batch Rscript failed (exit {proc.returncode}).\n"
                f"--- stderr ---\n{proc.stderr}\n--- stdout ---\n{proc.stdout}"
            )

        # When a TERM2NAME is supplied, key per-term by the display name
        # (enricher's Description); otherwise by the term id.
        result = _read_batch_result(
            out_file, set_names, key_by_name=bool(term2name)
        )
        provenance = _read_provenance_file(prov_file)
        _warn_on_partial_worker_errors(provenance)
        return result, provenance


def _warn_on_partial_worker_errors(provenance: Mapping[str, str]) -> None:
    """Surface partial per-set worker failures the R driver recorded.

    The R batch driver drops gene sets whose ``enricher`` call threw (a *partial*
    failure, distinct from the all-failed case which already raises). Those sets
    used to vanish silently, leaving holes in the matrix. The driver now records
    ``n_worker_errors`` in provenance; warn on it here so the hole is visible.
    """
    try:
        n_werr = int(provenance.get("n_worker_errors", "0") or "0")
    except (TypeError, ValueError):
        n_werr = 0
    if n_werr > 0:
        warnings.warn(
            f"clusterProfiler batch: {n_werr} gene set(s) failed in workers and were "
            f"dropped from the result (examples: "
            f"{provenance.get('worker_error_examples', '')}). The affected miRNAs are "
            f"absent from this matrix -- investigate before publishing.",
            stacklevel=3,
        )


def _read_provenance_file(prov_file: Path) -> dict[str, str]:
    """Read the R-written provenance JSON (versions from the run interpreter).

    L5: ``warnings.warn`` before returning ``{}`` on a missing / unparseable file
    so a lost version pin (clusterProfiler + R versions) does not vanish silently
    into an empty provenance dict.
    """
    if not prov_file.exists():
        warnings.warn(
            f"clusterProfiler provenance JSON not written ({prov_file}); the run's "
            f"clusterProfiler + R version pins are unavailable -- recording empty "
            f"provenance.",
            stacklevel=2,
        )
        return {}
    try:
        return json.loads(prov_file.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        warnings.warn(
            f"clusterProfiler provenance JSON unreadable ({prov_file}): {exc} -- "
            f"recording empty provenance (version pins lost).",
            stacklevel=2,
        )
        return {}


def _read_result(
    out_file: Path, *, key_by_name: bool = False
) -> dict[str, tuple[float, float, int]]:
    """Parse the R result TSV into ``{term: (pvalue, p_adjust, overlap)}``.

    Args:
        out_file: The R-written result TSV
            (``ID / Description / pvalue / p.adjust / Count``).
        key_by_name: If True, key the dict by the ``Description`` (display name);
            else by ``ID``. ``Description`` equals ``ID`` when no TERM2NAME was
            supplied, so this only matters when names were provided.
    """
    result: dict[str, tuple[float, float, int]] = {}
    if not out_file.exists():
        return result
    with out_file.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        key_col = "Description" if key_by_name else "ID"
        for row in reader:
            term = (row.get(key_col) or row.get("ID") or "").strip()
            if not term:
                continue
            result[term] = (
                float(row["pvalue"]),
                float(row["p.adjust"]),
                int(float(row["Count"])),
            )
    return result


def _read_batch_result(
    out_file: Path,
    set_names: Iterable[str],
    *,
    key_by_name: bool = False,
) -> dict[str, dict[str, tuple[float, float, int]]]:
    """Parse the combined batch TSV into ``{set_name: {term: (p, padj, overlap)}}``.

    The batch result TSV has a leading ``name`` column (the gene-set name) plus
    the per-term ``ID / Description / pvalue / p.adjust / Count``. Every input
    set name is pre-seeded with an empty dict so a set enricher scored no term
    for still appears in the result (callers can distinguish "no hits" from
    "missing").

    Args:
        out_file: The R-written combined result TSV.
        set_names: All gene-set names submitted (used to pre-seed empties).
        key_by_name: Key each term by ``Description`` (display name) if True,
            else by ``ID`` (matches :func:`_read_result`).
    """
    result: dict[str, dict[str, tuple[float, float, int]]] = {
        str(n): {} for n in set_names
    }
    if not out_file.exists():
        return result
    with out_file.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        key_col = "Description" if key_by_name else "ID"
        for row in reader:
            set_name = (row.get("name") or "").strip()
            if not set_name:
                continue
            term = (row.get(key_col) or row.get("ID") or "").strip()
            if not term:
                continue
            result.setdefault(set_name, {})[term] = (
                float(row["pvalue"]),
                float(row["p.adjust"]),
                int(float(row["Count"])),
            )
    return result


class ClusterProfilerEngine:
    """clusterProfiler engine adapter (implements :class:`EnrichmentEngine`).

    ``run`` requires an explicit TERM2GENE for the request's pathway source; for
    the controlled v2 reproduction the caller passes it via ``params["term2gene"]``
    (and optionally ``params["term2name"]``). Production GO TERM2GENE is
    materialised from ``org.*.eg.db`` by :mod:`mirpathdb_prod.pathways.go_genesets`
    (Reactome/WikiPathways/inferred-Reactome have their own materialisers); the
    production runners build that TERM2GENE and pass it to
    :func:`run_ora_batch_with_provenance` directly rather than through this generic
    adapter. Calling ``run`` without a TERM2GENE raises ``NotImplementedError`` to
    make the explicit-TERM2GENE contract clear.
    """

    name = ENGINE_CLUSTERPROFILER

    def run(self, request: EnrichmentRequest) -> EnrichmentResult:
        """Run GO/pathway over-representation via clusterProfiler.

        Args:
            request: The engine-agnostic request. ``request.params`` may carry
                ``term2gene`` (``term -> iterable(gene)``) and ``term2name``;
                ``request.background`` is the universe (required here — the
                controlled reproduction always supplies an explicit universe).

        Returns:
            An :class:`EnrichmentResult` with one row per scored term
            (``term`` / ``pvalue`` / ``p_adjust`` / ``overlap``) plus pinned
            ``provenance``.

        Raises:
            NotImplementedError: if no ``term2gene`` is supplied. This generic
                adapter does not materialise GO TERM2GENE itself; production
                callers build it via :mod:`mirpathdb_prod.pathways.go_genesets`
                and pass it in explicitly.
            ValueError: if no universe/background is supplied.
        """
        term2gene = request.params.get("term2gene")
        if term2gene is None:
            raise NotImplementedError(
                "ClusterProfilerEngine.run needs an explicit TERM2GENE in "
                "request.params['term2gene']. Production GO TERM2GENE is "
                "materialised by mirpathdb_prod.pathways.go_genesets and passed "
                "in by the runners; this generic adapter does not build it."
            )
        if request.background is None:
            raise ValueError(
                "ClusterProfilerEngine.run requires an explicit background/universe "
                "(request.background) for the over-representation test"
            )
        term2name = request.params.get("term2name")
        rscript_bin = request.params.get("rscript")
        if not isinstance(rscript_bin, str):
            rscript_bin = _resolve_rscript()

        raw, r_prov = _run_ora_full(
            request.genes,
            term2gene,  # type: ignore[arg-type]
            request.background,
            term2name=term2name if isinstance(term2name, Mapping) else None,
            rscript=rscript_bin,
        )
        rows: list[dict[str, object]] = [
            {"term": term, "pvalue": p, "p_adjust": padj, "overlap": overlap}
            for term, (p, padj, overlap) in sorted(raw.items())
        ]
        # Provenance from the SAME interpreter that ran the analysis (the
        # authoritative pin); fill the remaining pin slots from the request.
        provenance: dict[str, object] = dict.fromkeys(PROVENANCE_PINS)
        provenance["clusterprofiler_version"] = r_prov.get("clusterprofiler_version")
        provenance["r_version"] = r_prov.get("r_version")
        provenance["org_db_package"] = ORG_DB_PACKAGES.get(request.organism)
        return EnrichmentResult(
            engine=self.name,
            request=request,
            rows=rows,
            provenance=provenance,
        )


def make_term2name(terms: Sequence[str]) -> dict[str, str]:
    """Identity TERM2NAME for v2 term names (already display strings)."""
    return {t: t for t in terms}
