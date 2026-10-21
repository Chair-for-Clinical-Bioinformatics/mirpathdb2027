"""Production WikiPathways gene-set source (per-organism GMT -> TERM2GENE).

This is the **production** WikiPathways gene-set source for the uniform GO +
Reactome + WikiPathways enrichment build. WikiPathways publishes a *monthly data
release* of per-organism GMT files at ``https://data.wikipathways.org``; this
module downloads (or reuses a cached) ``wikipathways-<date>-gmt-<Species>.gmt``
and parses it into a single-family ``TERM2GENE`` map ``WPID -> set(gene SYMBOL)``
plus a ``TERM2NAME`` map ``WPID -> pathway name`` for clusterProfiler ``enricher``.

The WikiPathways GMT identifier line carries **Entrez gene IDs**, NOT symbols::

    <name>%WikiPathways_<date>%<WPID>%<Species>\t<url>\t<entrez1>\t<entrez2>\t...

so to keep the gene-ID space consistent with the GO source
(:mod:`mirpathdb_prod.pathways.go_genesets`, gene-SYMBOL space) and the target
sets, we map every Entrez member to SYMBOL **through the same ``org.*.eg.db``**
the GO sets came from (:func:`entrez_to_symbol`). Entrez IDs that do not resolve
are dropped (recorded in the per-organism census). The column-1 token also carries
the canonical WikiPathways ID (``WP100``) and the human-readable pathway name,
which become the TERM2GENE id and TERM2NAME respectively.

Organism coverage (the #1 thing to *check*, not assume): WikiPathways does NOT
ship every species. As of the 2026-06-10 release the model-organism collections
present are **Homo_sapiens / Mus_musculus / Rattus_norvegicus** (hsa/mmu/rno);
there is **no Macaca_mulatta** collection (the nearest primate is
*Pan troglodytes*). :data:`WP_SPECIES_FILE` therefore has no macaque entry and
:func:`species_file_for_taxid` raises for 9544 -- macaque gets NO WikiPathways,
matching the Reactome situation. If a future release adds Macaca_mulatta, add it
to :data:`WP_SPECIES_FILE` and it flows through unchanged.

ID space (the silent-failure risk): the WP member ids are Entrez; the org.db join
is the ONLY thing keeping WP in the GO/target SYMBOL space, so the join must use
the organism's own ``org.*.eg.db`` (org.Hs for hsa, org.Mm for mmu, org.Rn for
rno) -- never a cross-species map.

This module shells out to R for the Entrez->SYMBOL join (like
:mod:`mirpathdb_prod.pathways.go_genesets`); it does not import ``rpy2``. The R
binary is ``$MPD_RSCRIPT`` if set, else the first ``Rscript`` on ``PATH``. The GMT
download uses the stdlib (``urllib``) -- no third-party HTTP client.
"""

from __future__ import annotations

from mirpathdb_prod.pathways.label_normalization import (
    apply_name_repairs,
    assert_labels_clean,
    load_name_repairs,
    normalize_label,
)

import csv
import hashlib
import os
import shutil
import subprocess
import tempfile
import urllib.request
import warnings
from dataclasses import dataclass, field
from pathlib import Path

from mirpathdb_prod.pathways.go_genesets import ORG_DB_PACKAGES, org_db_for_taxid

# The single pathway family this source emits (matches the gmt.py token).
WIKIPATHWAYS_FAMILY = "WIKIPATHWAYS"

# WikiPathways monthly data release base + the release we pin for this build.
# The release date is recorded in provenance so the gene sets are auditable; bump
# WP_RELEASE to re-pin against a newer monthly release.
WP_DATA_BASE = "https://data.wikipathways.org"
WP_RELEASE = "20260610"

# Canonical taxid -> WikiPathways species token (the ``-gmt-<Species>.gmt`` suffix).
# ONLY species WikiPathways actually publishes are listed; macaque (9544) is
# deliberately ABSENT -- the release ships no Macaca_mulatta collection.
WP_SPECIES_FILE: dict[int, str] = {
    9606: "Homo_sapiens",
    10090: "Mus_musculus",
    10116: "Rattus_norvegicus",
}

RSCRIPT_ENV = "MPD_RSCRIPT"

# R driver: map a list of Entrez gene ids to SYMBOL through an org.*.eg.db package
# (the same annotation the GO sets came from), writing a two-column TSV. Identical
# shape to the GO/RefSeq drivers so the symbol space cannot drift.
_R_ENTREZ_DRIVER = r"""
suppressMessages(library(AnnotationDbi))
args      <- commandArgs(trailingOnly = TRUE)
pkg       <- args[[1]]
gene_file <- args[[2]]
out_file  <- args[[3]]
suppressMessages(library(pkg, character.only = TRUE))
db <- get(pkg)
ids <- readLines(gene_file)
ids <- unique(trimws(ids))
ids <- ids[nzchar(ids)]
res <- suppressMessages(AnnotationDbi::select(
  db, keys = ids, columns = "SYMBOL", keytype = "ENTREZID"))
res <- res[!is.na(res$SYMBOL), ]
res <- res[!duplicated(res$ENTREZID), ]
write.table(res[, c("ENTREZID", "SYMBOL")], file = out_file, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
prov <- sprintf('{"org_db_package": "%s", "org_db_version": "%s", "r_version": "%s"}',
                pkg, as.character(packageVersion(pkg)),
                paste(R.version$major, R.version$minor, sep = "."))
writeLines(prov, paste0(out_file, ".prov"))
"""


@dataclass
class WikiPathwaysGeneSets:
    """Materialised WikiPathways gene sets for one organism (single family).

    Attributes:
        taxid: Canonical NCBI taxid the gene sets were built for.
        species: The WikiPathways species token (e.g. ``Homo_sapiens``).
        release: The WikiPathways monthly release pinned (e.g. ``20260610``).
        org_db_package: The ``org.*.eg.db`` used for the Entrez->SYMBOL join.
        term2gene: ``{WPID: set(symbol)}`` (the single WIKIPATHWAYS family).
        term2name: ``{WPID: pathway_name}``.
        provenance: Version/release pins (release, org.db package + version, R).
    """

    taxid: int
    species: str
    release: str
    org_db_package: str
    term2gene: dict[str, set[str]]
    term2name: dict[str, str]
    provenance: dict[str, str] = field(default_factory=dict)

    def family(self, family: str = WIKIPATHWAYS_FAMILY) -> dict[str, set[str]]:
        """Return the ``WPID -> set(symbol)`` map (only WIKIPATHWAYS exists)."""
        if family != WIKIPATHWAYS_FAMILY:
            raise KeyError(
                f"no WikiPathways family {family!r}; only {WIKIPATHWAYS_FAMILY!r}"
            )
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


def _file_sha256(path: Path, *, chunk: int = 1 << 20) -> str:
    """Stream a sha256 of the (possibly large) GMT for snapshot pinning (H4)."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def species_file_for_taxid(taxid: int) -> str:
    """Return the WikiPathways species token for a canonical taxid.

    Raises:
        ValueError: if WikiPathways publishes no collection for ``taxid`` (e.g.
            macaque 9544) -- the caller must treat that organism as WP = none.
    """
    try:
        return WP_SPECIES_FILE[taxid]
    except KeyError as exc:
        raise ValueError(
            f"WikiPathways publishes no collection for taxid {taxid}; "
            f"available: {sorted(WP_SPECIES_FILE)} "
            f"(no Macaca_mulatta -> macaque gets no WikiPathways)"
        ) from exc


def gmt_filename(taxid: int, *, release: str = WP_RELEASE) -> str:
    """Return the WikiPathways GMT filename for ``taxid`` at ``release``."""
    return f"wikipathways-{release}-gmt-{species_file_for_taxid(taxid)}.gmt"


def gmt_url(taxid: int, *, release: str = WP_RELEASE) -> str:
    """Return the WikiPathways GMT download URL for ``taxid`` (H4).

    Pinned to the DATED archive path ``{WP_DATA_BASE}/{release}/gmt/...`` (NOT the
    floating ``/current/gmt/``) so a rebuild fetches the SAME monthly snapshot and
    the WikiPathways matrices stay rebuildable after the upstream monthly roll.
    ``release`` is the WikiPathways archive date token (e.g. ``20260610``), which
    matches the data-server archive layout ``https://data.wikipathways.org/<date>/
    gmt/``. Cached GMTs are reused if present (see :func:`fetch_gmt`), so this URL
    is only hit on a fresh download; a missing dated directory yields a loud HTTP
    error (not a silent wrong snapshot). Snapshot integrity is additionally pinned
    by the GMT's recorded ``sha256`` + ``size_bytes`` in
    :func:`build_wikipathways_genesets`'s provenance.
    """
    return f"{WP_DATA_BASE}/{release}/gmt/{gmt_filename(taxid, release=release)}"


def fetch_gmt(
    taxid: int, *, cache_dir: str | Path, release: str = WP_RELEASE, force: bool = False
) -> Path:
    """Download (or reuse a cached) WikiPathways GMT for ``taxid``.

    The GMT lands at ``<cache_dir>/wikipathways/<filename>``. A non-empty cached
    file is reused unless ``force``. Uses the stdlib ``urllib`` (no third-party
    HTTP dependency); the WikiPathways data server serves the *current* monthly
    release under ``/current/gmt/`` and the date is recorded in the filename. The
    GMT's sha256 + size are recorded in provenance (H4) so the snapshot is pinned
    even though the URL still points at the floating ``/current/`` (see
    :func:`gmt_url`'s TODO).

    Raises:
        ValueError: if ``taxid`` has no WikiPathways collection.
        urllib.error.URLError / OSError: on a failed download.
    """
    dest_dir = Path(cache_dir) / "wikipathways"
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / gmt_filename(taxid, release=release)
    if dest.exists() and dest.stat().st_size > 0 and not force:
        return dest
    url = gmt_url(taxid, release=release)
    with urllib.request.urlopen(url, timeout=300) as resp:  # noqa: S310 (https only)
        data = resp.read()
    dest.write_bytes(data)
    return dest


def parse_wp_gmt(path: str | Path) -> tuple[dict[str, set[str]], dict[str, str]]:
    """Parse a WikiPathways GMT into ``(WPID -> set(entrez), WPID -> name)``.

    The WikiPathways column-1 token is ``<name>%WikiPathways_<date>%<WPID>%<Species>``;
    we split on ``%`` to recover the pathway NAME (field 0) and the canonical WPID
    (field 2 -- e.g. ``WP100``). Members (columns 3..N) are **Entrez gene ids**
    (still unmapped here; :func:`entrez_to_symbol` maps them to SYMBOL).

    Raises:
        FileNotFoundError: if ``path`` does not exist.
    """
    gmt_path = Path(path)
    if not gmt_path.exists():
        raise FileNotFoundError(f"WikiPathways GMT not found: {gmt_path}")
    term2entrez: dict[str, set[str]] = {}
    term2name: dict[str, str] = {}
    with gmt_path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) < 3:  # need name + url + >=1 member to be useful
                continue
            token = cols[0]
            parts = token.split("%")
            # WPID is the %-field that starts with "WP"; fall back to field 2.
            wpid = ""
            for p in parts:
                if p.startswith("WP") and p[2:].isdigit():
                    wpid = p
                    break
            if not wpid and len(parts) >= 3:
                wpid = parts[2].strip()
            name = parts[0].strip() if parts else token.strip()
            if not wpid:
                continue
            members = {c.strip() for c in cols[2:] if c.strip()}
            if not members:
                continue
            term2entrez.setdefault(wpid, set()).update(members)
            term2name.setdefault(wpid, name)
    return term2entrez, term2name


def entrez_to_symbol(
    entrez_ids: set[str], *, org_db_package: str, rscript: str | None = None
) -> tuple[dict[str, str], dict[str, str]]:
    """Map Entrez gene ids to SYMBOL via an ``org.*.eg.db`` package.

    Shells out to R once (the same annotation the GO sets come from) so the
    WikiPathways symbol space cannot drift from the GO TERM2GENE space. Returns
    ``({entrez: symbol}, provenance)``; ids with no symbol are simply absent.

    Raises:
        RuntimeError: if the R subprocess fails.
    """
    ids = sorted({e.strip() for e in entrez_ids if e and e.strip()})
    if not ids:
        return {}, {}
    rscript_bin = _resolve_rscript(rscript)
    with tempfile.TemporaryDirectory(prefix="mpd_wp_entrez_") as tmp:
        tmp_path = Path(tmp)
        gene_file = tmp_path / "entrez.txt"
        out_file = tmp_path / "entrez_symbol.tsv"
        driver_file = tmp_path / "driver_entrez.R"
        gene_file.write_text("\n".join(ids) + "\n", encoding="utf-8")
        driver_file.write_text(_R_ENTREZ_DRIVER, encoding="utf-8")
        proc = subprocess.run(
            [rscript_bin, str(driver_file), org_db_package, str(gene_file), str(out_file)],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"WikiPathways Entrez->symbol Rscript failed for {org_db_package} "
                f"(exit {proc.returncode}).\n--- stderr ---\n{proc.stderr}"
            )
        mapping: dict[str, str] = {}
        if out_file.exists():
            with out_file.open("r", encoding="utf-8", newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    eid = (row.get("ENTREZID") or "").strip()
                    sym = (row.get("SYMBOL") or "").strip()
                    if eid and sym:
                        mapping[eid] = sym
        provenance: dict[str, str] = {}
        prov_file = Path(str(out_file) + ".prov")
        if prov_file.exists():
            import json

            try:
                provenance = json.loads(prov_file.read_text(encoding="utf-8"))
            except (OSError, ValueError) as exc:
                # L5: surface a lost org.db/R version pin instead of silently {}.
                warnings.warn(
                    f"WikiPathways Entrez->symbol provenance unreadable "
                    f"({prov_file}): {exc} -- recording empty provenance "
                    f"(org.db/R version pins lost).",
                    stacklevel=2,
                )
                provenance = {}
        else:
            # L5: the R driver should always write the .prov sidecar; warn if not.
            warnings.warn(
                f"WikiPathways Entrez->symbol provenance sidecar not written "
                f"({prov_file}); org.db/R version pins unavailable.",
                stacklevel=2,
            )
        return mapping, provenance


def _map_term2gene_to_symbols(
    term2entrez: dict[str, set[str]], entrez_symbol: dict[str, str]
) -> tuple[dict[str, set[str]], int, int]:
    """Re-key a ``WPID -> set(entrez)`` map to ``WPID -> set(symbol)``.

    Returns ``(term2gene, n_mapped_pairs, n_unmapped_pairs)``. Entrez members
    with no symbol are dropped; pathways left empty after mapping are dropped.
    """
    out: dict[str, set[str]] = {}
    mapped = 0
    unmapped = 0
    for wpid, entrez in term2entrez.items():
        syms: set[str] = set()
        for e in entrez:
            sym = entrez_symbol.get(e)
            if sym:
                syms.add(sym)
                mapped += 1
            else:
                unmapped += 1
        if syms:
            out[wpid] = syms
    return out, mapped, unmapped


def build_wikipathways_genesets(
    taxid: int,
    *,
    cache_dir: str | Path,
    release: str = WP_RELEASE,
    rscript: str | None = None,
    force_download: bool = False,
) -> WikiPathwaysGeneSets:
    """Materialise WikiPathways gene sets (SYMBOL space) for ``taxid``.

    Downloads (or reuses cached) the organism's WikiPathways GMT, parses the
    Entrez-member pathways, and maps every member to SYMBOL through the
    organism's ``org.*.eg.db`` so the result is in the same gene-ID space as the
    GO TERM2GENE and the target sets.

    Args:
        taxid: Canonical NCBI taxid (must be in :data:`WP_SPECIES_FILE`).
        cache_dir: Directory for the GMT cache (``<cache_dir>/wikipathways``).
        release: WikiPathways monthly release token (default :data:`WP_RELEASE`).
        rscript: Override the Rscript binary (else $MPD_RSCRIPT / PATH).
        force_download: Re-download even if a cached GMT exists.

    Returns:
        A :class:`WikiPathwaysGeneSets` with the symbol-space TERM2GENE/TERM2NAME
        and release/org.db/R version pins in ``provenance``.

    Raises:
        ValueError: if ``taxid`` has no WikiPathways collection (e.g. macaque).
        RuntimeError: if the Entrez->symbol R subprocess fails.
    """
    species = species_file_for_taxid(taxid)
    pkg = org_db_for_taxid(taxid)
    gmt_path = fetch_gmt(taxid, cache_dir=cache_dir, release=release, force=force_download)
    term2entrez, term2name = parse_wp_gmt(gmt_path)

    all_entrez: set[str] = set()
    for members in term2entrez.values():
        all_entrez |= members
    entrez_symbol, r_prov = entrez_to_symbol(all_entrez, org_db_package=pkg, rscript=rscript)
    term2gene, mapped, unmapped = _map_term2gene_to_symbols(term2entrez, entrez_symbol)
    # Keep only names for pathways that survived symbol mapping.
    term2name = {wpid: term2name.get(wpid, wpid) for wpid in term2gene}
    # DATA-001: the upstream GMT export mangles HTML entities in titles (the
    # dated archive itself ships "Parkinson 39 s disease pathway"; full trace in
    # docs/DATA-001_wikipathways_label_root_cause.md). Normalize every display
    # name, repair residue-flagged ids from the same release's GPML titles via
    # the checked-in map, and refuse to build if any residue survives.
    term2name = {wpid: normalize_label(name) for wpid, name in term2name.items()}
    term2name, n_label_repairs = apply_name_repairs(
        term2name, load_name_repairs(release=release)
    )
    assert_labels_clean(
        term2name.items(), context=f"WikiPathways {species} {release} TERM2NAME"
    )

    provenance: dict[str, str] = {
        "wikipathways_release": release,
        "wikipathways_species": species,
        "gmt_filename": gmt_path.name,
        # H4: pin the exact dated-archive URL + the GMT's sha256 so the snapshot is
        # verifiable on rebuild even if the upstream dated directory is later pruned.
        "gmt_url": gmt_url(taxid, release=release),
        "gmt_sha256": _file_sha256(gmt_path),
        "gmt_size_bytes": str(gmt_path.stat().st_size),
        "n_pathways_in_gmt": str(len(term2entrez)),
        "n_pathways_mapped": str(len(term2gene)),
        "n_entrez_pairs_mapped": str(mapped),
        "n_entrez_pairs_unmapped": str(unmapped),
        "n_label_repairs": str(n_label_repairs),
        **r_prov,
    }
    return WikiPathwaysGeneSets(
        taxid=taxid,
        species=species,
        release=release,
        org_db_package=pkg,
        term2gene=term2gene,
        term2name=term2name,
        provenance=provenance,
    )


def write_wikipathways_genesets_tsv(
    gene_sets: WikiPathwaysGeneSets, out_path: str | Path
) -> Path:
    """Persist a materialised WikiPathways TERM2GENE as a long-form TSV (cache/audit).

    Columns: ``family<TAB>wp_id<TAB>wp_name<TAB>symbol`` (one row per pair),
    deterministic order. :func:`read_wikipathways_genesets_tsv` reloads it.
    """
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", quoting=csv.QUOTE_NONE, escapechar="\\")
        writer.writerow(["family", "wp_id", "wp_name", "symbol"])
        for wpid in sorted(gene_sets.term2gene):
            name = gene_sets.term2name.get(wpid, wpid)
            for symbol in sorted(gene_sets.term2gene[wpid]):
                writer.writerow([WIKIPATHWAYS_FAMILY, wpid, name, symbol])
    return out


def read_wikipathways_genesets_tsv(
    path: str | Path, *, taxid: int, species: str, release: str, org_db_package: str
) -> WikiPathwaysGeneSets:
    """Reload a cached WikiPathways TERM2GENE TSV (written by the writer above)."""
    term2gene: dict[str, set[str]] = {}
    term2name: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            wpid = (row.get("wp_id") or "").strip()
            symbol = (row.get("symbol") or "").strip()
            if not wpid or not symbol:
                continue
            term2gene.setdefault(wpid, set()).add(symbol)
            name = (row.get("wp_name") or "").strip()
            if name:
                term2name.setdefault(wpid, name)
    return WikiPathwaysGeneSets(
        taxid=taxid,
        species=species,
        release=release,
        org_db_package=org_db_package,
        term2gene=term2gene,
        term2name=term2name,
    )


__all__ = [
    "ORG_DB_PACKAGES",
    "WIKIPATHWAYS_FAMILY",
    "WP_DATA_BASE",
    "WP_RELEASE",
    "WP_SPECIES_FILE",
    "WikiPathwaysGeneSets",
    "build_wikipathways_genesets",
    "entrez_to_symbol",
    "fetch_gmt",
    "gmt_filename",
    "gmt_url",
    "parse_wp_gmt",
    "read_wikipathways_genesets_tsv",
    "species_file_for_taxid",
    "write_wikipathways_genesets_tsv",
]
