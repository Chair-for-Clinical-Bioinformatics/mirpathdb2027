"""Single-pass, multi-cell TERM2GENE reconstruction from a GeneTrail2 tarball.

:func:`mirpathdb_prod.pathways.v2_genesets.load_v2_term2gene` streams the whole
tarball **once per (variant, family)**. That is fine for the 136 MB mmu tarball
but pathological for the 72 GB hsa tarball (gzip is not seekable, so each call
re-decompresses the entire archive — 10 cells would mean ten full passes).

This module reconstructs **all requested (variant, family) cells in ONE pass**.
It optionally pipes the archive through ``unpigz`` (parallel gunzip) into a
``tarfile`` stream, which on this host decompresses far faster than Python's
single-threaded :mod:`gzip`. The per-cell result objects are byte-for-byte the
same :class:`~mirpathdb_prod.pathways.v2_genesets.V2TermGenes` the single-cell
loader produces — so downstream code is identical.

Read-only; no network, no R.
"""

from __future__ import annotations

import shutil
import subprocess
import tarfile
from collections.abc import Iterable
from pathlib import Path

from mirpathdb_prod.pathways.gmt import KNOWN_FAMILIES
from mirpathdb_prod.pathways.v2_genesets import (
    V2TermGenes,
    _member_relpath_parts,
    _parse_result_text,
)


def _open_tar_stream(tar_path: Path) -> tuple[tarfile.TarFile, subprocess.Popen | None]:
    """Open a streaming TarFile, preferring ``unpigz`` for fast decompression.

    Returns ``(tarfile, popen_or_None)``. When ``unpigz`` is available the tar is
    opened in non-seeking stream mode (``r|``) over the decompressor's stdout;
    otherwise it falls back to Python's gzip (``r|gz``) over the raw file.
    """
    unpigz = shutil.which("unpigz") or shutil.which("pigz")
    if unpigz is not None:
        args = [unpigz, "-c", str(tar_path)] if unpigz.endswith("unpigz") else [
            unpigz, "-dc", str(tar_path)
        ]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        assert proc.stdout is not None
        tf = tarfile.open(fileobj=proc.stdout, mode="r|")
        return tf, proc
    tf = tarfile.open(tar_path, mode="r:gz")
    return tf, None


def load_v2_term2gene_multi(
    genetrail2_tarball: str | Path,
    *,
    organism: str,
    variants: Iterable[str],
    families: Iterable[str],
) -> dict[tuple[str, str], V2TermGenes]:
    """Reconstruct TERM2GENE for many (variant, family) cells in ONE pass.

    Args:
        genetrail2_tarball: Path to ``miRPathDB2_<org>_genetrail2_results.tar.gz``.
        organism: ``"hsa"`` or ``"mmu"`` (selects the in-tar directory marker).
        variants: GeneTrail2 subtree variant tokens (``miRTarBase_strong`` /
            ``miRTarBase_all``).
        families: Families to collect (subset of
            :data:`mirpathdb_prod.pathways.gmt.KNOWN_FAMILIES`).

    Returns:
        ``{(variant, family): V2TermGenes}`` for every requested combination
        (a combination with no members in the tarball maps to an empty
        :class:`~mirpathdb_prod.pathways.v2_genesets.V2TermGenes`).

    Raises:
        FileNotFoundError: if the tarball is missing.
        ValueError: if an unknown family is requested.
    """
    tar_path = Path(genetrail2_tarball)
    if not tar_path.exists():
        raise FileNotFoundError(f"GeneTrail2 tarball not found: {tar_path}")

    want_variants = set(variants)
    want_families = set(families)
    bad = want_families - set(KNOWN_FAMILIES)
    if bad:
        raise ValueError(f"unknown families requested: {sorted(bad)}")

    marker = f"miRPathDB2_{organism}_genetrail2_results"
    out: dict[tuple[str, str], V2TermGenes] = {
        (v, f): V2TermGenes() for v in want_variants for f in want_families
    }

    tf, proc = _open_tar_stream(tar_path)
    try:
        for member in tf:
            if not member.isfile():
                continue
            parsed = _member_relpath_parts(member.name, marker)
            if parsed is None:
                continue
            m_variant, mimat, m_family = parsed
            if m_variant not in want_variants or m_family not in want_families:
                continue
            extracted = tf.extractfile(member)
            if extracted is None:
                continue
            # Stream-mode (``r|``) extractfile objects are not seekable, so read
            # raw bytes and decode rather than wrapping in a TextIOWrapper.
            body = extracted.read().decode("utf-8", errors="replace")
            cell = out[(m_variant, m_family)]
            mirna_genes = cell.per_mirna_genes.setdefault(mimat, set())
            for term_name, genes in _parse_result_text(body):
                bucket = cell.term2gene.setdefault(term_name, set())
                bucket.update(genes)
                cell.term2name.setdefault(term_name, term_name)
                cell.universe.update(genes)
                mirna_genes.update(genes)
    finally:
        tf.close()
        if proc is not None:
            if proc.stdout is not None:
                proc.stdout.close()
            proc.wait()
    return out
