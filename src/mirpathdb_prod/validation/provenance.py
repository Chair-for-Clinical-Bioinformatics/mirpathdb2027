"""Provenance validation.

Verifies the upstream-source manifest (sizes / checksums of the acquired source
files) so the production build is reproducible and auditable. The manifest lives
at ``provenance/manifest.yaml``; this validator confirms the acquired artefacts on
disk match the recorded pins.

Two modes, both driven by :func:`validate_manifest`:

* **verify** — for every ``raw_files`` entry that carries a ``sha256`` pin, recompute
  the file's sha256 and require it to match (a mismatch or a missing pinned file is
  a hard failure).
* **record** — for entries marked ``checksum: on-demand`` (no recorded sha256), the
  on-disk sha256 + size are computed and returned under ``skipped`` so they can be
  transcribed into the manifest, but they are NOT failures.

:func:`validate_manifest` is the full audit (verify + record modes over every
declared source). :func:`verify_cached_sources` is the lighter **start-up gate**
the production runners call: it verifies only the sha256-pinned cached sources
that are *present* under the build cache and raises on a MISMATCH (swapped / stale
/ corrupt bytes), leaving an ABSENT file to each runner's own ``_cached()``
existence guard — so a portable build on a cleared cache is not spuriously failed
while tampering is still caught.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import yaml

# Repo root = three levels up (src/mirpathdb_prod/validation/provenance.py).
REPO_ROOT = Path(__file__).resolve().parents[3]


def file_sha256(path: Path, *, chunk: int = 1 << 20) -> str:
    """Stream a sha256 hex digest of a (possibly multi-GB) file."""
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _iter_raw_files(doc: dict[str, Any]):
    """Yield ``(source_name, entry)`` for every ``raw_files`` entry in the manifest.

    Walks the top-level ``sources`` mapping (the only place ``raw_files`` lists
    live); each source entry may carry a ``raw_files`` list of
    ``{path, bytes, sha256|checksum, note?}`` dicts.
    """
    sources = doc.get("sources", {}) or {}
    for name, entry in sources.items():
        if not isinstance(entry, dict):
            continue
        for raw in entry.get("raw_files", []) or []:
            if isinstance(raw, dict) and raw.get("path"):
                yield name, raw


def _resolve(path_str: str, data_dir: Path) -> Path:
    """Absolute paths used as-is; relative paths resolved against ``data_dir``."""
    p = Path(path_str).expanduser()
    return p if p.is_absolute() else (data_dir / p).resolve()


def validate_manifest(
    manifest_path: Path, data_dir: Path | None = None, *, require_present: bool = True
) -> dict:
    """Validate acquired source files against the provenance manifest.

    Args:
        manifest_path: Path to ``provenance/manifest.yaml``.
        data_dir: Base for resolving relative ``raw_files`` paths (default: repo
            root). Absolute paths
            (e.g. ``/tmp/mpd-data/...``) are used verbatim.
        require_present: When True (production default), a pinned file that is
            missing or whose sha256 mismatches is a hard failure. When False, a
            missing file is recorded as a failure but unpinned-missing is tolerated.

    Returns:
        ``{"ok": bool, "checked": [...], "failures": [...], "skipped": [...]}``
        where each list element is a dict with ``source``, ``path`` and (where
        computed) ``sha256``/``bytes``. ``ok`` is True iff there are no failures.
    """
    manifest_path = Path(manifest_path)
    base = Path(data_dir) if data_dir is not None else REPO_ROOT

    with manifest_path.open("r", encoding="utf-8") as fh:
        doc: dict[str, Any] = yaml.safe_load(fh) or {}

    checked: list[dict] = []
    failures: list[dict] = []
    skipped: list[dict] = []

    for source, raw in _iter_raw_files(doc):
        path = _resolve(str(raw["path"]), base)
        pinned = raw.get("sha256")
        rec = {"source": source, "path": str(path)}

        if not path.exists():
            # A directory entry (e.g. the kegg/ cache dir) is treated as present.
            if path.is_dir():
                skipped.append({**rec, "note": "directory entry (not hashed)"})
                continue
            rec["error"] = "missing"
            (failures if (require_present or pinned) else skipped).append(rec)
            continue

        if path.is_dir():
            skipped.append({**rec, "note": "directory entry (not hashed)"})
            continue

        digest = file_sha256(path)
        size = path.stat().st_size

        if pinned:
            if digest == pinned:
                checked.append({**rec, "sha256": digest, "bytes": size})
            else:
                failures.append(
                    {**rec, "error": "sha256_mismatch", "expected": pinned, "got": digest}
                )
        else:
            # checksum: on-demand — record the computed identity for transcription.
            skipped.append(
                {**rec, "sha256": digest, "bytes": size, "note": "unpinned (record-mode)"}
            )

    return {
        "ok": not failures,
        "checked": checked,
        "failures": failures,
        "skipped": skipped,
    }


def pinned_sha256_index(manifest_path: Path) -> dict[str, str]:
    """Build a ``{cache-relative key: sha256}`` index of the manifest's pinned sources.

    Keyed by each pin's position RELATIVE to its cache root (not by the absolute
    ``/tmp/mpd-data`` path recorded in the manifest), so a cached source stays
    locatable under the config-driven ``source_cache_dir`` on any host:

    * a pin sitting directly in a cache root keeps its bare ``basename``;
    * a pin nested ONE level deeper -- i.e. whose parent directory itself lives in
      a directory that holds other pins -- is keyed by its last TWO components,
      ``"<subdir>/<basename>"`` (e.g. ``wikipathways/wikipathways-...-gmt-Homo_sapiens.gmt``).

    Keying every pin by basename alone (the previous behaviour) made the three
    nested WikiPathways GMTs resolve to ``<cache>/<file>.gmt``, a path that never
    exists, so they were silently reported ABSENT and their sha256 pins were never
    actually verified. Only entries carrying a real ``sha256`` pin are included
    (``checksum: on-demand`` entries are skipped).
    """
    with Path(manifest_path).open("r", encoding="utf-8") as fh:
        doc: dict[str, Any] = yaml.safe_load(fh) or {}
    pairs: list[tuple[Path, str]] = []
    for _source, raw in _iter_raw_files(doc):
        pinned = raw.get("sha256")
        if pinned:
            pairs.append((Path(str(raw["path"])), str(pinned)))

    # Cache roots = every directory that directly holds a pinned file. A pin whose
    # parent's PARENT is also such a root sits one level down inside the cache and
    # needs the extra component to be addressable under source_cache_dir.
    roots = {p.parent for p, _pinned in pairs}
    index: dict[str, str] = {}
    for p, pinned in pairs:
        key = f"{p.parent.name}/{p.name}" if p.parent.parent in roots else p.name
        index[key] = pinned
    return index


def verify_cached_sources(
    cache_dir: Path, manifest_path: Path, *, names: list[str] | None = None
) -> dict:
    """Start-up gate: verify PRESENT cached sources against their manifest sha256.

    For every sha256-pinned source in the manifest that EXISTS under ``cache_dir``,
    recompute its sha256 and raise :class:`RuntimeError` on a mismatch (swapped /
    stale / corrupt bytes) — refusing to build published artefacts on tampered
    inputs. A pinned source that is ABSENT is deliberately NOT failed here (that is
    the runner's own ``_cached()`` existence guard's job), so a portable build on a
    host with a different / cleared cache is not spuriously blocked. Sources present
    but unpinned are ignored.

    Args:
        cache_dir: The build cache (``paths.source_cache_dir``); pins are matched
            against files under here by their cache-relative key (see
            :func:`pinned_sha256_index`) -- a bare basename, or ``<subdir>/<basename>``
            for a source nested one level inside the cache.
        manifest_path: ``provenance/manifest.yaml``.
        names: Optional restriction to specific basenames (else all pinned sources);
            matched against each key's basename so callers stay basename-addressable
            regardless of where the source sits inside the cache.

    Returns:
        ``{"verified": [...], "absent": [...], "cache_dir": str}``.

    Raises:
        RuntimeError: if any present pinned source's sha256 does not match.
    """
    cache_dir = Path(cache_dir)
    index = pinned_sha256_index(manifest_path)
    if names is not None:
        wanted = {Path(n).name for n in names}
        index = {k: s for k, s in index.items() if Path(k).name in wanted}

    verified: list[str] = []
    absent: list[str] = []
    mismatches: list[dict] = []
    for key, pinned in sorted(index.items()):
        p = cache_dir / key
        if not p.exists():
            absent.append(key)
            continue
        digest = file_sha256(p)
        if digest == pinned:
            verified.append(key)
        else:
            mismatches.append({"file": key, "expected": pinned, "got": digest})

    if mismatches:
        detail = "; ".join(
            f"{m['file']} (pinned {m['expected'][:12]}.., got {m['got'][:12]}..)"
            for m in mismatches
        )
        raise RuntimeError(
            f"provenance gate: {len(mismatches)} cached source(s) under {cache_dir} do "
            f"NOT match their pinned sha256 -- refusing to build on swapped/stale/"
            f"corrupt inputs. Re-fetch the correct source or update the pin in "
            f"provenance/manifest.yaml. Mismatches: {detail}"
        )
    return {"verified": verified, "absent": absent, "cache_dir": str(cache_dir)}
