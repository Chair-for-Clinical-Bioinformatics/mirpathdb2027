"""Pathway-source registry — reads ``config/pathway_sources.yaml``.

Loads the machine-readable pathway-source registry into typed
:class:`PathwaySource` records and exposes role/provenance/version/license +
redistribution queries over them. This is the seam the (later) ingest +
publish-bundle steps consult to answer: *which sources ship in the published DB,
which are license-gated and user-only, and which are still under evaluation?*

No pathway data is fetched or ingested here; this only reads the registry config.
Path resolution mirrors :mod:`mirpathdb_prod.config` (repo-root relative) so no
absolute user paths are hard-coded.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

# Repo root = three levels up (src/mirpathdb_prod/pathways/registry.py -> ... -> root).
REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REGISTRY_PATH = REPO_ROOT / "config" / "pathway_sources.yaml"

# Valid roles (mirror config/pathway_sources.yaml and docs/pathway_sources.md).
ROLE_PRIMARY = "primary"
ROLE_OPTIONAL_LICENSED = "optional_licensed"
ROLE_EVALUATE = "evaluate"
VALID_ROLES: frozenset[str] = frozenset(
    {ROLE_PRIMARY, ROLE_OPTIONAL_LICENSED, ROLE_EVALUATE}
)

# Redistribution flag values.
REDIST_ALLOWED = "allowed"
REDIST_BLOCKED = "blocked"


@dataclass(frozen=True)
class PathwaySource:
    """One pathway-source registry entry.

    Attributes:
        id: Stable registry id (e.g. "go", "kegg").
        name: Human-readable name.
        role: One of :data:`VALID_ROLES`.
        redistribution: Whether the source's data may be redistributed in a
            public bundle (``"allowed"`` / ``"blocked"``).
        license: License identifier/string.
        provenance: source / version / retrieval mapping (versions may be "TBC"
            until pinned at build time).
        description: Free-text role description.
        notes: Free-text notes (gates, Tuesday flags, ...).
    """

    id: str
    name: str
    role: str
    redistribution: str
    license: str
    provenance: dict[str, Any]
    description: str = ""
    notes: str = ""

    @property
    def is_primary(self) -> bool:
        """True if this is a primary (default, shipped) source."""
        return self.role == ROLE_PRIMARY

    @property
    def ships_in_published_db(self) -> bool:
        """True if this source may ship in the published DB bundle.

        A source ships only if its role is not the license-gated optional role
        AND its data is redistributable. KEGG (optional_licensed + blocked) is
        excluded; GO (primary + allowed) ships.
        """
        return self.role != ROLE_OPTIONAL_LICENSED and self.redistribution == REDIST_ALLOWED


@dataclass(frozen=True)
class PathwaySourceRegistry:
    """The parsed pathway-source registry."""

    sources: tuple[PathwaySource, ...]
    source_path: Path | None = None

    def get(self, source_id: str) -> PathwaySource:
        """Return the source with ``source_id`` (raises KeyError if absent)."""
        for s in self.sources:
            if s.id == source_id:
                return s
        raise KeyError(f"no pathway source with id {source_id!r}")

    def by_role(self, role: str) -> tuple[PathwaySource, ...]:
        """Return all sources with the given role."""
        return tuple(s for s in self.sources if s.role == role)

    def primary(self) -> PathwaySource:
        """Return the single primary pathway source (raises if not exactly one)."""
        prim = self.by_role(ROLE_PRIMARY)
        if len(prim) != 1:
            raise ValueError(
                f"expected exactly one primary pathway source, found {len(prim)}: "
                f"{[s.id for s in prim]}"
            )
        return prim[0]

    def publishable(self) -> tuple[PathwaySource, ...]:
        """Return sources that may ship in the published DB (redistribution gate)."""
        return tuple(s for s in self.sources if s.ships_in_published_db)


def load_pathway_sources(path: str | Path | None = None) -> PathwaySourceRegistry:
    """Load ``config/pathway_sources.yaml`` into a :class:`PathwaySourceRegistry`.

    Args:
        path: Explicit path to the registry YAML; defaults to
            ``<repo_root>/config/pathway_sources.yaml``.

    Returns:
        A frozen :class:`PathwaySourceRegistry`.

    Raises:
        FileNotFoundError: if the registry file is missing.
        ValueError: if an entry has no ``id`` or an unknown ``role``.
    """
    registry_path = Path(path).expanduser().resolve() if path else DEFAULT_REGISTRY_PATH
    if not registry_path.exists():
        raise FileNotFoundError(f"pathway-source registry not found: {registry_path}")

    with registry_path.open("r", encoding="utf-8") as fh:
        doc: dict[str, Any] = yaml.safe_load(fh) or {}

    raw_sources = doc.get("pathway_sources", []) or []
    sources: list[PathwaySource] = []
    for entry in raw_sources:
        source_id = entry.get("id")
        if not source_id:
            raise ValueError(f"pathway-source entry without an id: {entry!r}")
        role = entry.get("role", "")
        if role not in VALID_ROLES:
            raise ValueError(
                f"pathway source {source_id!r} has unknown role {role!r}; "
                f"valid roles: {sorted(VALID_ROLES)}"
            )
        sources.append(
            PathwaySource(
                id=str(source_id),
                name=str(entry.get("name", source_id)),
                role=role,
                redistribution=str(entry.get("redistribution", "")),
                license=str(entry.get("license", "")),
                provenance=dict(entry.get("provenance", {}) or {}),
                description=str(entry.get("description", "")).strip(),
                notes=str(entry.get("notes", "")).strip(),
            )
        )

    return PathwaySourceRegistry(sources=tuple(sources), source_path=registry_path)
