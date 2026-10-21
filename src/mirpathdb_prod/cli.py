"""Command-line interface for the miRPathDB 3.0 production database build.

Database-focused build CLI. Subcommands wire pure library code
(:mod:`mirpathdb_prod.sources`, :mod:`mirpathdb_prod.normalize`,
:mod:`mirpathdb_prod.validation`, :mod:`mirpathdb_prod.reporting`) to the
configured filesystem locations from ``config.yaml``. No prototype code is
imported.

Subcommands:
    inventory-prototype         Inventory the deprecated prototype snapshot.
    extract-legacy-organisms    Extract the prototype/legacy organism baseline.
    build-mirpathdb2-reference  Build the miRPathDB2 organism reference.
    build-new-organisms         Build the new production organism set.
    validate-organism-parity    Validate organism concordance (pass/fail).
    validate-census             Reconciliation gate over the built artefacts.
    validate-diana-methods      DIANA-TarBase method census + release gate.
    compare-go-propagation      GO ``all`` vs ``direct`` propagation sensitivity.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import yaml

from mirpathdb_prod import __version__
from mirpathdb_prod.config import Config, load_config
from mirpathdb_prod.models import Organism
from mirpathdb_prod.normalize.organisms import normalize_organisms
from mirpathdb_prod.normalize.taxonomy import canonical_scientific_name
from mirpathdb_prod.reporting.markdown import render_parity_report
from mirpathdb_prod.reporting.tables import render_organism_tsv
from mirpathdb_prod.sources import legacy, mirpathdb2
from mirpathdb_prod.validation.organism_parity import evaluate_organism_parity

# Source identifier stamped on the generated new-production organism set.
NEW_SOURCE = "new_version"
# Filenames within the configured validation dir.
NEW_ORGANISMS_FILENAME = "organisms_new.tsv"
PARITY_REPORT_MD = "organism_parity_report.md"
PARITY_REPORT_JSON = "organism_parity_report.json"


def _config(args: argparse.Namespace) -> Config:
    return load_config(getattr(args, "config", None))


def _new_source_version(cfg: Config) -> str:
    """Version string stamped on the generated new-production organisms."""
    return f"mirpathdb-prod@{__version__}"


def _parity_policy(cfg: Config) -> dict[str, list[int]]:
    """Read the config-driven floor+expand parity policy (taxid lists).

    Returns the three allowlists that drive the gate — ``hard_floor``,
    ``approved_extensions``, ``approved_removals`` — never silently hard-coded.
    """
    op = cfg.raw.get("validation", {}).get("organism_parity", {}) or {}

    def _ints(key: str) -> list[int]:
        return [int(x) for x in (op.get(key) or [])]

    return {
        "hard_floor": _ints("hard_floor"),
        "approved_extensions": _ints("approved_extensions"),
        "approved_removals": _ints("approved_removals"),
    }


def _build_new_organisms(cfg: Config) -> list[Organism]:
    """Build the new production organism set from config.yaml's organism list."""
    source_version = _new_source_version(cfg)
    records: list[dict[str, Any]] = []
    for o in cfg.organisms:
        tax_id = int(o["canonical_tax_id"])
        # Trust the curated taxonomy for the canonical scientific name; fall back
        # to the config-declared name only if the taxid is uncurated.
        sci = canonical_scientific_name(tax_id) or str(o.get("scientific_name", ""))
        records.append(
            {
                "canonical_tax_id": tax_id,
                "scientific_name": o.get("scientific_name", sci),
                "source_name": o.get("code", ""),
                "common_name": o.get("common_name", ""),
                "source_record_id": o.get("code", ""),
                "normalization_notes": (
                    f"in_mirpathdb2={o.get('in_mirpathdb2', False)}; "
                    f"depth_class={o.get('depth_class', '')}; "
                    f"native_sources={','.join(o.get('native_sources', []) or [])}"
                ),
            }
        )
    organisms = normalize_organisms(records, source=NEW_SOURCE, source_version=source_version)
    return sorted(organisms, key=lambda x: x.canonical_tax_id)


def _load_old(cfg: Config) -> list[Organism]:
    tsv = legacy.legacy_tsv_path(cfg.path("validation_dir"))
    rows = legacy.load_organisms(tsv)
    return normalize_organisms(
        rows, source=legacy.SOURCE, source_version=f"file:{tsv.name}"
    )


def _load_mirpathdb2(cfg: Config) -> list[Organism]:
    tsv = mirpathdb2.mirpathdb2_tsv_path(cfg.path("validation_dir"))
    rows = mirpathdb2.load_organisms(tsv)
    return normalize_organisms(
        rows, source=mirpathdb2.SOURCE, source_version=f"file:{tsv.name}"
    )


# --------------------------------------------------------------------------- #
# Command handlers
# --------------------------------------------------------------------------- #


def _cmd_inventory_prototype(args: argparse.Namespace) -> int:
    print("inventory-prototype is not yet implemented", file=sys.stderr)
    return 2


def _cmd_extract_legacy_organisms(args: argparse.Namespace) -> int:
    cfg = _config(args)
    organisms = _load_old(cfg)
    print(f"loaded {len(organisms)} legacy organism(s) from organisms_old.tsv:")
    for o in sorted(organisms, key=lambda x: x.canonical_tax_id):
        print(f"  {o.canonical_tax_id}\t{o.scientific_name}\t{o.source_name}")
    return 0


def _cmd_build_mirpathdb2_reference(args: argparse.Namespace) -> int:
    cfg = _config(args)
    organisms = _load_mirpathdb2(cfg)
    print(f"loaded {len(organisms)} miRPathDB2 organism(s) from organisms_mirpathdb2.tsv:")
    for o in sorted(organisms, key=lambda x: x.canonical_tax_id):
        print(f"  {o.canonical_tax_id}\t{o.scientific_name}\t{o.source_name}")
    return 0


def _cmd_build_new_organisms(args: argparse.Namespace) -> int:
    cfg = _config(args)
    organisms = _build_new_organisms(cfg)
    out = cfg.path("validation_dir") / NEW_ORGANISMS_FILENAME
    out.parent.mkdir(parents=True, exist_ok=True)
    tsv = render_organism_tsv(organisms)
    out.write_text(tsv, encoding="utf-8")
    print(f"wrote {len(organisms)} new-production organism(s) -> {out}")
    for o in organisms:
        print(f"  {o.canonical_tax_id}\t{o.scientific_name}\t{o.source_name}")
    return 0


def _cmd_validate_organism_parity(args: argparse.Namespace) -> int:
    cfg = _config(args)
    old = _load_old(cfg)
    new = _build_new_organisms(cfg)
    v2 = _load_mirpathdb2(cfg)
    policy = _parity_policy(cfg)

    result = evaluate_organism_parity(
        old,
        new,
        v2,
        hard_floor=policy["hard_floor"],
        approved_extensions=policy["approved_extensions"],
        approved_removals=policy["approved_removals"],
    )

    validation_dir = cfg.path("validation_dir")
    validation_dir.mkdir(parents=True, exist_ok=True)
    md_path = validation_dir / PARITY_REPORT_MD
    json_path = validation_dir / PARITY_REPORT_JSON
    md_path.write_text(render_parity_report(result), encoding="utf-8")
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"organism parity verdict: {result['verdict']}")
    for line in result["summary"]:
        print(f"  {line}")
    print(f"wrote report: {md_path}")
    print(f"wrote report: {json_path}")
    return 0 if result["passed"] else 1


def _cmd_validate_census(args: argparse.Namespace) -> int:
    """Reconciliation gate: matrices non-degenerate + census anchors match.

    Strict by default (``require_present``): an empty enrichment dir or a missing
    anchor manifest FAILS (exit 1). Pass ``--allow-empty`` to tolerate absent
    artefacts (e.g. a partial/dev run).
    """
    from mirpathdb_prod.validation import census

    cfg = _config(args)
    enrichment_dir = cfg.path("data_processed_dir") / "enrichment"
    provenance_dir = cfg.path("provenance_dir")
    contract = cfg.raw.get("census_contract")
    anchors = contract if isinstance(contract, dict) and contract else None
    require_present = not args.allow_empty

    mat = census.validate_matrices_nondegenerate(
        enrichment_dir, require_present=require_present
    )
    cen = census.validate_census(
        provenance_dir, anchors=anchors, require_present=require_present
    )
    print(
        f"census: matrices scanned={mat['n_matrices']} "
        f"degenerate={len(mat['degenerate'])}; "
        f"anchors checked={len(cen['checks'])} failed={len(cen['failures'])}"
    )
    if mat.get("empty") and require_present:
        print(f"  EMPTY: no enrichment matrices under {enrichment_dir}")
    for d in mat["degenerate"]:
        print(f"  DEGENERATE: {d['path']}")
    for f in cen["failures"]:
        print(
            f"  FAIL: {f['anchor']} expected={f.get('expected')} "
            f"got={f.get('got')} {f.get('reason', '')}".rstrip()
        )
    ok = mat["ok"] and cen["ok"]
    print(f"census verdict: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def _cmd_validate_diana_methods(args: argparse.Namespace) -> int:
    """DIANA method census + strict release gate (exit 1 iff an unknown method is HIGH).

    Loads the cached DIANA-TarBase v9 files with the production loaders, writes
    the method census under ``<validation_dir>/diana_methods/`` (plus the
    provenance block unless ``--no-provenance``) and returns 1 when at least one
    ``experimental_method`` outside the CLIP / HIGH / KNOWN_OTHER vocabulary
    contributes HIGH rows. ``--report-only`` is diagnostic: it always returns 0 and
    marks the output non-release-gating. There is deliberately NO bypass flag --
    an unknown HIGH method has to be added to the right vocabulary set on purpose.
    Any input problem -- a sha256 mismatch against ``provenance/manifest.yaml``, a
    missing / unreadable / truncated / corrupt cached file, or the production
    loader's own fail-loud ``ValueError`` (header drift, 0 strong rows, an empty
    species slice) -- is a gate failure (printed as ``diana method gate: FAIL
    (<reason>)``, exit 1), never a traceback. A usage error (an organism code not
    in config.yaml, nothing to audit because config.yaml lists no organisms, an
    ``--out-dir`` that collides with the provenance dir or lies inside
    ``data/processed``) returns 2 without running.

    The provenance block's ``report_path`` is the TSV's repo-relative path when
    ``--out-dir`` is under the repo (the default). For an ``--out-dir`` outside
    the repo it is the absolute path instead -- deterministic per host, but the
    one provenance field that is not byte-identical across hosts (CWD-relative
    would depend on where the command was run from, which is worse).
    """
    import zlib
    from pathlib import Path

    from mirpathdb_prod.config import REPO_ROOT
    from mirpathdb_prod.validation import diana_methods

    def usage_error(msg: str) -> int:
        print(f"validate-diana-methods: {msg}", file=sys.stderr)
        return 2

    cfg = _config(args)
    taxids = {str(o["code"]): int(o["canonical_tax_id"]) for o in cfg.organisms}
    requested = list(args.organisms) if args.organisms else list(taxids)
    unknown = sorted(set(requested) - set(taxids))
    if unknown:
        return usage_error(
            f"unknown organism code(s) {unknown}; config.yaml organisms: {sorted(taxids)}"
        )
    # Canonical config order regardless of how --organisms was spelled.
    organisms = [code for code in taxids if code in set(requested)]
    if not organisms:
        # Nothing was audited, so this is a configuration problem, not a gate FAIL.
        return usage_error("no organisms to audit (config.yaml organisms is empty)")

    # Free-form paths get the same ``~`` expansion config paths get (config._resolve).
    cache_dir = (
        Path(args.cache_dir).expanduser() if args.cache_dir else cfg.path("source_cache_dir")
    )
    out_dir = (
        Path(args.out_dir).expanduser()
        if args.out_dir
        else cfg.path("validation_dir") / "diana_methods"
    )
    provenance_dir = None if args.no_provenance else cfg.path("provenance_dir")
    manifest_path = cfg.path("provenance_dir") / "manifest.yaml"
    if provenance_dir is not None and out_dir.resolve() == provenance_dir.resolve():
        return usage_error(
            f"--out-dir {out_dir} is the provenance dir; the audit JSON and the "
            f"provenance block share the basename {diana_methods.PROVENANCE_JSON} "
            f"(use another --out-dir or --no-provenance)"
        )
    processed = cfg.paths.get("data_processed_dir")
    if processed is not None and out_dir.resolve().is_relative_to(processed.resolve()):
        return usage_error(
            f"--out-dir {out_dir} lies inside data_processed_dir ({processed}); the "
            f"audit never writes into production locations"
        )
    # Repo-relative label of the report for the provenance block; absolute (host-
    # stable, not CWD-relative) when the output lives outside the repo.
    repo_root = cfg.source_path.parent if cfg.source_path else REPO_ROOT
    tsv_path = out_dir / diana_methods.AUDIT_TSV
    try:
        report_path = tsv_path.resolve().relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        report_path = tsv_path.resolve().as_posix()

    try:
        result = diana_methods.run_diana_method_audit(
            cache_dir=cache_dir,
            manifest_path=manifest_path,
            organisms=organisms,
            taxids=taxids,
            out_dir=out_dir,
            provenance_dir=provenance_dir,
            report_only=bool(args.report_only),
            report_path=report_path,
            log=print,
        )
    # RuntimeError: sha256 gate. OSError: missing cached file (FileNotFoundError),
    # a non-gzip / unreadable file (BadGzipFile). EOFError: a truncated gzip
    # ("Compressed file ended before the end-of-stream marker"). zlib.error: a
    # corrupt deflate stream. Neither of the last two subclasses OSError, and an
    # unpinned / absent manifest means the sha256 gate cannot catch them first.
    # ValueError: the production loader's fail-loud guards (header drift, 0
    # strong rows, empty species slice, unmapped taxid).
    except (RuntimeError, ValueError, OSError, EOFError, zlib.error) as exc:
        print(f"diana method gate: FAIL ({exc})")
        return 1
    for line in diana_methods.verdict_lines(result):
        print(line)
    return result.exit_code


def _cmd_compare_go_propagation(args: argparse.Namespace) -> int:
    """GO ``all`` (production) vs ``direct`` propagation sensitivity comparison.

    Reuses the production ``all`` outputs when their build fingerprint matches,
    computes the ``direct`` side in an isolated directory, and writes the
    deterministic comparison artefacts + ``manifest.yaml``. Exit 0 on success, 2
    on a refusal (incompatible ``all`` without ``--recompute-all``; missing inputs
    under ``--reuse-only``), 1 on a hard error (message on stderr, no traceback).
    """
    from mirpathdb_prod.validation import go_propagation_driver as gpd

    cfg = _config(args)
    opts = gpd.ComparisonOptions(
        organisms=tuple(args.organisms),
        tiers=tuple(args.tiers),
        families=tuple(args.families),
        out_dir=Path(args.out_dir) if args.out_dir else None,
        rscript=args.rscript,
        fdr=args.fdr,
        secondary_fdrs=tuple(args.secondary_fdr) if args.secondary_fdr else None,
        structural_only=args.structural_only,
        recompute_all=args.recompute_all,
        reuse_only=args.reuse_only,
        max_cores=args.max_cores,
        argv=["mirpathdb-prod", *getattr(args, "argv", sys.argv[1:])],
    )
    try:
        return gpd.run_go_propagation_comparison(cfg, opts)
    except (ValueError, RuntimeError, OSError, KeyError, TypeError, yaml.YAMLError) as exc:
        # TypeError covers malformed manifest numerics (e.g. a null n_universe);
        # YAMLError an unrepresentable value while writing manifest.yaml.
        print(f"[go-sens] ERROR: {exc}", file=sys.stderr)
        return 1


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level argument parser with all subcommands."""
    parser = argparse.ArgumentParser(
        prog="mirpathdb-prod",
        description="miRPathDB 3.0 production database build (database-focused).",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "-c",
        "--config",
        default=None,
        help="Path to config.yaml (defaults to the repo-root config.yaml).",
    )

    subparsers = parser.add_subparsers(dest="command", metavar="<command>")
    subparsers.required = True

    p = subparsers.add_parser(
        "inventory-prototype",
        help="Inventory the deprecated prototype snapshot (reference only) (not yet implemented).",
    )
    p.set_defaults(func=_cmd_inventory_prototype)

    p = subparsers.add_parser(
        "extract-legacy-organisms",
        help="Extract the prototype/legacy organism baseline.",
    )
    p.set_defaults(func=_cmd_extract_legacy_organisms)

    p = subparsers.add_parser(
        "build-mirpathdb2-reference",
        help="Build the miRPathDB2 organism reference (primary validation source).",
    )
    p.set_defaults(func=_cmd_build_mirpathdb2_reference)

    p = subparsers.add_parser(
        "build-new-organisms",
        help="Build the new production organism set (hsa/mmu/rno/mml).",
    )
    p.set_defaults(func=_cmd_build_new_organisms)

    p = subparsers.add_parser(
        "validate-organism-parity",
        help="Validate organism concordance and emit an explicit pass/fail result.",
    )
    p.set_defaults(func=_cmd_validate_organism_parity)

    p = subparsers.add_parser(
        "validate-census",
        help="Reconciliation gate: enrichment matrices non-degenerate + census "
        "anchors match (strict by default; exit 1 on FAIL).",
    )
    p.add_argument(
        "--allow-empty",
        action="store_true",
        help="Tolerate an empty enrichment dir / missing anchor manifests "
        "instead of failing (default: absent artefacts FAIL the gate).",
    )
    p.set_defaults(func=_cmd_validate_census)

    p = subparsers.add_parser(
        "validate-diana-methods",
        help="DIANA-TarBase method census + release gate: FAILS (exit 1) iff any "
        "experimental_method outside the CLIP/HIGH/KNOWN_OTHER vocabulary "
        "contributes HIGH rows; WARN (exit 0) if unknown methods exist but none "
        "is HIGH; PASS otherwise. No bypass flag. Exit 2 = usage error.",
    )
    p.add_argument(
        "--report-only",
        action="store_true",
        help="Diagnostic run: write the census, print the verdict as WOULD FAIL / "
        "WARN / PASS and always exit 0 (output is marked NOT release-gating).",
    )
    p.add_argument(
        "--organisms",
        nargs="+",
        default=None,
        metavar="CODE",
        help="Organism codes to audit (default: every config.yaml organism, "
        "hsa mmu rno mml).",
    )
    p.add_argument(
        "--cache-dir",
        default=None,
        help="Cache holding the DIANA-TarBase v9 downloads "
        "(default: config.yaml paths.source_cache_dir).",
    )
    p.add_argument(
        "--out-dir",
        default=None,
        help="Where the census artefacts are written "
        "(default: <paths.validation_dir>/diana_methods).",
    )
    p.add_argument(
        "--no-provenance",
        action="store_true",
        help="Skip writing <paths.provenance_dir>/diana_method_audit.json.",
    )
    p.set_defaults(func=_cmd_validate_diana_methods)

    p = subparsers.add_parser(
        "compare-go-propagation",
        help="GO propagation sensitivity: compare the production `all` (GOALL) GO "
        "enrichment with an isolated `direct` (GO) run; never modifies "
        "data/processed/ or config.yaml (exit 2 = refusal, 1 = hard error).",
        description=(
            "Compare GO gene-set structure and ORA significance under go_propagation "
            "'all' (production default, GOALL/ONTOLOGYALL) vs 'direct' (GO/ONTOLOGY). "
            "The 'all' side reuses the production matrices when their build fingerprint "
            "(versions, target-set census, universes, cells) matches; the 'direct' side "
            "is computed with the production loaders + clusterProfiler kernel into an "
            "isolated directory and reused on re-runs. Artefacts: "
            "go_propagation_summary.json, go_structural_comparison.tsv, "
            "go_structural_terms.tsv.gz, go_cell_comparison.tsv, go_terms_gained.tsv, "
            "go_terms_lost.tsv (+ *_full.tsv.gz), GO_PROPAGATION_SENSITIVITY.md, "
            "manifest.yaml (the only file with timestamps). Side status vocabulary in "
            "manifest.yaml: reused_production | reused_isolated | computed | skipped | "
            "refused | mixed (per-organism statuses under status_by_organism)."
        ),
    )
    p.add_argument(
        "--organisms",
        nargs="+",
        default=["hsa", "mmu", "rno", "mml"],
        choices=["hsa", "mmu", "rno", "mml"],
        help="Organisms to compare (default: all four).",
    )
    p.add_argument(
        "--tiers",
        nargs="+",
        default=["HIGH", "MEDIUM", "LOW"],
        choices=["HIGH", "MEDIUM", "LOW"],
        help="Cumulative target tiers (default: HIGH MEDIUM LOW).",
    )
    p.add_argument(
        "--families",
        nargs="+",
        default=["GO_BP", "GO_CC", "GO_MF"],
        choices=["GO_BP", "GO_CC", "GO_MF"],
        help="GO families (default: GO_BP GO_CC GO_MF).",
    )
    p.add_argument(
        "--out-dir",
        default=None,
        help="Comparison root (default: <validation_dir>/go_propagation_sensitivity). "
        "Must not live under data_processed_dir.",
    )
    p.add_argument(
        "--rscript",
        default=None,
        help="Rscript binary (default: $MPD_RSCRIPT, else PATH). Exported as MPD_RSCRIPT "
        "for the run so the production target-set loaders use the same binary.",
    )
    p.add_argument(
        "--fdr",
        type=float,
        default=None,
        help="Primary significance threshold (p_adjust <=). Default: config key "
        "validation.go_propagation.fdr_significant if set, else the code constant "
        "reporting.significant.FDR_SIGNIFICANT (0.05, the production criterion).",
    )
    p.add_argument(
        "--secondary-fdr",
        nargs="+",
        type=float,
        default=None,
        help="Secondary sensitivity thresholds reported alongside (default: 0.005); "
        "NOT the production criterion.",
    )
    p.add_argument(
        "--structural-only",
        action="store_true",
        help="Skip the enrichment comparison; still write the structural artefacts + "
        "manifest.yaml.",
    )
    p.add_argument(
        "--recompute-all",
        action="store_true",
        help="Allow an isolated `all` side when the production `all` output is "
        "incompatible: a fingerprint-compatible run already under <out-dir>/all/ is "
        "reused, else one is computed there. Without this flag an incompatible `all` "
        "side is a refusal listing every mismatch (exit 2), even when such an isolated "
        "run exists.",
    )
    p.add_argument(
        "--reuse-only",
        action="store_true",
        help="Never compute: if a side is missing or incompatible, report exactly what "
        "and exit 2.",
    )
    p.add_argument(
        "--max-cores",
        type=int,
        default=None,
        help="Passed through as MPD_MAX_CORES to the clusterProfiler batch driver.",
    )
    p.set_defaults(func=_cmd_compare_go_propagation)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entrypoint. Returns a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)
    # The exact invocation is recorded by commands that write a run manifest.
    args.argv = list(argv) if argv is not None else sys.argv[1:]
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
