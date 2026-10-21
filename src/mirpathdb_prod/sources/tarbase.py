"""DIANA-TarBase target-source adapter (fetch/load, real implementation).

DIANA-TarBase v9 is a **PARALLEL VALIDATED resource to miRTarBase** — it carries
both verified (experimentally validated) interactions AND CLIP/AGO support; it is
NOT a CLIP-only source. Accordingly it feeds **BOTH HIGH and MEDIUM**: its
verified-strong evidence maps to HIGH (``val_strong``) and its CLIP/weaker
experimental evidence maps to MEDIUM (``val_weak`` / ``clip_supported``). The
exact verified-strong vs CLIP/weaker tier boundary is an open item flagged for
Tuesday — see ``docs/target_evidence.md``.

The resource spans **24 species INCLUDING macaque (mml, 9544)** — its hsa/mmu-only
ingest in the prototype was a **code limit, not a data limit** (the prototype
fetcher was hard-coded to hsa+mmu and rejected every other organism). For macaque,
DIANA's per-organism VOLUME is **exactly ONE validated MTI**
(``mml-miR-125b-5p -> ANGPT2``) — that lone interaction is **marked** into macaque's
HIGH band (MEDIUM empty); macaque is otherwise **predicted-dominant**, its
substantive evidence being the TargetScan-predicted LOW band. miRTarBase and miRDB
do NOT cover macaque.

File layout (verified against the live v9 downloads page,
``https://dianalab.e-ce.uth.gr/tarbasev9/downloads``): the resource ships FOUR
download files, NOT one-per-species — ``Homo_sapiens``, ``Mus_musculus``, a
``Viral_species`` bundle, and an ``Other_species`` bundle that holds **every other
species** (rat and macaque among them). So hsa/mmu load from their own per-species
files, while **rat (rno, 10116) and macaque (mml, 9544) load from
``Other_species_TarBase-v9.tsv.gz`` filtered by the ``species`` column** (analogous
to how TargetScan slices its multi-species file by Gene Tax ID). The production
extension therefore (a) maps each candidate taxid to its ``species`` string and (b)
adds a multi-species slice loader (:func:`iter_mti_multispecies` /
:func:`load_mti_multispecies`). In v9 the ``Other_species`` bundle carries ~1125
rat MTIs (26 miRNAs) but only a SINGLE validated macaque MTI
(``mml-miR-125b-5p -> ANGPT2``, a Luciferase Reporter Assay) — that is the real
extent of macaque's validated coverage; its bulk (predicted/LOW) comes from
TargetScan. (Dog/chicken were evaluated + excluded by owner choice.) Columns
(verified against the cached hsa/mmu files AND the ``Other_species`` bundle)::

    species | mirna_name | mirna_id | gene_name | gene_id | gene_location |
    transcript_name | transcript_id | chromosome | start | end | strand |
    experimental_method | regulation | tissue | cell_line | article_pubmed_id |
    confidence | interaction_group | cell_type | microt_score | comment

``mirna_name`` is the miRBase display name; ``gene_name`` is a gene SYMBOL and
``gene_id`` an Ensembl gene id. Fetch/load ONLY; no normalization/judgment; no
prototype code imported.
"""

from __future__ import annotations

import csv
import gzip
import warnings
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

from mirpathdb_prod.sources import validate_species_column

SOURCE = "diana_tarbase"
SOURCE_VERSION = "v9"

# Validated + CLIP (parallel to miRTarBase) -> val_strong (HIGH) + val_weak /
# clip_supported (MEDIUM). DIANA feeds BOTH HIGH and MEDIUM, not MEDIUM only.
PROVENANCE_TIERS: tuple[str, ...] = ("val_strong", "val_weak", "clip_supported")

# hsa/mmu are confirmed v9 (the prototype only ever INGESTED these — a code limit).
# DIANA-TarBase v9 covers 24 species; rat (10116) + macaque (9544) load from the
# multi-species Other_species bundle, sliced by the ``species`` column. Verified
# VOLUME: rat ~1125 MTIs / 26 miRNAs; macaque EXACTLY ONE validated MTI
# (mml-miR-125b-5p -> ANGPT2) — that lone interaction is marked into macaque's HIGH
# band (macaque is otherwise predicted-dominant, its substantive evidence being the
# TargetScan LOW band).
CONFIRMED_TAXIDS: frozenset[int] = frozenset({9606, 10090, 10116, 9544})
# Retained for back-compat with callers/tests that distinguished the formerly
# "to-be-extended" organisms; now EMPTY (rat + macaque are confirmed + wired).
CANDIDATE_TAXIDS: frozenset[int] = frozenset()

# Taxids served by the multi-species ``Other_species`` bundle (NOT their own
# per-species file). These are sliced by the ``species`` column.
OTHER_SPECIES_TAXIDS: frozenset[int] = frozenset({10116, 9544})

# Filename of the multi-species bundle that holds every non-hsa/mmu/viral species.
OTHER_SPECIES_FILENAME = "Other_species_TarBase-v9.tsv.gz"

# Canonical taxid -> the ``species`` column value used to slice the multi-species
# bundle (verified against the cached ``Other_species`` file).
_TAXID_TO_SPECIES: dict[int, str] = {
    10116: "Rattus norvegicus",
    9544: "Macaca mulatta",
}

# DIANA experimental methods that are CLIP/AGO-based (transcriptome-wide binding
# evidence) rather than gene-specific functional assays. These map to MEDIUM
# (``clip_supported``); gene-specific assays (Luciferase/Western/qPCR/...) map to
# HIGH (``val_strong``). Verified against the v9 method vocabulary.
CLIP_METHODS: frozenset[str] = frozenset(
    {
        "HITS-CLIP",
        "PAR-CLIP",
        "CLASH",
        "qCLASH",
        "iCLIP",
        "CLIP-Seq",
        "AGO-IP",
        "RIP-Seq",
        "Chimeric fragments",
        "IMPACT-Seq",
        "RPF-Seq",
        "Biotin-Microarrays",
        "Biotin-qPCR",
    }
)

# Gene-specific functional assays (NOT CLIP) that, with ``confidence == "1"``,
# roll up to HIGH (see :func:`rollup_tier`). Listed so :func:`validate_vocabulary`
# can catch a method-vocabulary drift (M1). Not exhaustive of DIANA's full method
# list -- an unseen method simply warns; it does not change the rollup.
HIGH_METHODS: frozenset[str] = frozenset(
    {
        "Luciferase Reporter Assay",
        "Western Blot",
        "qPCR",
        "Northern Blot",
        "GFP Reporter Assay",
        "Reporter Assay",
        "5RACE",
        "Immunohistochemistry",
        "ELISA",
        "Flow Cytometry",
        # TarBase v9 label resolved 2026-08-30 (owner-approved; see CHANGELOG):
        # 3LIFE -- HIGH: direct functional 3'UTR reporter assay (luminescence-
        # based 3'UTR screen of a defined target sequence). Was HIGH via the
        # fail-open rule; now curated.
        "3LIFE",
    }
)

# Known non-pair-specific methods observed in the real DIANA-TarBase v9 files
# (hsa/mmu/Other_species): transcriptome-wide readouts (RNA-Seq / Degradome /
# Microarrays / ...), proteome-wide (proteomic-response) readouts (pSILAC / SILAC /
# 2D-DIGE) and affinity-capture / binding-only evidence (Biotin-Seq / TRAP / EMSA).
# These roll up as MEDIUM/weak (NOT HIGH) but are part of the authoritative
# vocabulary, so they must NOT warn. Seeded from the distinct ``experimental_method``
# values in the cached v9 files; the seven TarBase v9 labels the 2026-08-29 release
# gate surfaced were resolved on 2026-08-30 (owner-approved, one comment each).
KNOWN_OTHER_METHODS: frozenset[str] = frozenset(
    {
        "RNA-Seq",
        "sRNA-Seq",
        "Degradome",
        "RT-qPCR",
        "Other",
        "Microarrays",
        "LAMP",
        "In Situ Hybridization",
        "Immunofluorescence",
        "Immunoblot",
        # --- TarBase v9 labels resolved 2026-08-30 (owner-approved, MEDIUM) ---
        # 2D-DIGE -- proteomic response (2D difference gel electrophoresis);
        # not pair-specific direct targeting.
        "2D-DIGE",
        # Biotin-Seq -- affinity-capture / association evidence (biotinylated-
        # miRNA pulldown + sequencing); kept in KNOWN_OTHER, no new class.
        "Biotin-Seq",
        # Electrophoretic Mobility Shift Assay -- single mmu row (mmu-miR-93-5p ->
        # Sp7; Yang et al. 2012, J Bone Miner Res, PMID 22467200): the EMSA there
        # supports the Sp7 -> miR-93 promoter (transcriptional) arm of the
        # feedback loop, i.e. protein-DNA binding, not direct miRNA-mRNA binding
        # for the pair; the pair keeps HIGH through its Western Blot + qPCR rows.
        # Fail-closed.
        "Electrophoretic Mobility Shift Assay",
        # Genetic Testing -- DIANA label for the 3'UTR-library positive/negative-
        # selection functional screen of Gäken et al. 2012, Nucleic Acids Res,
        # PMID 22323518 (miR-130a / miR-10a hits in MCF-7): library-scale genetic
        # selection, not pair-specific target-site mutation/rescue. Fail-closed.
        "Genetic Testing",
        # SILAC -- proteome-wide SILAC response; not pair-specific.
        "SILAC",
        # TRAP -- DIANA label for miTRAP (miRNA trapping by RNA in vitro affinity
        # purification; Braun et al. 2014, Nucleic Acids Res, PMID 24510096):
        # in-vitro affinity purification of miRNAs binding an in-vitro-transcribed
        # 3'UTR (MYC, ZEB2) in cell lysates -- binding / affinity capture, not
        # functional repression.
        "TRAP",
        # pSILAC -- proteome-wide pulsed-SILAC response; not pair-specific direct
        # targeting (was HIGH via the fail-open rule until 2026-08-30).
        "pSILAC",
    }
)

# The full EXPECTED experimental_method vocabulary = CLIP/AGO (MEDIUM) ∪ the
# gene-specific functional assays (HIGH) ∪ the known transcriptome-wide methods.
# Only a value outside ALL of these (genuine drift) warns (additive).
EXPECTED_METHODS: frozenset[str] = CLIP_METHODS | HIGH_METHODS | KNOWN_OTHER_METHODS

# The EXPECTED ``confidence`` column vocabulary. DIANA v9 uses "1" for verified
# interactions; non-"1"/empty maps to MEDIUM. A value outside this set warns.
EXPECTED_CONFIDENCE: frozenset[str] = frozenset({"1", "0", ""})

# Human-readable label for the HIGH/MEDIUM split, recorded in the target-set
# census so a manifest self-documents which tiering rule produced it (the rule
# changed with the HIGH-2 fix: transcriptome-wide methods are now MEDIUM).
#
# The method-census release gate (``mirpathdb-prod validate-diana-methods``,
# :mod:`mirpathdb_prod.validation.diana_methods`) records this very string as its
# ``classifier_rule`` and classifies rows with the SAME functions production
# uses -- :func:`normalize_method` / :func:`normalize_confidence` /
# :func:`rollup_tier` / :func:`classify_method` -- so the audit can never
# disagree with the band builders about which rows are HIGH.
#
# The rule text below is deliberately NOT versioned by vocabulary edits: the
# 2026-08-30 vocabulary revision (3LIFE -> HIGH_METHODS; pSILAC / SILAC / 2D-DIGE /
# Biotin-Seq / TRAP / Genetic Testing / Electrophoretic Mobility Shift Assay ->
# KNOWN_OTHER_METHODS) changed the SETS the rule consults, not the rule itself.
# It is recorded in CHANGELOG.md, provenance/BEHAVIOR_CHANGES.md and
# provenance/manifest.yaml (release_safeguards.diana_method_audit).
DIANA_TIER_RULE = (
    "HIGH = confidence==1 AND gene-specific assay "
    "(method not in CLIP_METHODS and not in KNOWN_OTHER_METHODS); else MEDIUM"
)

# Vocabulary classes of a NORMALIZED ``experimental_method`` value, as reported by
# :func:`classify_method`. These describe which of the three vocabulary sets a
# method belongs to; they are a CENSUS label, NOT the tier rule (an ``unknown``
# method at confidence 1 still rolls up HIGH -- fail-open -- which is exactly
# what the audit gate exists to surface).
VOCAB_KNOWN_HIGH = "known_high"
VOCAB_KNOWN_CLIP = "known_clip"
VOCAB_KNOWN_OTHER = "known_other"
VOCAB_UNKNOWN = "unknown"
VOCABULARY_CLASSES: tuple[str, ...] = (
    VOCAB_KNOWN_HIGH,
    VOCAB_KNOWN_CLIP,
    VOCAB_KNOWN_OTHER,
    VOCAB_UNKNOWN,
)


def normalize_method(value: Any) -> str:
    """Production normalization of a raw ``experimental_method`` cell.

    Exactly ``(value or "").strip()``: surrounding whitespace is removed, the
    case is NOT folded and no alias table is applied -- a value is compared
    against :data:`CLIP_METHODS` / :data:`HIGH_METHODS` /
    :data:`KNOWN_OTHER_METHODS` verbatim. ``None`` / empty become ``""``. This is
    the single normalization :func:`rollup_tier`, :func:`validate_vocabulary` and
    the method-census audit share, so a spelling variant (``"qPCR "``, ``"qpcr"``)
    is treated identically everywhere.
    """
    return (value or "").strip()


def normalize_confidence(value: Any) -> str:
    """Production normalization of a raw ``confidence`` cell (``(value or "").strip()``).

    DIANA v9 encodes a verified interaction as the string ``"1"``; the strong
    criterion compares the stripped value against exactly that literal.
    """
    return (value or "").strip()


def classify_method(method: str) -> str:
    """Vocabulary class of an ALREADY-normalized method (see :func:`normalize_method`).

    Precedence: :data:`CLIP_METHODS` -> :data:`VOCAB_KNOWN_CLIP`;
    :data:`HIGH_METHODS` -> :data:`VOCAB_KNOWN_HIGH`; :data:`KNOWN_OTHER_METHODS`
    -> :data:`VOCAB_KNOWN_OTHER`; anything else (including the empty string) ->
    :data:`VOCAB_UNKNOWN`. The three sets are pairwise disjoint today (asserted
    by the test-suite), so the precedence is documentation, not arbitration; it
    mirrors :func:`_is_high_evidence`, where the CLIP / KNOWN_OTHER exclusions
    outrank everything else.

    This is a census label only: it does NOT feed :func:`rollup_tier`, whose
    fail-open rule bands an ``unknown`` method HIGH at confidence 1.
    """
    if method in CLIP_METHODS:
        return VOCAB_KNOWN_CLIP
    if method in HIGH_METHODS:
        return VOCAB_KNOWN_HIGH
    if method in KNOWN_OTHER_METHODS:
        return VOCAB_KNOWN_OTHER
    return VOCAB_UNKNOWN


def is_unknown_method(method: str) -> bool:
    """True iff a normalized method is outside all three vocabulary sets."""
    return classify_method(method) == VOCAB_UNKNOWN


def _is_high_evidence(method: str, confidence: str) -> bool:
    """True iff a DIANA row is verified-strong = HIGH-tier evidence.

    HIGH = a verified (``confidence == "1"``) interaction from a **gene-specific
    functional assay**. Both the CLIP/AGO transcriptome-wide methods
    (:data:`CLIP_METHODS`) AND the documented non-pair-specific
    :data:`KNOWN_OTHER_METHODS` (transcriptome-wide RNA-Seq / sRNA-Seq /
    Degradome / RT-qPCR / Microarrays / ..., proteome-wide pSILAC / SILAC /
    2D-DIGE, affinity-capture Biotin-Seq / TRAP / EMSA) roll up to MEDIUM even
    at confidence 1.

    This closes the review's HIGH-2 defect: the previous rule was a CLIP-only
    *denylist* (``method not in CLIP_METHODS``), so every ``KNOWN_OTHER_METHODS``
    value -- which the module explicitly documents as MEDIUM/weak -- was banded
    HIGH at confidence 1, inflating the strongly-validated tier and contradicting
    both that comment and :func:`rollup_tier`'s docstring. Excluding
    ``KNOWN_OTHER_METHODS`` fixes exactly that. A *genuinely unknown* method
    (outside CLIP, HIGH and KNOWN_OTHER) is still treated as strong here
    (fail-open, unchanged) so a real-but-unlisted gene-specific assay is not
    silently demoted; :func:`validate_vocabulary` warns on it so the drift is
    surfaced and can be added to the right list deliberately. (A stricter pure
    ``method in HIGH_METHODS`` allowlist was considered but rejected: HIGH_METHODS
    is explicitly non-exhaustive, so it would risk demoting real evidence and
    would flip the fail-loud drift behaviour.)
    """
    return (
        confidence == "1"
        and bool(method)
        and method not in CLIP_METHODS
        and method not in KNOWN_OTHER_METHODS
    )


def validate_vocabulary(
    rows: Iterable[dict[str, Any]], *, require_strong: bool = False
) -> None:
    """Warn on unknown method/confidence values; raise if NO strong row exists (M1).

    Additive guard that does NOT change :func:`rollup_tier`: it scans the loaded
    rows' ``experimental_method`` + ``confidence`` columns and ``warnings.warn``s
    once per unknown value (so a renamed/added DIANA category surfaces rather than
    being silently bucketed MEDIUM).

    When ``require_strong`` is True (the production per-species + Other_species
    slices all carry verified-strong evidence) it ALSO raises ``ValueError`` if
    ZERO rows meet the strong criterion (``confidence == "1"`` AND a non-CLIP
    gene-specific method) among rows that carry a confidence value -- the signal
    that the verified-strong vocabulary drifted and would silently empty the HIGH
    band. Left False, the strong-absence is warn-only (a CLIP-only slice is a
    legitimate no-strong case). Macaque's lone strong MTI satisfies the check, so
    ``require_strong=True`` does not fire on that tiny-but-valid slice.
    """
    methods: set[str] = set()
    confidences: set[str] = set()
    n_conf_rows = 0
    n_strong = 0
    for row in rows:
        method = normalize_method(row.get("experimental_method"))
        confidence = normalize_confidence(row.get("confidence"))
        if method:
            methods.add(method)
        confidences.add(confidence)
        if confidence:
            n_conf_rows += 1
        if _is_high_evidence(method, confidence):
            n_strong += 1
    for value in sorted(methods - EXPECTED_METHODS):
        warnings.warn(
            f"DIANA-TarBase: unknown experimental_method {value!r} not in "
            f"EXPECTED_METHODS -- rolled up as MEDIUM (non-HIGH) unless it is a "
            f"gene-specific assay; confirm the method vocabulary did not drift.",
            stacklevel=2,
        )
    for value in sorted(confidences - EXPECTED_CONFIDENCE):
        warnings.warn(
            f"DIANA-TarBase: unknown confidence {value!r} not in "
            f"EXPECTED_CONFIDENCE {sorted(EXPECTED_CONFIDENCE)} -- treated as "
            f"non-verified (MEDIUM); confirm the confidence encoding did not drift.",
            stacklevel=2,
        )
    if require_strong and n_conf_rows and n_strong == 0:
        raise ValueError(
            "DIANA-TarBase: 0 rows meet the strong criterion (confidence == '1' "
            f"AND a non-CLIP gene-specific method) among {n_conf_rows} rows with a "
            "confidence value. The verified-strong vocabulary appears to have "
            "drifted -- refusing to silently empty the HIGH band."
        )


def species_for_taxid(taxid: int) -> str:
    """Return the ``species`` column value for a multi-species-bundle organism."""
    try:
        return _TAXID_TO_SPECIES[taxid]
    except KeyError as exc:
        raise ValueError(
            f"taxid {taxid} is not served by the DIANA Other_species bundle; "
            f"known: {sorted(_TAXID_TO_SPECIES)}"
        ) from exc


def rollup_tier(row: dict[str, Any]) -> str:
    """Map ONE DIANA MTI row to the user-facing HIGH/MEDIUM confidence band.

    Per the locked rollup, DIANA-TarBase is a PARALLEL VALIDATED resource that
    feeds BOTH HIGH and MEDIUM. The split here is:

    * ``HIGH`` -- a verified (``confidence == "1"``) interaction backed by a
      gene-specific functional assay (Luciferase reporter / Western blot / qPCR /
      3LIFE / ...): a method in NEITHER :data:`CLIP_METHODS` NOR the
      non-pair-specific :data:`KNOWN_OTHER_METHODS`.
    * ``MEDIUM`` -- CLIP/AGO transcriptome-wide support (:data:`CLIP_METHODS`),
      the documented non-pair-specific :data:`KNOWN_OTHER_METHODS` (RNA-Seq /
      Degradome / Microarrays / RT-qPCR / pSILAC / SILAC / 2D-DIGE / Biotin-Seq /
      TRAP / ...), OR a non-verified (``confidence != "1"``) row.

    See :func:`_is_high_evidence` for the exact rule (shared with the
    ``require_strong`` guard so the two cannot drift).

    Args:
        row: A raw DIANA MTI row dict keyed by :data:`TARBASE_COLUMNS`.

    Returns:
        ``"HIGH"`` or ``"MEDIUM"``.
    """
    method = normalize_method(row.get("experimental_method"))
    confidence = normalize_confidence(row.get("confidence"))
    return "HIGH" if _is_high_evidence(method, confidence) else "MEDIUM"

# Column header of the per-organism DIANA-TarBase v9 TSV (verified against the
# cached hsa/mmu files). A reshaped file fails loudly.
TARBASE_COLUMNS: tuple[str, ...] = (
    "species",
    "mirna_name",
    "mirna_id",
    "gene_name",
    "gene_id",
    "gene_location",
    "transcript_name",
    "transcript_id",
    "chromosome",
    "start",
    "end",
    "strand",
    "experimental_method",
    "regulation",
    "tissue",
    "cell_line",
    "article_pubmed_id",
    "confidence",
    "interaction_group",
    "cell_type",
    "microt_score",
    "comment",
)

# DIANA-TarBase v9 download filename per taxid (verified against the live
# downloads page). hsa/mmu have their OWN per-species file; rat + macaque are NOT
# shipped as per-species files — they live in the multi-species Other_species
# bundle, so they map to OTHER_SPECIES_FILENAME and must be loaded with the
# species-slicing loaders (:func:`load_mti_multispecies`).
_TAXID_TO_FILENAME: dict[int, str] = {
    9606: "Homo_sapiens_TarBase-v9.tsv.gz",
    10090: "Mus_musculus_TarBase-v9.tsv.gz",
    10116: OTHER_SPECIES_FILENAME,  # rat: in the multi-species bundle
    9544: OTHER_SPECIES_FILENAME,   # macaque: in the multi-species bundle
}


def filename_for_taxid(taxid: int) -> str:
    """Return the DIANA-TarBase v9 download filename for an organism.

    For rat (10116) and macaque (9544) this returns the multi-species
    ``Other_species`` bundle name (they have no per-species file); callers must
    then slice by the ``species`` column via :func:`load_mti_multispecies`.
    """
    try:
        return _TAXID_TO_FILENAME[taxid]
    except KeyError as exc:
        raise ValueError(
            f"no DIANA-TarBase v9 filename mapping for taxid {taxid}; "
            f"known: {sorted(_TAXID_TO_FILENAME)}"
        ) from exc


def load_mti(path: Path, *, taxid: int) -> list[dict[str, Any]]:
    """Load raw experimental/CLIP MTIs for one organism (streaming gzip TSV).

    Production extends coverage beyond the prototype's hsa+mmu (a code limit) to
    the candidate organisms in :data:`CANDIDATE_TAXIDS` — macaque (in the v9
    24-species set) and rat — once their v9 TSVs are confirmed; the prototype's
    hard rejection of non-hsa/mmu organisms is deliberately NOT reproduced (we only
    require the file to exist).

    Pure load: returns one dict per data row keyed by :data:`TARBASE_COLUMNS`
    (values verbatim — no normalization, no judgment).

    Args:
        path: Path to the DIANA-TarBase v9 per-organism gzip TSV.
        taxid: Canonical NCBI taxid to load (recorded for provenance). NOTE this
            pins the TARGET-GENE species only: a per-species DIANA file is
            partitioned by the species of the target gene, and its ``mirna_name``
            column legitimately carries miRNAs of OTHER species and of viruses
            (e.g. mmu-/rno-/ebv- names inside the human file) that were assayed
            against those genes. miRNA species and target species are two
            independent axes; do not read this argument as "single species".
            The miRNA axis is filtered downstream, at contract export
            (scripts/export_frontend_contracts.py, retain_row).

    Returns:
        Raw, un-normalized MTI row dicts (verbatim; no judgment).

    Raises:
        ValueError: if the header does not match :data:`TARBASE_COLUMNS`, or if
            the majority of rows' ``species`` values contradict ``taxid``
            (wrong-file detection; a minority mismatch warns and keeps rows).
        FileNotFoundError: if ``path`` does not exist.
    """
    rows = list(iter_mti(path, taxid=taxid))
    # Row-level species check: the file's own target-species annotations must
    # agree with the requested taxid (protection beyond the filename
    # convention; the miRNA axis stays unfiltered, see the docstring).
    validate_species_column(
        rows, column="species", taxid=taxid, path=path, source_label="DIANA-TarBase"
    )
    # M1: additive vocab check. hsa/mmu per-species DIANA carry verified-strong
    # evidence, so require it (a drift that empties HIGH raises).
    validate_vocabulary(rows, require_strong=True)
    return rows


def iter_mti(path: Path, *, taxid: int) -> list[dict[str, Any]]:
    """Eager loader shared by :func:`load_mti` (kept as a hook for streaming).

    Returns the same list :func:`load_mti` does; separated so a future caller can
    swap in a generator without changing the public signature. ``taxid`` is
    accepted for symmetry / future per-species validation.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"DIANA-TarBase v9 TSV not found: {path}")

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = tuple(reader.fieldnames or ())
        if header != TARBASE_COLUMNS:
            raise ValueError(
                f"{path}: unexpected DIANA-TarBase header\n"
                f"  expected: {TARBASE_COLUMNS}\n"
                f"  got:      {header}"
            )
        return [dict(row) for row in reader]


def iter_mti_multispecies(path: Path, *, taxid: int) -> Iterator[dict[str, Any]]:
    """Stream MTI rows for ONE organism out of the multi-species bundle.

    The ``Other_species`` bundle holds many species in a single TSV (same
    :data:`TARBASE_COLUMNS` header). This yields only the rows whose ``species``
    column matches the organism's canonical species string
    (:func:`species_for_taxid`) — the DIANA analogue of TargetScan's
    Gene-Tax-ID slice.

    Args:
        path: Path to the multi-species gzip TSV (``Other_species_TarBase-v9.tsv.gz``).
        taxid: Canonical NCBI taxid to slice out (must be in
            :data:`OTHER_SPECIES_TAXIDS`).

    Yields:
        Raw MTI row dicts (verbatim; no normalization, no judgment).

    Raises:
        ValueError: if ``taxid`` is not a bundle organism or the header is wrong.
        FileNotFoundError: if ``path`` does not exist.
    """
    if taxid not in OTHER_SPECIES_TAXIDS:
        raise ValueError(
            f"taxid {taxid} is not served by the DIANA Other_species bundle; "
            f"supported: {sorted(OTHER_SPECIES_TAXIDS)}"
        )
    wanted = species_for_taxid(taxid)
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"DIANA-TarBase v9 multi-species bundle not found: {path}")

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = tuple(reader.fieldnames or ())
        if header != TARBASE_COLUMNS:
            raise ValueError(
                f"{path}: unexpected DIANA-TarBase header\n"
                f"  expected: {TARBASE_COLUMNS}\n"
                f"  got:      {header}"
            )
        for row in reader:
            if (row.get("species") or "").strip() == wanted:
                yield dict(row)


def load_mti_multispecies(path: Path, *, taxid: int) -> list[dict[str, Any]]:
    """Eager wrapper of :func:`iter_mti_multispecies` (bundle -> one organism).

    Fail-loud (M2): raises ``ValueError`` if the species slice yields 0 rows for a
    requested production taxid -- a 0-row slice almost always means the species
    string drifted or the bundle was reshaped, not a real absence (rat carries
    ~1125 MTIs; macaque carries exactly one). The streaming
    :func:`iter_mti_multispecies` stays a generator (no eager guard) for callers
    that legitimately stream/probe.
    """
    rows = list(iter_mti_multispecies(Path(path), taxid=taxid))
    if not rows:
        raise ValueError(
            f"DIANA-TarBase Other_species slice for taxid {taxid} "
            f"({species_for_taxid(taxid)!r}) yielded 0 rows from {path}: expected "
            f"the species' MTIs (rat ~1125; macaque 1). Likely a drifted species "
            f"string or a reshaped bundle, not a real absence."
        )
    # M1: additive vocab check. The production Other_species slices (rat, macaque)
    # both carry >=1 verified-strong MTI (macaque's lone validated MTI is strong),
    # so require it -- a drift that empties HIGH raises.
    validate_vocabulary(rows, require_strong=True)
    return rows


def load_mti_for_taxid(path: Path, *, taxid: int) -> list[dict[str, Any]]:
    """Dispatch loader: per-species file (hsa/mmu) vs the multi-species bundle.

    For hsa/mmu this is :func:`load_mti`; for rat/macaque (in
    :data:`OTHER_SPECIES_TAXIDS`) it slices the bundle via
    :func:`load_mti_multispecies`. Lets callers say "load DIANA for taxid X"
    without knowing the file layout.
    """
    if taxid in OTHER_SPECIES_TAXIDS:
        return load_mti_multispecies(path, taxid=taxid)
    return load_mti(path, taxid=taxid)
