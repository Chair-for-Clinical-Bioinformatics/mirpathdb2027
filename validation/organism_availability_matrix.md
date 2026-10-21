# Organism availability matrix (target-source coverage)

This matrix is the documented **rationale** for the production organism set.
Production organisms = **4 species**: the miRPathDB2 floor (hsa, mmu) plus two
approved native-coverage extensions — rno (rat, **validated+predicted**) and mml
(rhesus macaque, **predicted** — substantively TargetScan-predicted, with a lone
DIANA-TarBase v9 validated MTI (`mml-miR-125b-5p->ANGPT2`) **marked**, not buried).
Dog (cfa) and chicken (gga) were
**evaluated and EXCLUDED** (owner decision 2026-06-24 — an owner choice to omit,
NOT a lack of native evidence). The structured form is
`organism_availability_matrix.tsv`; this Markdown file is the human-readable
companion. Coverage values are the verified per-source coverage; if the two ever
disagree, the TSV wins.

## The locked decision (organism set)

| code | scientific name | tax id | role | basis |
|---|---|---:|---|---|
| `hsa` | Homo sapiens | 9606 | **floor** (miRPathDB2) | must be present in production |
| `mmu` | Mus musculus | 10090 | **floor** (miRPathDB2) | must be present in production |
| `rno` | Rattus norvegicus | 10116 | **approved extension (rat)** | native validated (miRTarBase) + predicted (miRDB + TargetScan) |
| `mml` | Macaca mulatta | 9544 | **approved extension (macaque, predicted)** | TargetScan -> LOW (293 miRNAs, substantive); DIANA-TarBase v9 -> 1 validated MTI (`mml-miR-125b-5p->ANGPT2`), marked into HIGH; MEDIUM empty |
| `cfa` | Canis lupus familiaris | 9615 | **evaluated, EXCLUDED (reference row)** | covered by TargetScan + miRDB, but owner chose to omit dog |
| `gga` | Gallus gallus | 9031 | **evaluated, EXCLUDED (reference row)** | covered by TargetScan + miRDB, but owner chose to omit chicken |

`hard_floor = [9606, 10090]`. `approved_extensions = [10116, 9544]` (**rat +
macaque**). `approved_removals = []` (**empty** — macaque is restored as an
approved extension, so nothing is removed). Dog (9615) and chicken (9031) appear
in **neither** allowlist — they were never in miRPathDB2 or the prototype, so
they are simply **not production organisms** (evaluated and excluded by owner
choice). miRPathDB2 stays the primary external anchor (hsa/mmu); the old
prototype baseline (hsa/mmu/mml) stays as the historical reference for the
expansion audit — and since the prototype carried mml and production re-adds it,
there is no removal to audit.

## Coverage matrix

Coverage legend: **yes** = native, usable coverage · **thin** = present but
minor/tail (not depth-bearing) · **verify** = within the source's stated species
set but the prototype fetcher does not yet ingest it — production must verify the
per-species file exists and extend the fetcher · **no** = absent / no native rows.

The **four official target sources** are miRTarBase 9.0, DIANA-TarBase v9,
TargetScan v8, and **miRDB v6.0**. miRTarBase and DIANA-TarBase are **parallel
validated resources** (verified interactions + CLIP) — both feed **HIGH +
MEDIUM**. TargetScan and miRDB are predictors in the **LOW** tier. TargetScan
covers human/mouse/rat/dog/chicken/macaque and provides macaque's predicted/LOW
evidence (293 miRNAs) — macaque's **substantive** evidence band. DIANA-TarBase v9
genuinely covers **24 species including macaque**, but for macaque its VOLUME is
**exactly one validated MTI** (`mml-miR-125b-5p->ANGPT2`): that lone interaction
is **marked** into macaque's HIGH band; macaque's MEDIUM band is empty. So macaque
is **predicted-dominant** with the single validated MTI flagged, not broadly
validated. (miRTarBase and miRDB do not cover macaque.)

| code | tax id | class | miRTarBase 9.0 | DIANA-TarBase v9 | TargetScan v8 | miRDB v6.0 | depth_class | HIGH | MEDIUM | LOW |
|---|---:|---|---|---|---|---|---|:---:|:---:|:---:|
| `hsa` | 9606 | floor | yes (validated) | yes (validated + CLIP) | yes (anchor) | yes | **full** | ✓ | ✓ | ✓ |
| `mmu` | 10090 | floor | yes (validated) | yes (validated + CLIP) | yes | yes | **full** | ✓ | ✓ | ✓ |
| `rno` | 10116 | approved_extension | yes (3rd species) | verify | yes (conserved) | yes | **validated+predicted** | ✓ | ✓ | ✓ |
| `mml` | 9544 | approved_extension | no | 1 validated MTI (`mml-miR-125b-5p->ANGPT2`) | yes (predicted/LOW, 293 miRNAs) | no | **predicted** | ✓ (1 MTI) | ✗ | ✓ |
| `cfa` *(evaluated, EXCLUDED)* | 9615 | evaluated_excluded | thin (tail) | no (owner-omitted) | yes | yes | **n/a (excluded)** | ✗ | ✗ | ✗ |
| `gga` *(evaluated, EXCLUDED)* | 9031 | evaluated_excluded | thin (tail) | no (owner-omitted) | yes | yes | **n/a (excluded)** | ✗ | ✗ | ✗ |

Depth classes (production only): `hsa`/`mmu` = **full** (validated + CLIP + 2
predictors); `rno` = **validated+predicted**; `mml` = **predicted** (predicted-
dominant: TargetScan -> LOW for 293 miRNAs; DIANA-TarBase v9 -> HIGH = the single
marked validated MTI `mml-miR-125b-5p->ANGPT2` only; MEDIUM empty).

## Confidence tiers each production organism can actually populate

The user-facing HIGH/MEDIUM/LOW rollup is auditable back to source + method (see
`provenance/manifest.yaml` `target_sources.*.tier_rollup`):

- **HIGH** = miRTarBase strong/functional (reporter assay / Western blot / qPCR) +
  DIANA-TarBase verified-strong → internal `val_strong`.
- **MEDIUM** = miRTarBase weak (other experimental) + DIANA-TarBase CLIP/weaker
  experimental → internal `val_weak`, `clip_supported`.
- **LOW** = TargetScan + miRDB predicted (union; intersection is a stronger
  `low+`) → internal `pred_union`, `pred_intersection`.

(Exact verified-strong vs CLIP/weaker tier boundaries: **confirm Tuesday**.)

Per-organism consequence (production set):

- `hsa`, `mmu` — can populate **all three** tiers (miRTarBase + DIANA-TarBase feed
  HIGH/MEDIUM; TargetScan + miRDB feed LOW).
- `rno` — can populate all three; HIGH/MEDIUM are thinner than hsa/mmu (miRTarBase
  is the 3rd-largest species; DIANA pending fetcher extension), LOW is full via
  TargetScan + miRDB.
- `mml` — is **predicted-dominant**: TargetScan feeds a substantive **LOW** band
  (293 miRNAs). Its **HIGH** band is a SINGLE marked validated MTI
  (`mml-miR-125b-5p->ANGPT2`, DIANA-TarBase v9); its **MEDIUM** band is **empty**.
  miRTarBase and miRDB do not cover macaque.

## Why dog (`cfa`) and chicken (`gga`) were excluded

Dog and chicken were **evaluated** for the production set and **excluded** (owner
decision 2026-06-24). This is an **owner choice to omit**, NOT a lack of native
evidence — both ARE covered by native target sources:

- **TargetScan v8** — covers both dog and chicken (TargetScan is broad-vertebrate).
- **miRDB v6.0** — ships both dog and chicken.
- **miRTarBase 9.0** — only a thin/tail per-species set for both (not depth-bearing).
- **DIANA-TarBase v9** — contributes no HIGH/MEDIUM evidence for dog/chicken here
  (they are owner-omitted).

So dog and chicken DO have native predicted coverage (TargetScan + miRDB), but the
owner chose to **exclude** them from the production set. They were **never** in
miRPathDB2 nor in the v3 prototype, so they are in **neither** the
`approved_extensions` nor the `approved_removals` allowlist — they are simply not
production organisms. The reference rows above retain them so the exclusion (and
its rationale) is auditable.

## DIANA-TarBase v9 (validated + CLIP; 24 species — but macaque volume = 1 MTI)

DIANA-TarBase v9 is a **parallel validated resource to miRTarBase** — verified
interactions **and** CLIP/AGO support, **not** a CLIP-only source. It therefore
feeds **both HIGH and MEDIUM**. The resource **genuinely covers 24 species,
including macaque**. The rat + macaque fetcher extension has **landed** (both load
from the `Other_species_TarBase-v9.tsv.gz` bundle, sliced by the `species`
column), and the production build revealed the real per-organism VOLUME:

- For `rno` (rat) the bundle carries **~1125 validated MTIs across 26 miRNAs** —
  substantive validated evidence feeding rat's HIGH/MEDIUM. (rat is therefore
  `validated+predicted`.)
- For `mml` (macaque) the bundle carries **EXACTLY ONE validated MTI**:
  `mml-miR-125b-5p -> ANGPT2` (a Luciferase Reporter Assay). That single
  interaction is **marked** into macaque's **HIGH** band; macaque's **MEDIUM** band
  is **empty**. This is the real extent of macaque's experimentally-validated
  coverage — so macaque is **predicted-dominant** (its substantive input is the
  TargetScan-predicted LOW band), with the lone validated MTI flagged rather than
  buried under a "validated" label.

## Pathway-source coverage (GO / Reactome / WikiPathways)

Target-source coverage (above) drives the HIGH/MEDIUM/LOW input miRNAs; the
**pathway** sources are separate (see `config/pathway_sources.yaml`). Per-organism
pathway coverage (owner decision 2026-06-26):

| code | tax id | GO | Reactome | WikiPathways |
|---|---:|:---:|:---:|:---:|
| `hsa` | 9606 | native | native | native |
| `mmu` | 10090 | native | native | native |
| `rno` | 10116 | native | native | native |
| `mml` | 9544 | native | **INFERRED** (orthology-projected, limited) | **NONE** |

**Macaque Reactome is INFERRED, not curated.** Reactome ships **0 curated
`R-MML-` pathways** (confirmed: `reactome.db` v1.95.0 and
`Ensembl2Reactome_All_Levels.txt` both carry 0 macaque rows). Reactome itself,
however, lists *Macaca mulatta* as an **orthology-inferred** species. The
production build therefore replicates Reactome's own inference: every human
`R-HSA-` pathway gene set is projected onto macaque via the **NCBI
`gene_orthologs` 9606↔9544 strict 1:1 ortholog map** (17,586 pairs), mapped to
macaque SYMBOL through `org.Mmu.eg.db`, and relabelled `R-HSA-<n>` → `R-MML-<n>`.
Result: **2,777 inferred macaque Reactome pathways** over a **10,237-symbol**
macaque universe, marked INFERRED + LIMITED in
`provenance/macaque_inferred_reactome_manifest.json` (ortholog source pinned by
size + sha256). The inferred Reactome ORA matrices are
`data/processed/enrichment/mml/prod_REACTOME_{HIGH,MEDIUM,LOW}.csv` (LOW = 293
TargetScan miRNAs × 2,642 scored pathways; HIGH/MEDIUM = the lone validated MTI).

**Macaque WikiPathways is NONE.** WikiPathways genuinely ships **no
`Macaca_mulatta` collection** (no macaque GMT in the monthly release; *Macaca
mulatta* is not in `rWikiPathways::listOrganisms`). Unlike Reactome, WikiPathways
offers no orthology-inference path, so macaque WikiPathways stays NONE.

## Macaque (`mml`) — a predicted production organism (lone validated MTI marked)

`mml` (Macaca mulatta, 9544) is a **production organism** (no longer a removal).
Its depth is **predicted** (predicted-dominant): its **substantive** evidence is
the **TargetScan v8** predicted **LOW** band (293 miRNAs). Its only validated
evidence is a **single** DIANA-TarBase v9 MTI — `mml-miR-125b-5p -> ANGPT2` —
which is **marked** into the HIGH band (MEDIUM is empty). It has **no** miRTarBase
per-species file and is **not** shipped by miRDB v6.0 (5 species:
hsa/mmu/rat/dog/chicken), so those two sources contribute nothing for macaque.
Because the prototype baseline already
carried mml and production re-adds it, the parity gate records macaque as an
**approved extension** (predicted; lone validated MTI marked), not a removal —
`approved_removals` is empty.
