# miRPathDB 2027: pipeline of release v2027.1

miRPathDB links microRNAs to the pathways their target genes are enriched in. Release 2027.1 covers human, mouse,
rat and rhesus macaque. Its miRNA target genes come from miRTarBase 9.0, DIANA-TarBase v9, TargetScan 8 and miRDB
6.0, in cumulative High, Medium and Low evidence tiers. Each miRNA's targets are tested for over-representation in
Gene Ontology, Reactome and WikiPathways terms with clusterProfiler, with Benjamini–Hochberg adjustment over the full
pathway family. The database is online at https://web.ccb.uni-saarland.de/mirpathdb/.

This repository holds the analysis pipeline of the data release at tag `v2027.1` (commit `aeef54ba64aa818a04237fad837f9ae26b19ee69`), without version history.
The web application that serves https://web.ccb.uni-saarland.de/mirpathdb/ is not part of this repository.
The scripts of the paper's cross-species comparison run inside the build environment and are available from the authors on request.

| Path | Contents |
|---|---|
| `src/`, `scripts/`, `tests/` | the pipeline code and its tests |
| `config.yaml`, `config/`, `pyproject.toml` | configuration and dependencies |
| `provenance/` | the R environment (`environment.mpd_cp.yml`) and the production manifests |
| `release-manifests/2027.1/` | the checksum manifests of data release 2027.1 |
| `docs/reproducibility/checksums/` | the checksum manifests of the previous package, whose 57 enrichment strata are identical to those of release 2027.1 |
| `validation/`, `data/` | organism reference tables and WikiPathways name repairs that the pipeline reads |

The release data and the source files are not part of this repository. The data are downloadable from the site;
the source files have to be obtained from their providers (see "Rebuild one stratum").

## Download and verify release 2027.1

1. From the Download page of https://web.ccb.uni-saarland.de/mirpathdb/, download the files of release 2027.1
   (`enrichment_retained.tar.gz`, `pathway_labels.tsv`, `v2_reproduction_metrics.json` and `checksums.sha256`) into one folder.
2. Check them against the page's checksum list:

   ```bash
   sha256sum -c checksums.sha256
   ```

3. `release-manifests/2027.1/` holds the checksum manifests of the complete release package, which will be available for
   download at publication. The package's `RELEASE_INFO.json` names this pipeline commit, `aeef54ba64aa818a04237fad837f9ae26b19ee69`.

## Rebuild one stratum

A stratum is one table `enrichment_significant/<organism>/<FAMILY>_<TIER>.tsv` of the release. The steps below
recompute the rat GO Biological Process High stratum (13,560 rows) and compare it with the released file. They take
about five minutes on 8 cores. Three runners in `scripts/` compute the 57 strata: `run_production_pathways.py`
(36 strata: human and mouse, and rat Reactome and WikiPathways), `run_production_enrichment.py` (18: rat and macaque
GO) and `run_macaque_inferred_reactome.py` (3: macaque Reactome).

1. **Environments.** Python 3.11 with the pipeline's dependencies (`pip install -e '.[validation]'`), and R 4.5.3
   with the Bioconductor packages pinned in `provenance/environment.mpd_cp.yml`.
2. **Sources.** Obtain the source files the runners read (miRTarBase 9.0, DIANA-TarBase v9, TargetScan 8, miRDB
   6.0, and the NCBI and pathway files) and put them into one directory. The file names the runners expect are set
   in `src/mirpathdb_prod/sources/production_target_sets.py` and the source adapters next to it. The files are not
   part of this repository, and neither is the manifest that pins them by checksum, so the runner cannot verify
   them; you supply them.
3. **Paths.** At the repository root, write a `config.local.yaml` with your own paths. It overrides `config.yaml`
   and keeps the outputs out of the repository:

   ```yaml
   paths:
     data_processed_dir: <processed output directory>
     provenance_dir: <provenance output directory>
     source_cache_dir: <source directory>
   ```

   Relative paths in both files resolve against the repository root. Without an override, the outputs go to
   `data/processed/`. The location of the miRPathDB 2.0 archive, which only the 2.0 comparison code reads, defaults
   to `data/raw/mirpathdb2/version_2`.

4. **Run**, from the repository root:

   ```bash
   PYTHONPATH=src MPD_RSCRIPT=<R 4.5.3 environment>/bin/Rscript MPD_MAX_CORES=8 \
       python scripts/run_production_enrichment.py --organisms rno --tiers HIGH
   ```

5. **Compare**, keyed on miRNA and pathway: integers must be equal, and P-values must agree to a relative 1e-9.

   ```bash
   python scripts/reproducibility/compare_stratum.py \
       --rerun <processed output directory>/enrichment_significant/rno/GO_BP_HIGH.tsv \
       --package <release 2027.1 folder> --relpath enrichment_significant/rno/GO_BP_HIGH.tsv \
       --sums-copy release-manifests/2027.1/SHA256SUMS \
       --expected-rows 13560
   ```

   Exit status 0 means the stratum is reproduced.
   Run this way from this tag, one stratum per computation path reproduced the corresponding released tables: integers exactly, P-values within a relative tolerance of 1e-9.

## Tests

`pip install -e '.[dev]'`, then `pytest` from the repository root.

## Cite

If you use miRPathDB, please cite the miRPathDB 2027 article in Nucleic Acids Research. `CITATION.cff` holds the
citation metadata of this code.

## Licence

The code is released under the MIT licence (`LICENSE`).

## Use of AI assistance

The code was developed and tested with Claude (Anthropic), used through Claude Code, which also helped verify the data, as stated in the paper. The authors made all scientific and curatorial decisions, checked AI-assisted outputs against source data, code and literature, and take responsibility for the code.
