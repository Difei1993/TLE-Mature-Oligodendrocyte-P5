# TLE Mature Oligodendrocyte P5 — analysis code

This repository contains the **52 author-written scientific analysis and figure scripts** associated with the manuscript **“Donor-aware characterization of a mature oligodendrocyte remodeling program in human temporal lobe epilepsy”** (planned submission to *Molecular Neurobiology*). The scientific reference is the frozen 2026-10-08 V3 manuscript package.

## Scientific scope and interpretation

In the GSE190452 human neocortical dataset (four TLE and four control donors), P5 was characterized as a mature oligodendrocyte transcriptional program associated with myelin-related, lysosomal/proteostasis and lipid processing functions. The donor-level directional association **did not survive multiple testing** across the 21-program discovery family (**FDR = 0.60**). GSE140393 and SEA-AD offer cross-dataset context but not an independent controlled disease replication. GSE295237 consists of FCD I/control temporal neocortical spatial samples, **not TLE samples**; no consistent spatial reference support was observed. APOE–LRP1 remains an exploratory association, **not evidence for causal signaling**. Senescence-related overlap does not establish senescent oligodendrocytes.

## Repository contents

- `scripts/`: all **52 original scientific scripts**, unchanged from the reviewed public candidate (R, Python and PowerShell).
- `RUN_ORDER.md`: stage order and historical checkpoint boundaries.
- `metadata/source_manifest.tsv` and `metadata/script_to_result_map.tsv`: script provenance, original and candidate SHA256, documented inputs and archived Figure/table associations.
- `metadata/datasets.tsv` and `docs/data_inputs.md`: official data acquisition and input requirements.
- `environment/`: available historical versions and recovered session records (not a complete reproducible lockfile).
- `docs/source_gaps.md`, `docs/known_limitations.md`: explicit reconstruction limitations.
- `LICENSE`: MIT license for the author-owned scripts; **does not license external datasets, gene-set catalogs, regulator priors or third-party content**.

**No scientific input/result datasets or figure SVG assets are redistributed in this code-only package.** This includes 35 frozen small results/figures from a private pre-publication candidate, two excluded mixed-source signature catalogs, human expression matrices, SEA-AD H5ADs and donor/spot-level records. Original input/checkpoint names in manifests describe historical dependencies, *not* files available in this repository. Original hashes for unbundled inputs are recorded where previously available; users must lawfully obtain or reconstruct matching inputs. This is intentional, not a claim of turnkey reproduction.

## Public data origins

| Dataset | Role | Official access |
|---|---|---|
| GSE190452 | Human neocortex discovery, 4 TLE and 4 controls | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190452) |
| GSE140393 | External human TLE cell-context study | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140393) |
| SEA-AD | Adult MTG oligodendrocyte reference | [SEA-AD](https://brain-map.org/consortia/sea-ad/our-data) |
| GSE295237 | Four FCD I/control spatial reference samples | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295237) |

Follow each provider's access, citation and reuse conditions. Access to data does not itself authorize public redistribution of copies or derivatives.

## Historical environments and execution

Historical evidence includes **R 4.4.2**, **Python 3.10.20 / cNMF 1.7.1**, and a **Python 3.12.10 / decoupler 2.2.0** regulator environment. These are separate stage environments. Other stages and historical transitive dependencies are not fully specified. Consult `environment/README.md` and `RUN_ORDER.md` before any re-execution. The original code may require dated checkpoint objects and specific directory names; do not substitute newly regenerated inputs without verifying alignment to the frozen analysis.

### Important limits

1. The original producer for `human_harmony_r0p4_level1_annotated_with_clinical_metadata_repaired.rds` was **not found**. The historical checkpoint exists outside this code-only distribution.
2. The independent original exporter for `34A5C_FROZEN_HUMAN_TOP100.csv` was **not found**. The archived Top100 input was verified in the private stage but is **not distributed** here.
3. Several regulator/figure stages depend on archived latent states, consensus metadata, donor/spot results and fixed timestamp bindings unavailable in this code-only package.
4. **No scientific pipeline was re-executed when preparing this repository.** Script/hash checks and documentation audits are not proof of independent end-to-end reproducibility.

See [historical source gaps](docs/source_gaps.md) and [limitations](docs/known_limitations.md). The absence of a full one-command rerun should not be taken to change the reported archived study results.

## License and citation

**Copyright © 2026 Difei Wang.** The original research code in `scripts/` is licensed under the [MIT License](LICENSE). Source-data and third-party materials retain their respective terms. Cite the companion manuscript when it is published and the software repository using [`CITATION.cff`](CITATION.cff) if appropriate. No manuscript DOI, software DOI, Zenodo deposit or formal GitHub release is claimed at preparation time.
