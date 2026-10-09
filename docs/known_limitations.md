# Historical and reproducibility limitations

- **SOURCE NOT FOUND:** no standalone producer was recovered for the repaired clinical metadata Seurat object required before balanced cNMF.
- **SOURCE NOT FOUND:** no standalone historical exporter was recovered for `34A5C_FROZEN_HUMAN_TOP100.csv`.
- **No frozen result/input artifacts are distributed in this code-only release.** Previously staged derived tables, gene/weight checkpoints and figure SVGs are withheld pending separate author/provider review. Two mixed-source 26-signature catalogs were already withheld for unresolved redistribution rights.
- The regulator and glial scripts may depend on archived latent factors/consensus checkpoints not included here. Historical figure rendering relies on donor/spot aggregates, original SVG assets and other missing outputs.
- The recorded R, cNMF and regulator software environments are historical evidence, not a fully tested, installable multi-language environment lock. External H5AD and figure environments remain partly undocumented.
- Timestamp selectors and fixed checkpoint directory names may prevent a fresh run from reproducing archived output paths. Do not silently replace original checkpoints or claim exact re-execution.
- Full end-to-end code execution and independent scientific replication are **NOT TESTED**. The original 52 scientific scripts have not been modified for this code-only preparation.

See `source_gaps.md`, `../RUN_ORDER.md`, `../metadata/source_manifest.tsv` and `../metadata/script_to_result_map.tsv` for precise stage evidence and provenance.
