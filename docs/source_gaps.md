# Historical source gaps

Two independent historical producer steps remain **SOURCE NOT FOUND**. This is a provenance boundary; neither a recovered implementation nor an independent complete rerun is claimed.

| Boundary | Located evidence and retained input | Remaining limitation |
|---|---|---|
| `human_harmony_r0p4_level1_annotated_with_clinical_metadata_repaired.rds` | Archived repaired global Seurat checkpoint exists between early global annotation and balanced cNMF. It contains human data and is not distributed here. | Its complete independent clinical-metadata repair producer was not located. The early global annotation script and a myeloid-specific repair script do not establish this missing transformation. |
| `34A5C_FROZEN_HUMAN_TOP100.csv` | The exact frozen human P5 gene/weight checkpoint is retained with original SHA256. The approved stage compared it at file level with final manuscript Top100: identical ordered membership of 100 genes, with four last-printed-digit weight differences of at most 1e-14. Both originals were preserved separately, without substitution. | Its standalone historical export script was not located. External scoring can consume the checkpoint, but its original export step cannot be independently reconstructed from a recovered script. |

These gaps do not alter archived scientific results. No new exporter or metadata repair is invented. Neither gap is evidence that the P5 biological finding is invalid. Full scientific execution is NOT TESTED. See `metadata/frozen_gene_checkpoint_integrity.json`, `metadata/source_manifest.tsv` and `REPOSITORY_UNRESOLVED.md`.
