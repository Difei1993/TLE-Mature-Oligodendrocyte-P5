# Run order for the original P5 analysis chain

These are documented scientific execution instructions, not a record of reruns performed during repository preparation. **Every analysis command below is NOT TESTED.** No analysis runner was invoked. Supply provider inputs and archived checkpoints in a NEW reconstruction work tree, set P5_WORK_ROOT (and TLE_PROJECT, if used) to that tree, and change the working directory there. Restore the cNMF executable environment at `04_tables/python_envs/cnmf310/Scripts/`. Use the relevant historical R/Python environment, not the current syntax checker.

## Candidate input availability

Two original mixed-source signature catalogs are withheld for unresolved redistribution rights; see `metadata/excluded_input_manifest.tsv` and [third-party notices](docs/third_party_notices.md). Catalog-dependent stages require lawful restoration of those exact originals; they cannot be executed from this candidate alone. The independent metadata-repair producer and frozen-human-Top100 exporter remain SOURCE NOT FOUND, as described in [source gaps](docs/source_gaps.md).

This code-only release does not contain the 35 previously staged frozen input files, the three SVG aliases, or the copy-only input-restoration helper. Exact historical inputs must be lawfully obtained outside the repository. No scientific script was changed.

## Frozen checkpoint binding

The archived balanced input is 32B0_BALANCED_CNMf_INPUTS_20260728_215801; primary/sensitivity H5AD work trees are 32B_CNMf_WORK_20260728_215801 and 32B2A_CNMf_SENSITIVITY_20260728_215801; final factors are 32B2B_FINAL_CNMf_20260728_215801. The P5/program/statistics freeze is 32B3A_UNBIASED_PROGRAM_FREEZE_20260729_010256. Later content/catalog/pathway/classification/34A score, donor, resource and regulator freezes must match `metadata/archived_evidence.tsv` and the literal paths in the sources.

Fresh stages often generate current timestamps; fixed downstream bindings and `latest()` selectors do not guarantee that a fresh run will consume the intended checkpoint. Bind reconstructed outputs explicitly to the archived layout in a separate work tree before proceeding. Do not silently mix new archives with old freezes and do not describe this unresolved binding as turnkey reproduction.

## Stage sequence

| Stage | Scripts in execution order | Required input and endpoint | Figure/table |
|---|---|---|---|
| 01 Discovery compartment | 01 → 02 → 03 → 04 → 05 → 07 → 09 → 11 → 12 → 13 → 14 → 15 (full filenames in source manifest) | GSE190452 processed counts → global Harmony/Level1 annotation. Step 04 writes QC flags CSV, not an RDS; 05 reads raw metrics plus that CSV. Step 09 writes human_qc_singlets.rds. | Fig. 1; S1 |
| Metadata checkpoint boundary | **SOURCE NOT FOUND** | Supply archived repaired clinical-metadata annotated global Seurat object before 32B0; original global cluster 0 contains 15,947 mature OL nuclei. Do not substitute a later lineage reclustering or the unrelated myeloid repair script. | Fig. 1; S1 |
| 02 Balanced cNMF/P5 | 32B0 R → 32B0 H5AD converter → 32B1 H5AD pilot → 32B2A sensitivity converter → 32B2A K confirmation → 32B2B final100 → 32B2C filtered audit | Four glial lineage inputs; two independently balanced seeds; mature OL 450/donor, K=6. Unfiltered density=2.0 in 32B2B is diagnostic. P5 freezing uses the filtered plan (OL threshold 0.01), not the unfiltered diagnostic spectrum. | Fig. 1; S2/S3 |
| 03 Donor statistics / safe freeze | 32B3A → 32B3B0; then 32B3B1A/B (functional stage); then 32B3C → 34A0 | Frozen spectra/usage/H5AD metadata → matched programs, Top100/full spectrum, 21-program family effects/LODO. BH remains across 7 astrocyte+5 myeloid+3 OPC+6 OL primary programs. | Fig. 1; S2/S3/S4 |
| 04 External contextualization | SEA-AD extractor → SEA-AD R exporter; 35C_01 → 35C_03 → 35C_04 → 35C_05 → 35C_06 → 35C_07 → 35C_08B | Provider H5AD raw.X and GSE140393 counts/manual cluster labels; 35C_03 and 35C_04 use the original 35C_01 candidate object directly. Superseded 35C_02 is not needed. Adult TLE/SEA-AD mature OL scores and separate OLIG2-positive nuclear context. | Fig. 2,S1; S1/S5 |
| 05 Functional annotation / P5 evidence | 32B3B1A → 32B3B1B → 32B3C/34A0 as above; then P5_RESTART1 → P5_RESTART2 | 26-signature catalog, frozen g:Profiler response/background, P5 Top100/full spectrum and external outputs. RESTART1 needs the preceding external archive; RESTART2 consumes RESTART1 plus catalogs. | Fig. 3; S2/S4/S5 |
| 06 Regulators | 34A1B2 → 34A2 → 34A2B → 34A3 → P5_RESTART3 | Primary H5AD, archived contrastiveVI latent/consensus prerequisites, gene-score weights, donor scores, dated CollecTRI. The unreported contrastiveVI training is not staged. Targets overlapping P5 Top100 are excluded in primary ULM; within-condition permutation and network-specific FDR remain unchanged. | Fig. S1D; S6/S7 |
| 07 Glial context / LR | P5_RESTART3_FIX1 → P5_RESTART4_FIX8_FAST | FIX1 donor composite correlations supersede the older context output from base RESTART3. FIX8 embeds its own H5AD exporter, then uses TMM/log2 CPM and frozen LR/ligand-target priors. | Fig. 4; S6/S8 |
| 08 Spatial reference | P5_RESTART5_FIX2_REGENERATED | Four GSE295237 matrices/positions, frozen P5 Top100 and LR shortlist; produce descriptive within-sample residual/neighborhood correlations and retained negative qualification. | Fig. S2; S9 |
| Evidence gate | P5_RESTART6 | Frozen RESTART1–5 outputs, after FIX1/FIX8/FIX2 effective results exist. | Manuscript evidence freeze |
| 09 Data-driven figures | draw_final_polished.R → redraw_figure3_typography.R | Frozen manuscript tables plus excluded donor/spot tables and archived Figure4C_donor_coordinates.tsv. This coordinate input is produced in the earlier draw_p5_figures.R; its provenance is documented, but obsolete full plotting code is not needed for the retained final figures. | Figs. 1–4,S1,S2 |
| 09 Final formatting/export | format_submission_figures.R → render_final_technical_figures.R → finalize_final_technical_svg_attributes.R; then render_p5_review_figures.R | The three layers of unchanged archived SVG inputs are represented by 15 distinct files and three byte-identical restoration aliases. Main V3 TIFF/EPS exporter uses final submission SVGs. Historical review output folder name is retained; no V3 file was overwritten or re-exported. | V3 final figures |
| Archived PDF box patch | finalize_final_technical_pdf_boxes.R | **Not a generic rerun stage:** exact archived PDFs, page dictionaries/object IDs/xref offsets and dimension metadata required. Never apply old offsets to newly rendered PDFs. | Existing Figure1/S1/S2 PDF technical formatting |

**Ordering hazard:** never run base RESTART3 after FIX1, because it would overwrite the final sender/context shortlist. Do not run obsolete LR FIX1–FIX7 or the obsolete nuclear bridge. No new tests or scientific recalculations are part of this repository preparation.

## Executable entry patterns

After prerequisites and historical environments have been restored, individual script entry patterns are `Rscript <original-relative-script-path>`, `<historical-python> <original-relative-script-path>`, or `powershell -File <original-relative-PS1-path>`. The staged helper copies these original logical locations. Use `metadata/source_manifest.tsv` for every exact filename, input and output. Extraction/export scripts have command-line arguments; their original argparse definitions are authoritative. No commands here have been executed during staging.

## Endpoint example for Figure 1

Start with GSE190452 preprocessing/annotation and the repaired metadata checkpoint, generate both balanced cNMF inputs, obtain the filtered K=6 matched program factors, and run 32B3A. Its frozen donor-level effect table contains CLR/Hedges g/exact allocation/BH21/LODO. RESTART1 retains P5 evidence; `S3_all_discovery_effects.tsv` and `S3_P5_effects.tsv` are archived manuscript bindings, not newly recomputed tables. The final draw script additionally needs donor usage tables that are excluded from redistribution. Static source recovery therefore does not demonstrate complete regeneration of Figure 1.
