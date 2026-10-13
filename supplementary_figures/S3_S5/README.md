# Display-only code for Supplementary Figures S3–S5

Version 1.0.1 adds the display helpers described below. The existing v1.0.0 analysis release and its version DOI remain available as the prior version. This directory adds three author-owned display helpers; it does not modify or replace the 52 scientific scripts under `scripts/`.

## What the scripts actually do

| Script | Scope |
|---|---|
| `build_p5_final_supplementary_figures.py` | Recompose archived plot layers, plot already frozen descriptive tables and stored UMAP coordinates, and export the eight S3–S5 artwork leaves as PDF/SVG/600dpi RGB PNG/TIFF plus individual S3/S4/S5 PDFs. |
| `p5_display_image_layers.py` | Locate existing image axes, tick positions and legend layers for the builder. Its optional CLI only inspects geometry of three archived DotPlots. It never recomputes marker expression. |
| `export_p5_saved_external_umap.R` | Read saved GSE140393 metadata/UMAP slots from two existing RDS objects and export matching labels and exact little-endian double coordinates for private local drawing. It neither calculates UMAP nor changes annotations. |

No QC, normalization, Harmony, clustering, cNMF fitting, differential expression, enrichment, NicheNet, spatial analysis or statistical testing is performed here. P5 FDR remains 0.60. The six Top8 lists are read in frozen order, not selected afresh. The historical pilot threshold 0.5 is distinct from final consensus threshold 0.01; no K selection is rerun.

## Archived inputs are not distributed

`provenance/required_inputs.tsv` lists all 17 files the builder requires. Names and hashes describe historical local dependencies; they are not downloads offered by this repository. The official data-origin links in the root repository explain provider access, but do not supply matching historical checkpoints, derived TSVs or archived figures.

| Panels | Required evidence and limitation |
|---|---|
| S3a–d | Four exact archived GSE190452 UMAP/DotPlot PNGs, including the r=0.4 panel in a historical grid. Pixel crops assume those original layouts; arbitrary new figures are not interchangeable inputs. |
| S3e | Frozen donor-count TSV derived from Table S1; it describes availability and the two 450-per-donor sampling schemes, not disease-related cell abundance. |
| S4a–c | Saved primary K4–14 and sensitivity K6/K7 statistics in one archived TSV. No model is fitted or optimized by this script. |
| S4d | Frozen 48-row, six-program Top8 TSV. Neuronal-associated genes are retained. No gene catalog or substitute gene ranking is bundled. |
| S4e/f | Two exact final K6/dt0.01 consensus-diagnostic PNGs. These are component spectra, not new cell-cluster diagnostics. |
| S4g | Frozen 21-program primary effects/FDR TSV derived from Table S3. Zero-to-effect segments are not confidence intervals, and CLR effects are not absolute expression changes. |
| S5a/b | Two sets of private saved-label CSVs, binary double coordinates and palette CSVs, read from the already annotated objects if lawfully available. They retain 3,666 candidate/2,138 oligodendroglial cells and the original categories. These include cell keys and are not approved here for public redistribution. |
| S5c | Original marker DotPlot PNG with frozen percent-expressed and scaled-mean-expression layers. No marker averages or detection percentages are recomputed. |

S5 represents the five GSE140393 whole-cell TLE samples, not a joint adult TLE/SEA-AD UMAP or independent controlled disease replication. Original input versions, label/barcode alignment and permitted reuse remain the user's responsibility. Do not substitute unverified regenerated objects and describe the output as the frozen study result.

These three helpers do not independently reconstruct the full scientific analysis, recover missing historical producers, or assemble the main manuscript/ESM1 cover and caption pages. A public repository checkout alone cannot fully redraw S3–S5 because its required frozen inputs are deliberately absent. Full independent scientific reproduction is **NOT TESTED**.

## Recorded display environment

`DISPLAY_ENVIRONMENT.json` preserves the runtime recorded for the approved figure preparation: Python 3.12.14, matplotlib 3.10.9, NumPy 2.3.5, Pillow 12.3.0, SciPy 1.18.0 and pypdf 6.10.0; the read-only R conversion used R 4.6.1. `requirements-display-recorded.txt` records only the five Python packages imported by these helpers. This is a recorded display environment, not a portable transitive lockfile or the historical biological-analysis environment.

The R reader needs `methods` and compatible saved-object classes/packages; class/environment compatibility outside the original local installation has not been tested. Arial is requested by the unchanged drawing code. Font binaries are not redistributed. A different installed font or dependency implementation may change typography or image output; cross-platform and pixel-identical redraws are **NOT TESTED**.

## Configure paths without local-machine paths in source

Copy `paths.example.json` to `paths.local.json` and set four paths to your lawful local archives and a new empty output directory. JSON-relative paths resolve against the configuration file's directory. A matching CLI option overrides its JSON field and is resolved relative to the current working directory. Configuration may live outside the repository; do not publish your machine-specific configuration.

Check availability without opening scientific data or creating outputs:

```sh
python supplementary_figures/S3_S5/build_p5_final_supplementary_figures.py --config /path/to/paths.local.json --check-inputs-only
```

When authorized and matching frozen inputs exist, draw only the archived evidence:

```sh
python supplementary_figures/S3_S5/build_p5_final_supplementary_figures.py --config /path/to/paths.local.json
```

All four fields may instead be provided directly:

```sh
python supplementary_figures/S3_S5/build_p5_final_supplementary_figures.py --analysis-root /path/to/archive --frozen-table-dir /path/to/display-tables --umap-export-dir /path/to/private-exports --output-dir /path/to/new-output
```

If coordinate exports do not exist but the exact original RDS objects are lawfully available, first check their presence; this does not read the objects:

```sh
Rscript supplementary_figures/S3_S5/export_p5_saved_external_umap.R --objects-dir /path/to/frozen-objects --output-dir /path/to/new-private-exports --check-inputs-only
```

Remove `--check-inputs-only` only when intending the private read-only conversion. It writes sensitive cell-key/coordinate exports and must use a new empty local directory. The builder then uses that directory as `umap_export_dir`. Keep both scientific inputs and generated output directories outside the public repository. Existing non-empty outputs are refused to prevent overwriting frozen/submission evidence. Input availability checks do not establish content validity or scientific reproducibility.

## Provenance, execution and license

`provenance/source_manifest.tsv` records approved original script SHA256s and the adapted SHA256s. The public diffs remove original absolute local paths for privacy and are information-only rather than apply-ready patches; exact byte-preserving originals and full diffs remain in the author's local review area. `figure_source_map.tsv` maps each panel to its original generating code/archive version without exposing machine paths. `execution_status.tsv` distinguishes the original approved local display execution from candidate syntax, help and availability checks; candidate rendering and biological analysis are **NOT TESTED**.

The repository's existing [MIT license](../../LICENSE) covers these author-owned scripts and authored documentation, preserving its copyright/license notice. Importing third-party packages does not redistribute their code; those packages retain their licenses. The MIT license does not apply to provider data, Seurat/H5AD objects, cell-level exports, third-party catalogs, archived input images/databases or fonts. No such inputs are included here. See `provenance/excluded_materials.tsv` and the existing [third-party notices](../../docs/third_party_notices.md).
