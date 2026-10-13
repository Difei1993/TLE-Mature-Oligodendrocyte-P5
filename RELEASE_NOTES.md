# P5 analysis and supplementary display code v1.0.1

This version adds the three display-only helpers for Supplementary Figures S3–S5 in `supplementary_figures/S3_S5/`, with external path configuration, recorded display dependencies, panel-source mappings and script provenance. The original 52 scientific scripts and the three reviewed display scripts are unchanged by the metadata repair.

README and publication notes now correctly identify the formally released v1.0.0 and its Zenodo version DOI 10.5281/zenodo.23266083. CITATION.cff identifies version 1.0.1. The v1.0.0 release and DOI remain available as the prior version; that DOI must not be presented as the DOI of v1.0.1.

The added helpers reuse archived image layers, saved UMAP coordinates/manual annotations and frozen tables. No scientific analysis, figure algorithm, parameter, gene list, program ranking, statistical effect or FDR was changed or recomputed. The repository does not redistribute human expression matrices, Seurat/H5AD objects, private cell-level coordinate exports, third-party catalogs or internal review packages. The MIT license covers author-owned software, not excluded third-party data.

Syntax, dependency and input-availability checks are documented in the display-code provenance. Independent full scientific reproduction and cross-platform redraw remain NOT TESTED; the public code-only distribution lacks several historical inputs and producers needed for complete regeneration. The metadata-repair task did not run the scientific or display scripts.

A new version-specific Zenodo DOI is not supplied in these prepared notes. Record it only after a real deposit is created and the archived version is verified. Do not reuse the prior version DOI as the identifier of this version.
