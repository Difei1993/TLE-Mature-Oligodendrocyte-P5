###############################################################################
# 32B-3C
# Manuscript-safe labels, claim boundaries and validation priorities
#
# No new statistical test.
# No new enrichment.
# No program reselection.
###############################################################################

from __future__ import annotations

import hashlib
import json
import zipfile
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


###############################################################################
# 0. Paths
###############################################################################

PROJECT = Path(
    __import__("os").environ["P5_WORK_ROOT"]
)

RELEASE_ROOT = PROJECT / "99_frozen_releases"

STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")

OUTPUT_ROOT = (
    RELEASE_ROOT
    / f"32B3C_MANUSCRIPT_SAFE_PROGRAM_FREEZE_{STAMP}"
)

OUTPUT_ROOT.mkdir(
    parents=True,
    exist_ok=False,
)

(OUTPUT_ROOT / "tables").mkdir()
(OUTPUT_ROOT / "figures").mkdir()
(OUTPUT_ROOT / "manuscript_text").mkdir()
(OUTPUT_ROOT / "freeze").mkdir()
(OUTPUT_ROOT / "audit").mkdir()


###############################################################################
# 1. Locate approved inputs
###############################################################################

def latest_directory(prefix: str) -> Path:
    candidates = [
        path
        for path in RELEASE_ROOT.glob(f"{prefix}*")
        if path.is_dir()
    ]

    if not candidates:
        raise FileNotFoundError(
            f"No release directory found for: {prefix}"
        )

    return max(
        candidates,
        key=lambda path: path.stat().st_mtime,
    )


B3A_ROOT = latest_directory(
    "32B3A_UNBIASED_PROGRAM_FREEZE_"
)

B3B0_ROOT = latest_directory(
    "32B3B0_PROGRAM_CONTENT_AUDIT_"
)

B3B1B_ROOT = latest_directory(
    "32B3B1B_PATHWAY_AND_26E_INTERPRETATION_"
)


###############################################################################
# 2. Required files and gates
###############################################################################

INTERPRETATION_PATH = (
    B3B1B_ROOT
    / "freeze"
    / "32B3B1B_PROGRAM_INTERPRETATION.csv"
)

ANALYSIS_LOCK_PATH = (
    B3B1B_ROOT
    / "freeze"
    / "32B3B1B_ANALYSIS_LOCK.json"
)

GENERAL_PATHWAYS_PATH = (
    B3B1B_ROOT
    / "general_pathways"
    / "32B3B1B_GENERAL_PATHWAYS.csv"
)

EFFECT_PATH = (
    B3A_ROOT
    / "donor_effects"
    / "32B3A_PROGRAM_TLE_EFFECT_SUMMARY.csv"
)

TOP_GENE_PATH = (
    B3A_ROOT
    / "program_spectra"
    / "32B3A_TOP100_GENES.csv"
)

CONTENT_AUDIT_PATH = (
    B3B0_ROOT
    / "32B3B0_PROGRAM_CONTENT_AUDIT.csv"
)

for required_path in (
    INTERPRETATION_PATH,
    ANALYSIS_LOCK_PATH,
    GENERAL_PATHWAYS_PATH,
    EFFECT_PATH,
    TOP_GENE_PATH,
    CONTENT_AUDIT_PATH,
):
    if not required_path.is_file():
        raise FileNotFoundError(
            f"Required file missing: {required_path}"
        )


with ANALYSIS_LOCK_PATH.open(
    "r",
    encoding="utf-8",
) as handle:
    upstream_lock = json.load(handle)

if upstream_lock.get("status") != "PASS":
    raise RuntimeError(
        "32B3B1B is not in PASS status."
    )

if upstream_lock.get(
    "formal_multiaxis_candidates"
) != []:
    raise RuntimeError(
        "Unexpected formal multiaxis candidate was found."
    )

if upstream_lock.get(
    "exploratory_partial_or_mixed_candidates"
) != ["Oligodendrocytes__P5"]:
    raise RuntimeError(
        "Unexpected exploratory candidate universe."
    )


###############################################################################
# 3. Load inputs
###############################################################################

interpretation = pd.read_csv(
    INTERPRETATION_PATH
)

general_pathways = pd.read_csv(
    GENERAL_PATHWAYS_PATH
)

effects = pd.read_csv(
    EFFECT_PATH
)

top_genes = pd.read_csv(
    TOP_GENE_PATH
)

content_audit = pd.read_csv(
    CONTENT_AUDIT_PATH
)


EXPECTED_PROGRAMS = {
    "Astrocytes__P2",
    "Astrocytes__P4",
    "Astrocytes__P5",
    "Astrocytes__P6",
    "Astrocytes__P7",
    "Microglia_macrophage__P2",
    "Microglia_macrophage__P3",
    "Microglia_macrophage__P4",
    "OPC__P3",
    "Oligodendrocytes__P5",
    "Oligodendrocytes__P6",
}

if set(interpretation["program_id"]) != EXPECTED_PROGRAMS:
    raise RuntimeError(
        "The eligible program universe does not match the reviewed 11 programs."
    )


###############################################################################
# 4. Manually reviewed manuscript-safe labels
###############################################################################

LABEL_PLAN = {
    "Astrocytes__P2": {
        "safe_label":
            "Astrocyte-associated neural development and adhesion program",

        "evidence_tier":
            "LINEAGE_CONTEXT",

        "manuscript_role":
            "Reference state; not senescence-associated",

        "allowed_wording":
            (
                "An astrocyte-associated developmental and "
                "cell-adhesion transcriptional program."
            ),

        "prohibited_wording":
            "Astrocyte senescence or a TLE-associated astrocyte program.",

        "spatial_priority": 1,
        "experimental_priority": 1,
    },

    "Astrocytes__P4": {
        "safe_label":
            "Astrocyte chaperone-autophagy and proteostasis program",

        "evidence_tier":
            "CONTEXTUAL_SENESCENCE_FEATURES",

        "manuscript_role":
            "Senescence-associated features without disease association",

        "allowed_wording":
            (
                "A chaperone-autophagy/proteostasis program that "
                "overlapped several senescence-associated signatures "
                "but was not associated with TLE."
            ),

        "prohibited_wording":
            "TLE-associated astrocyte senescence.",

        "spatial_priority": 2,
        "experimental_priority": 1,
    },

    "Astrocytes__P5": {
        "safe_label":
            "Astrocyte inflammatory and immediate-response signaling program",

        "evidence_tier":
            "COMPETING_RESPONSE",

        "manuscript_role":
            "Competing inflammatory/stimulus response",

        "allowed_wording":
            (
                "An astrocytic inflammatory and stimulus-response "
                "program dominated by competing biological responses."
            ),

        "prohibited_wording":
            "Astrocyte senescence program.",

        "spatial_priority": 1,
        "experimental_priority": 1,
    },

    "Astrocytes__P6": {
        "safe_label":
            "Astrocyte transport, projection and barrier-associated program",

        "evidence_tier":
            "CONTEXTUAL_SENESCENCE_FEATURES",

        "manuscript_role":
            "Contextual state without disease association",

        "allowed_wording":
            (
                "An astrocyte transport and cell-projection program "
                "with senescence-associated features but no TLE effect."
            ),

        "prohibited_wording":
            "A TLE-driven astrocyte senescence state.",

        "spatial_priority": 2,
        "experimental_priority": 1,
    },

    "Astrocytes__P7": {
        "safe_label":
            "Astrocyte stimulus-response signaling program",

        "evidence_tier":
            "COMPETING_RESPONSE",

        "manuscript_role":
            "General stimulus-response context",

        "allowed_wording":
            "A general astrocytic stimulus-response signaling program.",

        "prohibited_wording":
            "Astrocyte senescence.",

        "spatial_priority": 1,
        "experimental_priority": 1,
    },

    "Microglia_macrophage__P2": {
        "safe_label":
            "Myeloid Rho-GTPase, motility and phagocytic signaling program",

        "evidence_tier":
            "COMPETING_RESPONSE",

        "manuscript_role":
            "Motility/phagocytic response context",

        "allowed_wording":
            (
                "A myeloid Rho-GTPase, migration and phagocytic "
                "signaling program."
            ),

        "prohibited_wording":
            "Senescent microglia.",

        "spatial_priority": 1,
        "experimental_priority": 1,
    },

    "Microglia_macrophage__P3": {
        "safe_label":
            "Myeloid lipid, oxygen-stress and phagolysosomal program",

        "evidence_tier":
            "CONTEXTUAL_SENESCENCE_FEATURES",

        "manuscript_role":
            "Mixed senescence-featured program without disease association",

        "allowed_wording":
            (
                "A myeloid lipid/oxygen-stress and phagolysosomal "
                "program with senescence-associated features but no "
                "TLE association."
            ),

        "prohibited_wording":
            "TLE-associated senescent microglia.",

        "spatial_priority": 2,
        "experimental_priority": 1,
    },

    "Microglia_macrophage__P4": {
        "safe_label":
            "Myeloid immune-receptor, stress and proliferative-response program",

        "evidence_tier":
            "COMPETING_RESPONSE",

        "manuscript_role":
            "Immune/stress competing response",

        "allowed_wording":
            (
                "A myeloid immune-receptor and stress-response program "
                "with a proliferative competing explanation."
            ),

        "prohibited_wording":
            "Stable cell-cycle exit or microglial senescence.",

        "spatial_priority": 1,
        "experimental_priority": 1,
    },

    "OPC__P3": {
        "safe_label":
            "OPC glutamatergic, adhesion and junction-response program",

        "evidence_tier":
            "COMPETING_RESPONSE",

        "manuscript_role":
            "Glutamatergic/adhesion response context",

        "allowed_wording":
            (
                "An OPC-associated glutamatergic, cell-adhesion and "
                "junction-response program."
            ),

        "prohibited_wording":
            "OPC senescence program.",

        "spatial_priority": 1,
        "experimental_priority": 1,
    },

    "Oligodendrocytes__P5": {
        "safe_label": (
            "TLE-associated oligodendroglial myelin–lysosomal/"
            "chaperone-autophagy–lipid stress program with partial "
            "senescence-associated transcriptional features"
        ),

        "evidence_tier":
            "HIGH_PRIORITY_EXPLORATORY",

        "manuscript_role":
            "Primary exploratory disease-linked program",

        "allowed_wording": (
            "A directionally supported TLE-associated oligodendroglial "
            "myelin–lysosomal/chaperone-autophagy–lipid stress program "
            "with partial senescence-associated transcriptional features."
        ),

        "prohibited_wording": (
            "Senescent oligodendrocytes drive TLE; an oligodendrocyte "
            "senescence program; causal oligodendrocyte senescence."
        ),

        "spatial_priority": 3,
        "experimental_priority": 3,
    },

    "Oligodendrocytes__P6": {
        "safe_label":
            "Oligodendrocyte differentiation and myelination program",

        "evidence_tier":
            "LINEAGE_CONTEXT",

        "manuscript_role":
            "Oligodendroglial identity/reference program",

        "allowed_wording":
            "An oligodendrocyte differentiation and myelination program.",

        "prohibited_wording":
            "Oligodendrocyte senescence program.",

        "spatial_priority": 2,
        "experimental_priority": 1,
    },
}


###############################################################################
# 5. Extract pathway anchors and top genes
###############################################################################

ANCHOR_TERMS = {
    "Oligodendrocytes__P5": [
        "oligodendrocyte differentiation",
        "central nervous system myelination",
        "Chaperone Mediated Autophagy",
        "lysosomal transport",
        "response to oxidative stress",
        "lipid transport",
        "iron ion transport",
    ],

    "Astrocytes__P4": [
        "Chaperone Mediated Autophagy",
        "Autophagy",
        "Aggrephagy",
    ],

    "Microglia_macrophage__P3": [
        "response to lipid",
        "response to oxygen-containing compound",
        "apoptotic process",
    ],
}


label_rows = []

for program_id in sorted(EXPECTED_PROGRAMS):
    label = LABEL_PLAN[program_id]

    interpretation_row = interpretation.loc[
        interpretation["program_id"] == program_id
    ].iloc[0]

    effect_row = effects.loc[
        effects["program_id"] == program_id
    ].iloc[0]

    pathway_subset = general_pathways.loc[
        general_pathways["program_id"] == program_id
    ].sort_values("adjusted_p_value")

    top_pathway_names = (
        pathway_subset["term_name"]
        .drop_duplicates()
        .head(8)
        .tolist()
    )

    requested_anchors = ANCHOR_TERMS.get(
        program_id,
        [],
    )

    detected_anchors = []

    for anchor in requested_anchors:
        match = pathway_subset.loc[
            pathway_subset["term_name"].str.lower()
            == anchor.lower()
        ]

        if not match.empty:
            detected_anchors.append(anchor)

    gene_subset = (
        top_genes.loc[
            top_genes["program_id"] == program_id
        ]
        .sort_values("rank")
    )

    top20 = (
        gene_subset["gene"]
        .astype(str)
        .head(20)
        .tolist()
    )

    label_rows.append(
        {
            "program_id": program_id,
            "lineage": interpretation_row["lineage"],
            "safe_label": label["safe_label"],
            "evidence_tier": label["evidence_tier"],
            "manuscript_role": label["manuscript_role"],
            "upstream_interpretation":
                interpretation_row["interpretation"],

            "primary_effect":
                effect_row["primary_effect"],

            "primary_hedges_g":
                effect_row["primary_hedges_g"],

            "primary_exact_p":
                effect_row["primary_exact_p"],

            "primary_fdr":
                effect_row["primary_fdr"],

            "sensitivity_effect":
                effect_row["sensitivity_effect"],

            "same_direction":
                effect_row["same_direction"],

            "top_general_pathways":
                ";".join(top_pathway_names),

            "reviewed_anchor_pathways":
                ";".join(detected_anchors),

            "top20_genes":
                ";".join(top20),

            "allowed_wording":
                label["allowed_wording"],

            "prohibited_wording":
                label["prohibited_wording"],

            "spatial_priority":
                label["spatial_priority"],

            "experimental_priority":
                label["experimental_priority"],
        }
    )


labels = pd.DataFrame(
    label_rows
)

labels.to_csv(
    OUTPUT_ROOT
    / "tables"
    / "32B3C_MANUSCRIPT_SAFE_PROGRAM_LABELS.csv",
    index=False,
)


###############################################################################
# 6. Claims and evidence boundaries
###############################################################################

claim_matrix = pd.DataFrame(
    [
        {
            "claim_id": "C1",
            "claim": (
                "No unbiased cNMF program met the predefined formal "
                "multiaxis senescence-associated threshold."
            ),
            "status": "SUPPORTED_NEGATIVE",
            "allowed_in_results": True,
            "causal": False,
        },

        {
            "claim_id": "C2",
            "claim": (
                "Oligodendrocytes__P5 was the only program with stable "
                "directional TLE support and partial senescence-associated "
                "transcriptional features."
            ),
            "status": "SUPPORTED_EXPLORATORY",
            "allowed_in_results": True,
            "causal": False,
        },

        {
            "claim_id": "C3",
            "claim": (
                "Oligodendrocytes__P5 reflects a mixed myelin, lysosomal/"
                "autophagy, lipid and stress response rather than a pure "
                "senescence program."
            ),
            "status": "SUPPORTED",
            "allowed_in_results": True,
            "causal": False,
        },

        {
            "claim_id": "C4",
            "claim": (
                "Astrocytic and myeloid programs showed senescence-associated "
                "features without a reproducible TLE effect."
            ),
            "status": "SUPPORTED",
            "allowed_in_results": True,
            "causal": False,
        },

        {
            "claim_id": "C5",
            "claim": (
                "Senescent oligodendrocytes drive or cause chronic TLE."
            ),
            "status": "PROHIBITED_CAUSAL_OVERCLAIM",
            "allowed_in_results": False,
            "causal": True,
        },
    ]
)

claim_matrix.to_csv(
    OUTPUT_ROOT
    / "tables"
    / "32B3C_CLAIM_MATRIX.csv",
    index=False,
)


###############################################################################
# 7. Validation priority
###############################################################################

validation_priority = pd.DataFrame(
    [
        {
            "priority": 1,
            "program_id": "Oligodendrocytes__P5",
            "validation_layer": "Spatial projection",
            "recommended_readout": (
                "Project the frozen meta-spectrum without refitting; "
                "test animal/donor-level regional effects and myelinated-tract "
                "localization."
            ),
        },

        {
            "priority": 2,
            "program_id": "Oligodendrocytes__P5",
            "validation_layer": "Lineage-localized tissue validation",
            "recommended_readout": (
                "OLIG2/SOX10 or CC1 with CRYAB, LAMP1/2, SCD, APOD, "
                "PLP1/CNP; quantify per animal, not per cell."
            ),
        },

        {
            "priority": 3,
            "program_id": "Oligodendrocytes__P5",
            "validation_layer": "Orthogonal senescence boundary",
            "recommended_readout": (
                "Test CDKN1A/p21, LMNB1, gamma-H2AX and lipofuscin/"
                "SA-beta-gal together with apoptosis and proliferation "
                "exclusion. A single marker is insufficient."
            ),
        },

        {
            "priority": 4,
            "program_id": (
                "Astrocytes__P4;Astrocytes__P6;"
                "Microglia_macrophage__P3"
            ),
            "validation_layer": "Context programs",
            "recommended_readout": (
                "Use only as contextual spatial modules; do not describe "
                "them as disease-associated without independent evidence."
            ),
        },
    ]
)

validation_priority.to_csv(
    OUTPUT_ROOT
    / "tables"
    / "32B3C_VALIDATION_PRIORITY.csv",
    index=False,
)


###############################################################################
# 8. Manuscript text
###############################################################################

results_text = """\
# Manuscript-safe Results text

Unbiased cNMF analysis did not identify a program that fulfilled the
predefined criteria for a formal multiaxis senescence-associated state.
Among the content-QC–eligible programs, Oligodendrocytes__P5 was the only
program showing consistent directional enrichment in TLE across two
independent balanced cell samples. The program was enriched for
oligodendrocyte differentiation and myelination, chaperone-mediated
autophagy and lysosomal transport, oxidative-stress responses, and lipid
and iron transport. It also overlapped several senescence-associated
transcriptional signatures, but concurrent enrichment of phagolysosomal/
lipid-processing and immediate-response programs indicated a mixed
biological state. Accordingly, we describe this signal as a
TLE-associated oligodendroglial myelin–lysosomal/chaperone-autophagy–lipid
stress program with partial senescence-associated transcriptional
features, rather than as evidence of cellular senescence.
"""

discussion_text = """\
# Manuscript-safe Discussion text

The results suggest that chronic TLE may involve transcriptional
remodelling within the oligodendroglial compartment rather than a stable
expansion of a discrete senescent cell population. The leading program
combined myelin and oligodendrocyte-identity functions with lysosomal,
chaperone-autophagy, lipid-handling and oxidative-stress features. Its
overlap with several senescence-associated signatures is hypothesis
generating, but is not sufficient to establish irreversible cellular
senescence. Competing injury and phagolysosomal/lipid-processing
responses remain plausible explanations. Spatial localization and
orthogonal multi-axis tissue validation will therefore be required
before assigning a senescence-related cellular phenotype.
"""

(
    OUTPUT_ROOT
    / "manuscript_text"
    / "32B3C_RESULTS_TEXT.md"
).write_text(
    results_text,
    encoding="utf-8",
)

(
    OUTPUT_ROOT
    / "manuscript_text"
    / "32B3C_DISCUSSION_TEXT.md"
).write_text(
    discussion_text,
    encoding="utf-8",
)


###############################################################################
# 9. Evidence-priority figure
###############################################################################

plot_table = labels.sort_values(
    [
        "experimental_priority",
        "primary_effect",
    ]
).reset_index(drop=True)

positions = np.arange(
    len(plot_table)
)

figure, axis = plt.subplots(
    figsize=(10, 7)
)

axis.axvline(
    0,
    linestyle="--",
    linewidth=1,
)

axis.scatter(
    plot_table["primary_effect"],
    positions,
    label="Primary balanced sample",
)

axis.scatter(
    plot_table["sensitivity_effect"],
    positions,
    marker="x",
    label="Independent balanced sample",
)

axis.set_yticks(positions)

axis.set_yticklabels(
    [
        f"{program}\n{tier}"
        for program, tier in zip(
            plot_table["program_id"],
            plot_table["evidence_tier"],
        )
    ],
    fontsize=8,
)

axis.set_xlabel(
    "Donor-level CLR program effect: TLE − Normal"
)

axis.set_title(
    "Frozen unbiased programs and manuscript-safe evidence tiers"
)

axis.legend()
figure.tight_layout()

figure.savefig(
    OUTPUT_ROOT
    / "figures"
    / "32B3C_PROGRAM_EVIDENCE_TIERS.png",
    dpi=350,
)

figure.savefig(
    OUTPUT_ROOT
    / "figures"
    / "32B3C_PROGRAM_EVIDENCE_TIERS.pdf",
)

plt.close(figure)


###############################################################################
# 10. Machine-readable lock
###############################################################################

lock = {
    "stage": "32B3C",
    "created_at": datetime.now().isoformat(),
    "status": "PASS",
    "formal_multiaxis_senescence_programs": [],
    "primary_exploratory_program":
        "Oligodendrocytes__P5",
    "primary_exploratory_safe_label":
        LABEL_PLAN[
            "Oligodendrocytes__P5"
        ]["safe_label"],
    "causal_claims_authorized": False,
    "senescent_cell_label_authorized": False,
    "spatial_projection_authorized": [
        "Oligodendrocytes__P5",
        "Astrocytes__P4",
        "Astrocytes__P6",
        "Microglia_macrophage__P3",
    ],
    "next_stage":
        "33A_SPATIAL_AND_CROSS_SPECIES_PROGRAM_VALIDATION",
}

with (
    OUTPUT_ROOT
    / "freeze"
    / "32B3C_MANUSCRIPT_SAFE_LOCK.json"
).open(
    "w",
    encoding="utf-8",
) as handle:
    json.dump(
        lock,
        handle,
        indent=2,
        ensure_ascii=False,
    )


###############################################################################
# 11. Executive summary
###############################################################################

summary = f"""\
# 32B-3C Manuscript-safe Program Freeze

## Status: PASS

## Primary exploratory program

`Oligodendrocytes__P5`

Safe label:

> {LABEL_PLAN["Oligodendrocytes__P5"]["safe_label"]}

## Evidence boundary

- No formal multiaxis senescence-associated program was identified.
- The primary candidate has directional but not FDR-significant TLE support.
- The program contains both senescence-associated and competing
  injury/lipid-processing features.
- Cellular senescence and causal language remain unauthorized.

## Next stage

Spatially and cross-species validate the frozen program without refitting
or changing its gene membership. Use animal or donor as the inferential unit.
"""

(
    OUTPUT_ROOT
    / "32B3C_EXECUTIVE_SUMMARY.md"
).write_text(
    summary,
    encoding="utf-8",
)


###############################################################################
# 12. Manifest and ZIP
###############################################################################

def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


output_files = [
    path
    for path in OUTPUT_ROOT.rglob("*")
    if path.is_file()
]

manifest_rows = []

for path in output_files:
    manifest_rows.append(
        {
            "relative_path":
                path.relative_to(
                    OUTPUT_ROOT
                ).as_posix(),
            "size_bytes":
                path.stat().st_size,
            "sha256":
                sha256_file(path),
        }
    )

pd.DataFrame(
    manifest_rows
).to_csv(
    OUTPUT_ROOT
    / "32B3C_FILE_MANIFEST.csv",
    index=False,
)

zip_path = Path(
    str(OUTPUT_ROOT)
    + ".zip"
)

with zipfile.ZipFile(
    zip_path,
    mode="w",
    compression=zipfile.ZIP_DEFLATED,
    compresslevel=9,
) as archive:

    for path in sorted(
        OUTPUT_ROOT.rglob("*")
    ):
        if path.is_file():
            archive.write(
                path,
                arcname=path.relative_to(
                    OUTPUT_ROOT
                ),
            )

zip_hash = sha256_file(
    zip_path
)

Path(
    str(zip_path)
    + ".sha256.txt"
).write_text(
    f"{zip_hash}  {zip_path.name}\n",
    encoding="ascii",
)


print()
print("=" * 72)
print("32B-3C completed successfully")
print("=" * 72)
print("Status: PASS")
print("Primary exploratory program: Oligodendrocytes__P5")
print()
print("Output:")
print(OUTPUT_ROOT)
print()
print("ZIP:")
print(zip_path)
print()
print("SHA256:")
print(zip_hash)
print("=" * 72)