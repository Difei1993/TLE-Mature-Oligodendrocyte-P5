###############################################################################
# 32B-3B0
# General biological-content and lineage-compatibility audit of frozen cNMF
# programs.
#
# No senescence gene sets are used.
# No new disease-effect test is performed.
###############################################################################

from __future__ import annotations

import hashlib
import json
import re
import zipfile
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import hypergeom


###############################################################################
# 0. Paths
###############################################################################

PROJECT = Path(
    __import__("os").environ["P5_WORK_ROOT"]
)

RELEASE_ROOT = PROJECT / "99_frozen_releases"

INPUT_CANDIDATES = [
    path
    for path in RELEASE_ROOT.glob(
        "32B3A_UNBIASED_PROGRAM_FREEZE_*"
    )
    if path.is_dir()
]

if not INPUT_CANDIDATES:
    raise FileNotFoundError(
        "No 32B3A_UNBIASED_PROGRAM_FREEZE directory was found."
    )

INPUT_ROOT = max(
    INPUT_CANDIDATES,
    key=lambda path: path.stat().st_mtime,
)

STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")

OUTPUT_ROOT = (
    RELEASE_ROOT
    / f"32B3B0_PROGRAM_CONTENT_AUDIT_{STAMP}"
)

OUTPUT_ROOT.mkdir(
    parents=True,
    exist_ok=False,
)

FIGURE_ROOT = OUTPUT_ROOT / "figures"
FIGURE_ROOT.mkdir()


###############################################################################
# 1. Input files
###############################################################################

LOCK_PATH = (
    INPUT_ROOT
    / "freeze"
    / "32B3A_FINAL_PROGRAM_LOCK.csv"
)

EFFECT_PATH = (
    INPUT_ROOT
    / "donor_effects"
    / "32B3A_PROGRAM_TLE_EFFECT_SUMMARY.csv"
)

TOP_GENE_PATH = (
    INPUT_ROOT
    / "program_spectra"
    / "32B3A_TOP100_GENES.csv"
)

for required_path in (
    LOCK_PATH,
    EFFECT_PATH,
    TOP_GENE_PATH,
):
    if not required_path.is_file():
        raise FileNotFoundError(
            f"Required file not found: {required_path}"
        )

program_lock = pd.read_csv(LOCK_PATH)
effects = pd.read_csv(EFFECT_PATH)
top_genes = pd.read_csv(TOP_GENE_PATH)


###############################################################################
# 2. Curated general marker modules
#
# These are used only for content/QC classification—not senescence testing.
###############################################################################

MARKER_SETS = {
    "NEURONAL_SYNAPTIC": set("""
        RBFOX1 RBFOX3 SYT1 SNAP25 NRXN1 NRXN2 NRXN3 CNTNAP2
        PTPRD DLGAP1 DLGAP2 DLGAP3 NLGN1 NLGN2 NRG3 ROBO2
        FGF14 KCNIP4 KCNQ5 KCND2 KCND3 RIMS1 RIMS2 GRM3 GRM5
        ERBB4 GABRB1 IL1RAPL1 LRRTM4 DPP10 HOMER1 NPTX2
        PLCB1 CADM2 ANKS1B FAM155A SCN11A DNM1 DNM3 ENC1
        CAMK1D
    """.split()),

    "ASTRO_HOMEOSTASIS": set("""
        AQP4 ALDH1L1 SLC1A2 SLC1A3 GLUL GJA1 ATP1A2 ATP1B2
        SOX9 NFIA NFIB APOE CLU MAOB CFI GRM3 ADGRV1 GPM6A
        DTNA PTPRZ1
    """.split()),

    "ASTRO_REACTIVE_ECM": set("""
        GFAP VIM CD44 OSMR STAT3 STAT1 SERPINA3 SERPINE1
        LCN2 C3 C4B ANXA2 TNC IL1R1 EMP1 WWTR1 NAMPT
        ADAMTS9 CLU APOE PLSCR1 COLEC12
    """.split()),

    "MICRO_HOMEOSTATIC": set("""
        P2RY12 P2RY13 CX3CR1 TMEM119 SALL1 MEF2C MEF2A
        GPR34 HEXB CSF1R TREM2 TYROBP DOCK8 LRMDA FRMD4A
        SORL1 BLNK APBB1IP
    """.split()),

    "MICRO_ACTIVATED": set("""
        SPP1 APOE LPL CST7 TLR2 CD83 RGS1 CCL4 CCL4L2
        IL1B TNF SERPINE1 FKBP5 CEBPD FOS JUN NR4A3 IRAK3
        TGFBI ADGRE5 BACH1 RUNX1
    """.split()),

    "OPC_PROGENITOR": set("""
        PDGFRA CSPG4 VCAN OLIG1 OLIG2 SOX10 PTPRZ1 BCAN
        NKX2-2 DLL3 GPR17 TNR
    """.split()),

    "OLIGO_MYELIN": set("""
        PLP1 MBP MOBP MAG CNP MOG CLDN11 ERMN FA2H QKI
        SOX10 OLIG1 OLIG2 MAL TF SELENOP PTGDS APOD SCD
        PPP1R14A SLC44A1
    """.split()),

    "VASCULAR": set("""
        PECAM1 VWF EMCN KDR FLT1 CLDN5 RAMP2 RGCC APOLD1
        CA4 PLVAP RBP7 COL4A1 COL4A2 PDGFRB CSPG4 RGS5
        ACTA2 TAGLN
    """.split()),

    "CELL_CYCLE": set("""
        MKI67 TOP2A CENPF CENPK RRM2 NUSAP1 TACC3 KIF15
        ASPM CDC25C NCAPG BRCA2 BRIP1 DIAPH3 KNL1 SPC24
        SMC4 POC1A UBE2C BIRC5 CCNB1 CCNB2 CDK1
    """.split()),

    "IMMEDIATE_EARLY_GC_STRESS": set("""
        FOS JUN JUNB EGR1 EGR2 DUSP1 NR4A1 NR4A2 NR4A3
        FKBP5 SGK1 CEBPD PDK4 RGS1
    """.split()),

    "INTERFERON": set("""
        STAT1 STAT2 IRF1 IRF7 IFIT1 IFIT2 IFIT3 ISG15
        MX1 OAS1 OAS2 GBP1 GBP2 PLSCR1 B2M HLA-A HLA-B
    """.split()),

    "CHAPERONE_OXIDATIVE_STRESS": set("""
        HSP90AA1 HSPA1A HSPA1B HSPB1 CRYAB FTH1 FTL MT3
        UBB UBC CLU NDUFA4 GPX1 SOD1 SOD2 HMOX1
    """.split()),

    "LIPID_MYELIN_METABOLISM": set("""
        SCD APOE APOD PLP1 CNP SELENOP PTGDS FASN HMGCR
        SREBF1 LPL ABCA1 ABCG1 SOAT1 FABP5 GALC UGT8 FA2H
    """.split()),
}


EXPECTED_MODULES = {
    "Astrocytes": {
        "ASTRO_HOMEOSTASIS",
        "ASTRO_REACTIVE_ECM",
    },

    "Microglia_macrophage": {
        "MICRO_HOMEOSTATIC",
        "MICRO_ACTIVATED",
    },

    "OPC": {
        "OPC_PROGENITOR",
        "OLIGO_MYELIN",
    },

    "Oligodendrocytes": {
        "OLIGO_MYELIN",
        "LIPID_MYELIN_METABOLISM",
    },
}


INCOMPATIBLE_MODULES = {
    "Astrocytes": {
        "NEURONAL_SYNAPTIC",
        "MICRO_HOMEOSTATIC",
        "MICRO_ACTIVATED",
        "OLIGO_MYELIN",
        "VASCULAR",
    },

    "Microglia_macrophage": {
        "NEURONAL_SYNAPTIC",
        "ASTRO_HOMEOSTASIS",
        "OLIGO_MYELIN",
    },

    "OPC": {
        "NEURONAL_SYNAPTIC",
        "MICRO_HOMEOSTATIC",
        "MICRO_ACTIVATED",
        "VASCULAR",
    },

    "Oligodendrocytes": {
        "NEURONAL_SYNAPTIC",
        "ASTRO_HOMEOSTASIS",
        "MICRO_HOMEOSTATIC",
        "MICRO_ACTIVATED",
        "VASCULAR",
    },
}


###############################################################################
# 3. Utility functions
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


def benjamini_hochberg(values: pd.Series) -> pd.Series:
    values_array = values.to_numpy(dtype=float)
    number = len(values_array)

    order = np.argsort(values_array)
    ordered = values_array[order]

    adjusted_ordered = (
        ordered
        * number
        / np.arange(1, number + 1)
    )

    adjusted_ordered = np.minimum.accumulate(
        adjusted_ordered[::-1]
    )[::-1]

    adjusted_ordered = np.minimum(
        adjusted_ordered,
        1.0,
    )

    adjusted = np.empty(number)
    adjusted[order] = adjusted_ordered

    return pd.Series(
        adjusted,
        index=values.index,
    )


def is_lnc_or_unannotated(gene: str) -> bool:
    return bool(
        re.match(
            r"^(AC|AL|AP|LINC|MIR|SNHG)",
            gene,
        )
        or re.search(
            r"-AS\d*$",
            gene,
        )
    )


###############################################################################
# 4. Marker-overlap audit
###############################################################################

overlap_rows = []
audit_rows = []

for lock_row in program_lock.itertuples(
    index=False
):
    program_id = lock_row.program_id
    lineage = lock_row.lineage

    ranked_table = (
        top_genes.loc[
            top_genes["program_id"] == program_id
        ]
        .sort_values("rank")
    )

    if len(ranked_table) != 100:
        raise RuntimeError(
            f"{program_id}: expected 100 genes, found "
            f"{len(ranked_table)}"
        )

    genes_20 = set(
        ranked_table.loc[
            ranked_table["rank"] <= 20,
            "gene",
        ].astype(str)
    )

    genes_50 = set(
        ranked_table.loc[
            ranked_table["rank"] <= 50,
            "gene",
        ].astype(str)
    )

    genes_100 = set(
        ranked_table["gene"].astype(str)
    )

    spectrum_path = (
        INPUT_ROOT
        / "program_spectra"
        / f"{program_id}__meta_spectrum.csv"
    )

    spectrum = pd.read_csv(
        spectrum_path
    )

    universe = set(
        spectrum["gene"].astype(str)
    )

    marker_counts_50 = {}

    for marker_name, marker_genes in MARKER_SETS.items():
        marker_in_universe = (
            marker_genes & universe
        )

        overlap_20 = (
            genes_20 & marker_genes
        )

        overlap_50 = (
            genes_50 & marker_genes
        )

        overlap_100 = (
            genes_100 & marker_genes
        )

        population_size = len(universe)
        success_population = len(
            marker_in_universe
        )
        sample_size = len(genes_100)
        observed = len(overlap_100)

        if (
            population_size > 0
            and success_population > 0
        ):
            p_value = float(
                hypergeom.sf(
                    observed - 1,
                    population_size,
                    success_population,
                    sample_size,
                )
            )
        else:
            p_value = 1.0

        marker_counts_50[
            marker_name
        ] = len(overlap_50)

        overlap_rows.append(
            {
                "lineage": lineage,
                "program_id": program_id,
                "marker_module": marker_name,
                "overlap_top20": len(overlap_20),
                "overlap_top50": len(overlap_50),
                "overlap_top100": len(overlap_100),
                "overlapping_genes_top100": ";".join(
                    sorted(overlap_100)
                ),
                "universe_genes": population_size,
                "marker_genes_in_universe":
                    success_population,
                "hypergeometric_p": p_value,
            }
        )

    dominant_module = max(
        marker_counts_50,
        key=marker_counts_50.get,
    )

    dominant_count = marker_counts_50[
        dominant_module
    ]

    expected_count = max(
        marker_counts_50[module]
        for module in EXPECTED_MODULES[lineage]
    )

    incompatible_count = max(
        marker_counts_50[module]
        for module in INCOMPATIBLE_MODULES[
            lineage
        ]
    )

    lnc_fraction_top20 = (
        sum(
            is_lnc_or_unannotated(gene)
            for gene in genes_20
        )
        / len(genes_20)
    )

    cell_cycle_count = marker_counts_50[
        "CELL_CYCLE"
    ]

    stress_count = marker_counts_50[
        "CHAPERONE_OXIDATIVE_STRESS"
    ]

    if cell_cycle_count >= 5:
        content_flag = (
            "ACTIVE_CELL_CYCLE_COMPETING_STATE"
        )

    elif (
        incompatible_count >= 5
        and incompatible_count
        > expected_count + 1
    ):
        content_flag = (
            "LIKELY_LINEAGE_INCOMPATIBLE_RNA_PROGRAM"
        )

    elif lnc_fraction_top20 >= 0.50:
        content_flag = (
            "LOW_ANNOTATABILITY_LNCRNA_DOMINANT"
        )

    elif stress_count >= 4:
        content_flag = (
            "STRESS_CHAPERONE_PROGRAM"
        )

    else:
        content_flag = (
            "BIOLOGICALLY_PLAUSIBLE_PENDING_ENRICHMENT"
        )

    eligible_for_senescence_overlap = bool(
        lock_row.evidence_status
        == "CORE_REPLICATED"
        and content_flag not in {
            "ACTIVE_CELL_CYCLE_COMPETING_STATE",
            "LIKELY_LINEAGE_INCOMPATIBLE_RNA_PROGRAM",
            "LOW_ANNOTATABILITY_LNCRNA_DOMINANT",
        }
    )

    eligible_for_competing_analysis = bool(
        lock_row.evidence_status
        != "EXCLUDED_UNSTABLE"
    )

    audit_rows.append(
        {
            "lineage": lineage,
            "program_id": program_id,
            "original_evidence_status":
                lock_row.evidence_status,

            "cross_sampling_correlation":
                lock_row.cross_sampling_correlation,

            "dominant_general_module":
                dominant_module,

            "dominant_overlap_top50":
                dominant_count,

            "expected_lineage_overlap_top50":
                expected_count,

            "incompatible_overlap_top50":
                incompatible_count,

            "lnc_fraction_top20":
                lnc_fraction_top20,

            "cell_cycle_overlap_top50":
                cell_cycle_count,

            "stress_overlap_top50":
                stress_count,

            "content_flag":
                content_flag,

            "eligible_for_senescence_overlap":
                eligible_for_senescence_overlap,

            "eligible_for_competing_analysis":
                eligible_for_competing_analysis,
        }
    )


overlap_table = pd.DataFrame(
    overlap_rows
)

overlap_table[
    "marker_overlap_fdr"
] = benjamini_hochberg(
    overlap_table[
        "hypergeometric_p"
    ]
)

audit_table = pd.DataFrame(
    audit_rows
)


###############################################################################
# 5. Merge disease-effect information
###############################################################################

effect_columns = [
    "program_id",
    "primary_effect",
    "primary_hedges_g",
    "primary_exact_p",
    "primary_fdr",
    "primary_lodo",
    "sensitivity_effect",
    "sensitivity_exact_p",
    "sensitivity_lodo",
    "same_direction",
]

audit_table = audit_table.merge(
    effects[effect_columns],
    on="program_id",
    how="left",
    validate="one_to_one",
)


def assign_review_priority(
    row: pd.Series,
) -> str:

    if row[
        "original_evidence_status"
    ] == "EXCLUDED_UNSTABLE":
        return "EXCLUDED_UNSTABLE"

    if row[
        "content_flag"
    ] == "LIKELY_LINEAGE_INCOMPATIBLE_RNA_PROGRAM":
        return "QC_EXCLUDE_LINEAGE_INCOMPATIBLE"

    if row[
        "content_flag"
    ] == "LOW_ANNOTATABILITY_LNCRNA_DOMINANT":
        return "MANUAL_REVIEW_LOW_ANNOTATABILITY"

    if row[
        "content_flag"
    ] == "ACTIVE_CELL_CYCLE_COMPETING_STATE":
        return "COMPETING_ACTIVE_PROLIFERATION"

    exploratory_disease_support = bool(
        row["eligible_for_senescence_overlap"]
        and row["same_direction"]
        and row["primary_exact_p"] < 0.10
        and row["primary_lodo"] >= 0.75
        and row["sensitivity_lodo"] >= 0.75
    )

    if exploratory_disease_support:
        return "HIGH_PRIORITY_EXPLORATORY_DISEASE_PROGRAM"

    if row[
        "eligible_for_senescence_overlap"
    ]:
        return "STANDARD_BIOLOGICAL_ANNOTATION"

    return "EXPLORATORY_ONLY"


audit_table[
    "review_priority"
] = audit_table.apply(
    assign_review_priority,
    axis=1,
)


###############################################################################
# 6. Save tables
###############################################################################

overlap_table.to_csv(
    OUTPUT_ROOT
    / "32B3B0_MARKER_OVERLAP_LONG.csv",
    index=False,
)

audit_table.to_csv(
    OUTPUT_ROOT
    / "32B3B0_PROGRAM_CONTENT_AUDIT.csv",
    index=False,
)

review_queue = audit_table.sort_values(
    [
        "review_priority",
        "primary_exact_p",
        "program_id",
    ]
)

review_queue.to_csv(
    OUTPUT_ROOT
    / "32B3B0_PROGRAM_REVIEW_QUEUE.csv",
    index=False,
)


###############################################################################
# 7. Marker-overlap heatmap
###############################################################################

heatmap_table = (
    overlap_table.pivot(
        index="program_id",
        columns="marker_module",
        values="overlap_top50",
    )
    .loc[
        audit_table.sort_values(
            [
                "lineage",
                "program_id",
            ]
        )["program_id"]
    ]
)

figure_height = max(
    8,
    0.32 * len(heatmap_table),
)

figure, axis = plt.subplots(
    figsize=(12, figure_height)
)

image = axis.imshow(
    heatmap_table.to_numpy(),
    aspect="auto",
)

axis.set_xticks(
    range(heatmap_table.shape[1])
)

axis.set_xticklabels(
    heatmap_table.columns,
    rotation=45,
    ha="right",
)

axis.set_yticks(
    range(heatmap_table.shape[0])
)

axis.set_yticklabels(
    heatmap_table.index
)

axis.set_title(
    "General marker overlap among frozen unbiased cNMF programs"
)

figure.colorbar(
    image,
    ax=axis,
    label="Marker genes among top 50 program genes",
)

figure.tight_layout()

figure.savefig(
    FIGURE_ROOT
    / "32B3B0_MARKER_OVERLAP_HEATMAP.png",
    dpi=300,
)

figure.savefig(
    FIGURE_ROOT
    / "32B3B0_MARKER_OVERLAP_HEATMAP.pdf",
)

plt.close(figure)


###############################################################################
# 8. Effect figure after content audit
###############################################################################

plot_table = audit_table.sort_values(
    [
        "lineage",
        "primary_effect",
    ]
).reset_index(drop=True)

positions = np.arange(
    len(plot_table)
)

figure, axis = plt.subplots(
    figsize=(
        10,
        max(7, 0.32 * len(plot_table)),
    )
)

axis.axvline(
    0,
    linestyle="--",
    linewidth=1,
)

axis.scatter(
    plot_table["primary_effect"],
    positions,
    label="Primary sampling",
)

axis.scatter(
    plot_table["sensitivity_effect"],
    positions,
    marker="x",
    label="Sensitivity sampling",
)

axis.set_yticks(positions)

axis.set_yticklabels(
    [
        (
            f"{program} | {flag}"
        )
        for program, flag in zip(
            plot_table["program_id"],
            plot_table["content_flag"],
        )
    ],
    fontsize=7,
)

axis.set_xlabel(
    "Donor CLR effect: TLE − Normal"
)

axis.set_title(
    "Disease effects interpreted after program-content QC"
)

axis.legend()
figure.tight_layout()

figure.savefig(
    FIGURE_ROOT
    / "32B3B0_EFFECTS_AFTER_CONTENT_QC.png",
    dpi=300,
)

figure.savefig(
    FIGURE_ROOT
    / "32B3B0_EFFECTS_AFTER_CONTENT_QC.pdf",
)

plt.close(figure)


###############################################################################
# 9. Executive summary
###############################################################################

incompatible_programs = audit_table.loc[
    audit_table[
        "content_flag"
    ] == "LIKELY_LINEAGE_INCOMPATIBLE_RNA_PROGRAM",
    "program_id",
].tolist()

cell_cycle_programs = audit_table.loc[
    audit_table[
        "content_flag"
    ] == "ACTIVE_CELL_CYCLE_COMPETING_STATE",
    "program_id",
].tolist()

eligible_programs = audit_table.loc[
    audit_table[
        "eligible_for_senescence_overlap"
    ],
    "program_id",
].tolist()

high_priority = audit_table.loc[
    audit_table[
        "review_priority"
    ] == "HIGH_PRIORITY_EXPLORATORY_DISEASE_PROGRAM",
    "program_id",
].tolist()

summary = f"""\
# 32B-3B0 Program Content and Lineage-Compatibility Audit

## Status: COMPLETED

Input:
`{INPUT_ROOT}`

## Core findings

- Frozen programs audited: {len(audit_table)}
- Likely lineage-incompatible RNA programs: {len(incompatible_programs)}
- Active cell-cycle competing programs: {len(cell_cycle_programs)}
- Programs eligible for later senescence-overlap testing: {len(eligible_programs)}

## Likely lineage-incompatible programs

{", ".join(incompatible_programs) if incompatible_programs else "None"}

## Active proliferation/cell-cycle programs

{", ".join(cell_cycle_programs) if cell_cycle_programs else "None"}

## High-priority exploratory disease-associated programs after content QC

{", ".join(high_priority) if high_priority else "None"}

## Interpretation boundary

- Reproducibility does not establish lineage compatibility.
- A lineage-incompatible program may reflect ambient RNA, doublets,
  anatomical admixture or donor-specific contamination.
- Active proliferation is a competing explanation and cannot be interpreted
  as stable cell-cycle exit.
- Stress/chaperone programs remain eligible for later analysis but must be
  evaluated against oxidative-stress and injury-response explanations.
- No senescence gene set was used in this stage.

## Next stage

Run pathway annotation on programs that pass content QC, then test overlap
with the previously frozen senescence programs and competing-response modules.
"""

(
    OUTPUT_ROOT
    / "32B3B0_EXECUTIVE_SUMMARY.md"
).write_text(
    summary,
    encoding="utf-8",
)


###############################################################################
# 10. Lock, manifest and ZIP
###############################################################################

lock = {
    "stage": "32B3B0",
    "created_at": datetime.now().isoformat(),
    "input_root": str(INPUT_ROOT),
    "senescence_gene_sets_used": False,
    "new_disease_tests_performed": False,
    "lineage_incompatible_programs":
        incompatible_programs,
    "active_cell_cycle_programs":
        cell_cycle_programs,
    "eligible_for_senescence_overlap":
        eligible_programs,
    "high_priority_exploratory_programs":
        high_priority,
    "next_stage":
        "32B3B1_PATHWAY_AND_FROZEN_SENESCENCE_INTERPRETATION",
}

with (
    OUTPUT_ROOT
    / "32B3B0_CONTENT_AUDIT_LOCK.json"
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

output_files = [
    path
    for path in OUTPUT_ROOT.rglob("*")
    if path.is_file()
]

manifest_rows = []

for path in output_files:
    manifest_rows.append(
        {
            "relative_path": str(
                path.relative_to(
                    OUTPUT_ROOT
                )
            ),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    )

pd.DataFrame(
    manifest_rows
).to_csv(
    OUTPUT_ROOT
    / "32B3B0_FILE_MANIFEST.csv",
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
print("32B-3B0 completed successfully")
print("=" * 72)
print(f"Input: {INPUT_ROOT}")
print(f"Programs audited: {len(audit_table)}")
print(
    "Lineage-incompatible programs:",
    incompatible_programs,
)
print(
    "High-priority exploratory programs:",
    high_priority,
)
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