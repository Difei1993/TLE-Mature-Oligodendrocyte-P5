###############################################################################
# 32B-3B1B
# General pathway annotation followed by exact frozen-26E interpretation.
#
# This stage:
#   - uses only programs that passed 32B3B0 content QC;
#   - performs no new donor-level disease test;
#   - keeps primary, external, deconfounded and competing layers separate;
#   - uses program-specific gene universes.
###############################################################################

from __future__ import annotations

import hashlib
import json
import math
import re
import time
import zipfile
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from scipy.stats import hypergeom, mannwhitneyu


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
    / f"32B3B1B_PATHWAY_AND_26E_INTERPRETATION_{STAMP}"
)

OUTPUT_ROOT.mkdir(
    parents=True,
    exist_ok=False,
)

(OUTPUT_ROOT / "general_pathways").mkdir()
(OUTPUT_ROOT / "frozen_26E").mkdir()
(OUTPUT_ROOT / "figures").mkdir()
(OUTPUT_ROOT / "audit").mkdir()
(OUTPUT_ROOT / "freeze").mkdir()


###############################################################################
# 1. Locate the latest approved upstream releases
###############################################################################

def latest_directory(prefix: str) -> Path:
    candidates = [
        path
        for path in RELEASE_ROOT.glob(f"{prefix}*")
        if path.is_dir()
    ]

    if not candidates:
        raise FileNotFoundError(
            f"No directory found for prefix: {prefix}"
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

CATALOG_ROOT = latest_directory(
    "32B3B1A_FROZEN_26E_CATALOG_AUDIT_"
)


###############################################################################
# 2. Required inputs and status gates
###############################################################################

CONTENT_PATH = (
    B3B0_ROOT
    / "32B3B0_PROGRAM_CONTENT_AUDIT.csv"
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

CATALOG_PATH = (
    CATALOG_ROOT
    / "freeze"
    / "32B3B1A_FROZEN_26E_GENESETS_LONG.csv"
)

CATALOG_LOCK_PATH = (
    CATALOG_ROOT
    / "32B3B1A_CATALOG_LOCK.json"
)

required_files = [
    CONTENT_PATH,
    EFFECT_PATH,
    TOP_GENE_PATH,
    CATALOG_PATH,
    CATALOG_LOCK_PATH,
]

for required_path in required_files:
    if not required_path.is_file():
        raise FileNotFoundError(
            f"Required file missing: {required_path}"
        )

with CATALOG_LOCK_PATH.open(
    "r",
    encoding="utf-8",
) as handle:
    catalog_lock = json.load(handle)

if (
    catalog_lock.get("status")
    != "READY_EXACT_26E_CATALOG_FROZEN"
):
    raise RuntimeError(
        "The exact 26E catalog is not in READY status."
    )


###############################################################################
# 3. Load data
###############################################################################

content_audit = pd.read_csv(
    CONTENT_PATH
)

effects = pd.read_csv(
    EFFECT_PATH
)

top_genes = pd.read_csv(
    TOP_GENE_PATH
)

catalog = pd.read_csv(
    CATALOG_PATH
)

eligible_programs = (
    content_audit.loc[
        content_audit[
            "eligible_for_senescence_overlap"
        ].astype(bool),
        "program_id",
    ]
    .astype(str)
    .tolist()
)

if len(eligible_programs) != 11:
    raise RuntimeError(
        "Expected 11 content-QC-eligible programs, "
        f"found {len(eligible_programs)}."
    )

if catalog["program"].nunique() != 26:
    raise RuntimeError(
        "Frozen catalog does not contain exactly 26 programs."
    )


###############################################################################
# 4. Utilities
###############################################################################

def now_string() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value

    return str(value).strip().lower() in {
        "true",
        "t",
        "1",
        "yes",
    }


def bh_adjust(values: pd.Series) -> pd.Series:
    numeric = values.to_numpy(dtype=float)

    number = len(numeric)
    order = np.argsort(numeric)
    ordered = numeric[order]

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


def safe_json_text(value: object) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
    )


def normalized_layer(value: str) -> str:
    value = str(value)

    if value == "Primary senescence dimension":
        return "PRIMARY_SENESCENCE_AXIS"

    if value == "External published signature":
        return "EXTERNAL_SIGNATURE"

    if value == "Deconfounded sensitivity dimension":
        return "DECONFOUNDED_SENSITIVITY"

    if value == "Inverse supplemental dimension":
        return "INVERSE_SUPPLEMENTAL"

    if value == "Competing biological program":
        return "COMPETING_PROGRAM"

    return "OTHER"


def spectrum_path(program_id: str) -> Path:
    return (
        B3A_ROOT
        / "program_spectra"
        / f"{program_id}__meta_spectrum.csv"
    )


###############################################################################
# 5. Record input provenance
###############################################################################

provenance = {
    "stage": "32B3B1B",
    "created_at": now_string(),
    "B3A_root": str(B3A_ROOT),
    "B3B0_root": str(B3B0_ROOT),
    "catalog_root": str(CATALOG_ROOT),
    "catalog_source":
        catalog_lock.get("selected_catalog"),
    "eligible_programs": eligible_programs,
    "new_disease_tests_performed": False,
    "program_discovery_reopened": False,
}

with (
    OUTPUT_ROOT
    / "audit"
    / "32B3B1B_INPUT_PROVENANCE.json"
).open(
    "w",
    encoding="utf-8",
) as handle:
    json.dump(
        provenance,
        handle,
        indent=2,
        ensure_ascii=False,
    )


###############################################################################
# 6. General pathway annotation through g:Profiler
###############################################################################

GPROFILER_PROFILE_URL = (
    "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
)

GPROFILER_VERSION_URL = (
    "https://biit.cs.ut.ee/gprofiler/api/util/data_versions/"
)

session = requests.Session()

session.headers.update(
    {
        "User-Agent":
            "TLE-senescence-32B3B1B/1.0",
        "Content-Type": "application/json",
        "Accept": "application/json",
    }
)

general_rows = []
general_errors = []

try:
    version_response = session.get(
        GPROFILER_VERSION_URL,
        params={"organism": "hsapiens"},
        timeout=120,
    )

    version_response.raise_for_status()

    version_data = version_response.json()

    with (
        OUTPUT_ROOT
        / "general_pathways"
        / "gprofiler_data_versions.json"
    ).open(
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            version_data,
            handle,
            indent=2,
            ensure_ascii=False,
        )

except Exception as error:
    version_data = {
        "error": repr(error)
    }

    general_errors.append(
        {
            "program_id": "DATA_VERSION",
            "error": repr(error),
        }
    )


for program_id in eligible_programs:
    ranked = (
        top_genes.loc[
            top_genes["program_id"] == program_id
        ]
        .sort_values("rank")
    )

    query_genes = (
        ranked.loc[
            ranked["rank"] <= 100,
            "gene",
        ]
        .astype(str)
        .str.upper()
        .drop_duplicates()
        .tolist()
    )

    spectrum = pd.read_csv(
        spectrum_path(program_id)
    )

    background_genes = (
        spectrum["gene"]
        .astype(str)
        .str.upper()
        .drop_duplicates()
        .tolist()
    )

    payload = {
        "organism": "hsapiens",
        "query": query_genes,
        "sources": [
            "GO:BP",
            "REAC",
        ],
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
        "domain_scope": "custom",
        "background": background_genes,
        "no_evidences": True,
    }

    try:
        response = session.post(
            GPROFILER_PROFILE_URL,
            json=payload,
            timeout=180,
        )

        response.raise_for_status()
        result = response.json()

        raw_output_path = (
            OUTPUT_ROOT
            / "general_pathways"
            / f"{program_id}__gprofiler_raw.json"
        )

        with raw_output_path.open(
            "w",
            encoding="utf-8",
        ) as handle:
            json.dump(
                result,
                handle,
                indent=2,
                ensure_ascii=False,
            )

        for item in result.get(
            "result",
            []
        ):
            general_rows.append(
                {
                    "program_id": program_id,
                    "source": item.get("source"),
                    "term_id": item.get("native"),
                    "term_name": item.get("name"),
                    "adjusted_p_value":
                        item.get("p_value"),
                    "term_size":
                        item.get("term_size"),
                    "query_size":
                        item.get("query_size"),
                    "intersection_size":
                        item.get(
                            "intersection_size"
                        ),
                    "effective_domain_size":
                        item.get(
                            "effective_domain_size"
                        ),
                    "precision":
                        item.get("precision"),
                    "recall":
                        item.get("recall"),
                    "intersections":
                        safe_json_text(
                            item.get(
                                "intersections",
                                [],
                            )
                        ),
                    "background_mode":
                        "PROGRAM_SPECIFIC_CUSTOM",
                }
            )

    except Exception as error:
        general_errors.append(
            {
                "program_id": program_id,
                "error": repr(error),
            }
        )

    time.sleep(0.5)


general_pathways = pd.DataFrame(
    general_rows
)

general_errors_table = pd.DataFrame(
    general_errors
)

general_pathways.to_csv(
    OUTPUT_ROOT
    / "general_pathways"
    / "32B3B1B_GENERAL_PATHWAYS.csv",
    index=False,
)

general_errors_table.to_csv(
    OUTPUT_ROOT
    / "general_pathways"
    / "32B3B1B_GENERAL_PATHWAY_ERRORS.csv",
    index=False,
)


###############################################################################
# 7. Exact 26E overlap and ranked-spectrum tests
###############################################################################

catalog["gene"] = (
    catalog["gene"]
    .astype(str)
    .str.upper()
)

gene_sets = {
    program: set(
        subset["gene"]
    )
    for program, subset
    in catalog.groupby("program")
}

catalog_metadata = (
    catalog[
        [
            "program",
            "analysis_class",
            "evidence_layer",
            "expected_direction",
            "direction_multiplier",
        ]
    ]
    .drop_duplicates()
    .set_index("program")
)

overlap_rows = []


for program_id in eligible_programs:
    ranked_top = (
        top_genes.loc[
            top_genes["program_id"] == program_id
        ]
        .sort_values("rank")
        .copy()
    )

    spectrum = pd.read_csv(
        spectrum_path(program_id)
    )

    spectrum["gene"] = (
        spectrum["gene"]
        .astype(str)
        .str.upper()
    )

    spectrum = (
        spectrum
        .drop_duplicates(
            subset=["gene"]
        )
        .set_index("gene")
    )

    universe = set(
        spectrum.index
    )

    all_scores = spectrum[
        "meta_program_score"
    ].astype(float)

    top100 = set(
        ranked_top.loc[
            ranked_top["rank"] <= 100,
            "gene",
        ]
        .astype(str)
        .str.upper()
    )

    top50 = set(
        ranked_top.loc[
            ranked_top["rank"] <= 50,
            "gene",
        ]
        .astype(str)
        .str.upper()
    )

    top20 = set(
        ranked_top.loc[
            ranked_top["rank"] <= 20,
            "gene",
        ]
        .astype(str)
        .str.upper()
    )

    for gene_set_name, gene_set in gene_sets.items():
        metadata_row = catalog_metadata.loc[
            gene_set_name
        ]

        genes_in_universe = (
            gene_set & universe
        )

        non_set_genes = (
            universe - genes_in_universe
        )

        if len(genes_in_universe) >= 3:
            set_values = all_scores.loc[
                sorted(genes_in_universe)
            ].to_numpy(dtype=float)

            non_set_values = all_scores.loc[
                sorted(non_set_genes)
            ].to_numpy(dtype=float)

            test = mannwhitneyu(
                set_values,
                non_set_values,
                alternative="greater",
                method="asymptotic",
            )

            rank_p = float(
                test.pvalue
            )

            rank_auc = float(
                test.statistic
                / (
                    len(set_values)
                    * len(non_set_values)
                )
            )

        else:
            rank_p = 1.0
            rank_auc = float("nan")

        overlap_100 = (
            top100 & gene_set
        )

        overlap_p = float(
            hypergeom.sf(
                len(overlap_100) - 1,
                len(universe),
                len(genes_in_universe),
                len(top100),
            )
        ) if genes_in_universe else 1.0

        overlap_rows.append(
            {
                "program_id": program_id,
                "gene_set": gene_set_name,
                "analysis_class":
                    metadata_row[
                        "analysis_class"
                    ],
                "evidence_layer":
                    metadata_row[
                        "evidence_layer"
                    ],
                "layer_group":
                    normalized_layer(
                        metadata_row[
                            "evidence_layer"
                        ]
                    ),
                "expected_direction":
                    metadata_row[
                        "expected_direction"
                    ],
                "direction_multiplier":
                    int(
                        metadata_row[
                            "direction_multiplier"
                        ]
                    ),
                "universe_size":
                    len(universe),
                "gene_set_size_catalog":
                    len(gene_set),
                "gene_set_size_in_universe":
                    len(genes_in_universe),
                "overlap_top20":
                    len(
                        top20 & gene_set
                    ),
                "overlap_top50":
                    len(
                        top50 & gene_set
                    ),
                "overlap_top100":
                    len(overlap_100),
                "overlapping_genes_top100":
                    ";".join(
                        sorted(overlap_100)
                    ),
                "hypergeometric_p":
                    overlap_p,
                "rank_auc":
                    rank_auc,
                "rank_mannwhitney_p":
                    rank_p,
            }
        )


overlap_table = pd.DataFrame(
    overlap_rows
)

overlap_table[
    "hypergeometric_fdr_global"
] = bh_adjust(
    overlap_table[
        "hypergeometric_p"
    ]
)

overlap_table[
    "rank_fdr_global"
] = bh_adjust(
    overlap_table[
        "rank_mannwhitney_p"
    ]
)

overlap_table[
    "hypergeometric_fdr_within_program"
] = (
    overlap_table
    .groupby(
        "program_id",
        group_keys=False,
    )["hypergeometric_p"]
    .apply(bh_adjust)
)

overlap_table[
    "rank_fdr_within_program"
] = (
    overlap_table
    .groupby(
        "program_id",
        group_keys=False,
    )["rank_mannwhitney_p"]
    .apply(bh_adjust)
)

overlap_table[
    "formal_gene_set_hit"
] = (
    (
        (
            overlap_table[
                "hypergeometric_fdr_global"
            ] < 0.10
        )
        & (
            overlap_table[
                "overlap_top100"
            ] >= 2
        )
    )
    |
    (
        (
            overlap_table[
                "rank_fdr_global"
            ] < 0.10
        )
        & (
            overlap_table[
                "rank_auc"
            ] >= 0.60
        )
    )
)

overlap_table[
    "exploratory_gene_set_hit"
] = (
    ~overlap_table[
        "formal_gene_set_hit"
    ]
    & (
        (
            (
                overlap_table[
                    "hypergeometric_p"
                ] < 0.05
            )
            & (
                overlap_table[
                    "overlap_top100"
                ] >= 2
            )
        )
        |
        (
            (
                overlap_table[
                    "rank_mannwhitney_p"
                ] < 0.05
            )
            & (
                overlap_table[
                    "rank_auc"
                ] >= 0.58
            )
        )
    )
)

overlap_table.to_csv(
    OUTPUT_ROOT
    / "frozen_26E"
    / "32B3B1B_26E_OVERLAP_LONG.csv",
    index=False,
)


###############################################################################
# 8. Merge disease directions without repeating disease tests
###############################################################################

effect_columns = [
    "program_id",
    "evidence_status",
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

effect_table = effects[
    effect_columns
].copy()

effect_table["same_direction"] = (
    effect_table[
        "same_direction"
    ].apply(as_bool)
)

overlap_table = overlap_table.merge(
    effect_table,
    on="program_id",
    how="left",
    validate="many_to_one",
)

overlap_table[
    "aligned_primary_effect"
] = (
    overlap_table[
        "primary_effect"
    ]
    * overlap_table[
        "direction_multiplier"
    ]
)

overlap_table[
    "aligned_sensitivity_effect"
] = (
    overlap_table[
        "sensitivity_effect"
    ]
    * overlap_table[
        "direction_multiplier"
    ]
)

overlap_table[
    "disease_direction_aligned"
] = (
    overlap_table[
        "aligned_primary_effect"
    ] > 0
) & (
    overlap_table[
        "aligned_sensitivity_effect"
    ] > 0
)

overlap_table.to_csv(
    OUTPUT_ROOT
    / "frozen_26E"
    / "32B3B1B_26E_OVERLAP_WITH_DISEASE_DIRECTION.csv",
    index=False,
)


###############################################################################
# 9. Program-level interpretation
###############################################################################

summary_rows = []

for program_id in eligible_programs:
    subset = overlap_table.loc[
        overlap_table[
            "program_id"
        ] == program_id
    ].copy()

    effect_row = effect_table.loc[
        effect_table[
            "program_id"
        ] == program_id
    ].iloc[0]

    formal_hits = subset.loc[
        subset[
            "formal_gene_set_hit"
        ]
    ]

    exploratory_hits = subset.loc[
        subset[
            "exploratory_gene_set_hit"
        ]
    ]

    aligned_formal_hits = formal_hits.loc[
        formal_hits[
            "disease_direction_aligned"
        ]
    ]

    primary_hits = aligned_formal_hits.loc[
        aligned_formal_hits[
            "layer_group"
        ] == "PRIMARY_SENESCENCE_AXIS"
    ]

    external_hits = aligned_formal_hits.loc[
        aligned_formal_hits[
            "layer_group"
        ] == "EXTERNAL_SIGNATURE"
    ]

    inverse_hits = aligned_formal_hits.loc[
        aligned_formal_hits[
            "layer_group"
        ] == "INVERSE_SUPPLEMENTAL"
    ]

    competing_hits = formal_hits.loc[
        formal_hits[
            "layer_group"
        ] == "COMPETING_PROGRAM"
    ]

    deconfounded_hits = aligned_formal_hits.loc[
        aligned_formal_hits[
            "layer_group"
        ] == "DECONFOUNDED_SENSITIVITY"
    ]

    disease_directional_support = bool(
        effect_row[
            "same_direction"
        ]
        and effect_row[
            "primary_exact_p"
        ] < 0.10
        and effect_row[
            "primary_lodo"
        ] >= 0.75
        and effect_row[
            "sensitivity_lodo"
        ] >= 0.75
    )

    disease_formal_support = bool(
        disease_directional_support
        and pd.notna(
            effect_row[
                "primary_fdr"
            ]
        )
        and effect_row[
            "primary_fdr"
        ] < 0.10
    )

    strong_competing = bool(
        competing_hits[
            "gene_set"
        ].isin(
            [
                "Proliferation",
                "Immediate_early_response",
                "IFN_response",
                "NFkB_inflammation_chemokine",
                "Phagolysosome_lipid_processing",
            ]
        ).any()
    )

    independent_senescence_axes = (
        len(primary_hits)
        + len(inverse_hits)
    )

    external_support = len(
        external_hits
    )

    if (
        disease_formal_support
        and independent_senescence_axes >= 2
        and external_support >= 1
        and not strong_competing
    ):
        interpretation = (
            "FORMAL_MULTIAXIS_SENESCENCE_ASSOCIATED_PROGRAM"
        )

    elif (
        disease_directional_support
        and (
            independent_senescence_axes >= 2
            or (
                independent_senescence_axes >= 1
                and external_support >= 1
            )
        )
        and not strong_competing
    ):
        interpretation = (
            "EXPLORATORY_PARTIAL_SENESCENCE_CONSISTENT"
        )

    elif (
        disease_directional_support
        and (
            independent_senescence_axes >= 1
            or external_support >= 1
        )
        and strong_competing
    ):
        interpretation = (
            "EXPLORATORY_MIXED_SENESCENCE_AND_COMPETING_RESPONSE"
        )

    elif (
        independent_senescence_axes >= 1
        or external_support >= 1
        or len(deconfounded_hits) >= 1
    ):
        interpretation = (
            "SENESCENCE_FEATURES_WITHOUT_DISEASE_ASSOCIATION"
        )

    elif len(competing_hits) >= 1:
        interpretation = (
            "COMPETING_RESPONSE_DOMINANT"
        )

    else:
        interpretation = (
            "INSUFFICIENT_26E_ALIGNMENT"
        )

    summary_rows.append(
        {
            "program_id": program_id,
            "disease_formal_support":
                disease_formal_support,
            "disease_directional_support":
                disease_directional_support,
            "primary_senescence_axes_aligned":
                len(primary_hits),
            "external_signatures_aligned":
                len(external_hits),
            "inverse_dimensions_aligned":
                len(inverse_hits),
            "deconfounded_hits_aligned":
                len(deconfounded_hits),
            "competing_program_hits":
                len(competing_hits),
            "strong_competing_explanation":
                strong_competing,
            "formal_hit_names":
                ";".join(
                    formal_hits[
                        "gene_set"
                    ].tolist()
                ),
            "exploratory_hit_names":
                ";".join(
                    exploratory_hits[
                        "gene_set"
                    ].tolist()
                ),
            "interpretation":
                interpretation,
        }
    )


program_summary = pd.DataFrame(
    summary_rows
)

program_summary = program_summary.merge(
    content_audit[
        [
            "program_id",
            "lineage",
            "content_flag",
            "review_priority",
        ]
    ],
    on="program_id",
    how="left",
    validate="one_to_one",
)

program_summary.to_csv(
    OUTPUT_ROOT
    / "freeze"
    / "32B3B1B_PROGRAM_INTERPRETATION.csv",
    index=False,
)


###############################################################################
# 10. Heatmap
###############################################################################

heatmap_input = overlap_table.copy()

heatmap_input[
    "display_value"
] = np.where(
    heatmap_input[
        "formal_gene_set_hit"
    ],
    -np.log10(
        np.maximum(
            heatmap_input[
                "hypergeometric_fdr_global"
            ],
            1e-300,
        )
    ),
    0,
)

heatmap = (
    heatmap_input.pivot(
        index="program_id",
        columns="gene_set",
        values="display_value",
    )
    .loc[eligible_programs]
)

figure, axis = plt.subplots(
    figsize=(16, 7)
)

image = axis.imshow(
    heatmap.to_numpy(),
    aspect="auto",
)

axis.set_xticks(
    range(heatmap.shape[1])
)

axis.set_xticklabels(
    heatmap.columns,
    rotation=70,
    ha="right",
    fontsize=7,
)

axis.set_yticks(
    range(heatmap.shape[0])
)

axis.set_yticklabels(
    heatmap.index,
    fontsize=8,
)

axis.set_title(
    "Frozen 26E overlap after unbiased cNMF discovery"
)

figure.colorbar(
    image,
    ax=axis,
    label="-log10 global FDR; zero denotes no formal hit",
)

figure.tight_layout()

figure.savefig(
    OUTPUT_ROOT
    / "figures"
    / "32B3B1B_26E_OVERLAP_HEATMAP.png",
    dpi=350,
)

figure.savefig(
    OUTPUT_ROOT
    / "figures"
    / "32B3B1B_26E_OVERLAP_HEATMAP.pdf",
)

plt.close(figure)


###############################################################################
# 11. Completion status and report
###############################################################################

general_pathway_success = (
    len(general_errors) == 0
)

STATUS = (
    "PASS"
    if general_pathway_success
    else "CONDITIONAL_PASS_GENERAL_PATHWAY_API_INCOMPLETE"
)

mixed_candidates = program_summary.loc[
    program_summary[
        "interpretation"
    ].str.contains(
        "MIXED|PARTIAL",
        regex=True,
    ),
    "program_id",
].tolist()

formal_candidates = program_summary.loc[
    program_summary[
        "interpretation"
    ]
    == "FORMAL_MULTIAXIS_SENESCENCE_ASSOCIATED_PROGRAM",
    "program_id",
].tolist()

report = f"""\
# 32B-3B1B General Pathway and Frozen-26E Interpretation

## Status: {STATUS}

## Inputs

- B3A program freeze: `{B3A_ROOT}`
- B3B0 content audit: `{B3B0_ROOT}`
- Exact 26E catalog: `{CATALOG_ROOT}`

## Design safeguards

- cNMF programs were frozen before this analysis.
- Only 11 programs that passed content QC were tested.
- No new donor-level disease test was performed.
- Primary, external, deconfounded and competing layers were kept separate.
- Deconfounded gene sets were not counted as independent senescence axes.
- Program-specific gene universes were used.
- General GO/Reactome annotation preceded final senescence interpretation.

## Formal multiaxis senescence-associated programs

{", ".join(formal_candidates) if formal_candidates else "None."}

## Exploratory partial or mixed programs

{", ".join(mixed_candidates) if mixed_candidates else "None."}

## Interpretation boundary

A single SASP gene, external signature, lysosomal signal or oxidative-stress
program is insufficient to define cellular senescence. Competing programs
such as immediate-early response, inflammation, proliferation and
phagolysosomal/lipid processing remain explicit alternative explanations.

## Next stage

Freeze manuscript-safe labels, identify the highest-priority non-causal
programs, and decide which programs merit spatial mapping or experimental
validation.
"""

(
    OUTPUT_ROOT
    / "32B3B1B_EXECUTIVE_SUMMARY.md"
).write_text(
    report,
    encoding="utf-8",
)


###############################################################################
# 12. Lock, manifest and ZIP
###############################################################################

lock = {
    "stage": "32B3B1B",
    "created_at": now_string(),
    "status": STATUS,
    "eligible_program_count":
        len(eligible_programs),
    "general_pathway_api_complete":
        general_pathway_success,
    "formal_multiaxis_candidates":
        formal_candidates,
    "exploratory_partial_or_mixed_candidates":
        mixed_candidates,
    "new_disease_tests_performed":
        False,
    "deconfounded_sets_counted_as_independent_axes":
        False,
    "next_stage":
        "32B3C_MANUSCRIPT_SAFE_PROGRAM_FREEZE",
}

with (
    OUTPUT_ROOT
    / "freeze"
    / "32B3B1B_ANALYSIS_LOCK.json"
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
            "relative_path":
                str(
                    path.relative_to(
                        OUTPUT_ROOT
                    )
                ),
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
    / "32B3B1B_FILE_MANIFEST.csv",
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
print("32B-3B1B completed")
print("=" * 72)
print(f"Status: {STATUS}")
print(f"Programs tested: {len(eligible_programs)}")
print("Formal candidates:", formal_candidates)
print("Exploratory partial/mixed candidates:", mixed_candidates)
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