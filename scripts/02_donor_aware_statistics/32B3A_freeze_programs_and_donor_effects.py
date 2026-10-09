###############################################################################
# 32B-3A
# Freeze reproducible unbiased cNMF programs and estimate donor-level TLE effects
#
# Discovery is already complete and locked.
# Disease labels are used only now, after program discovery.
# Senescence gene sets are NOT used in this stage.
###############################################################################

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import shutil
import time
import zipfile
from datetime import datetime
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


###############################################################################
# 0. Paths
###############################################################################

PROJECT = Path(
    __import__("os").environ["P5_WORK_ROOT"]
)

CNMF_ROOT = (
    PROJECT
    / "04_tables"
    / "32B2B_FINAL_CNMf_20260728_215801"
)

PRIMARY_H5AD_ROOT = (
    PROJECT
    / "04_tables"
    / "32B_CNMf_WORK_20260728_215801"
    / "_h5ad_inputs"
)

SENSITIVITY_H5AD_ROOT = (
    PROJECT
    / "04_tables"
    / "32B2A_CNMf_SENSITIVITY_20260728_215801"
    / "_h5ad_inputs"
)

OUTPUT_ROOT = (
    PROJECT
    / "99_frozen_releases"
    / (
        "32B3A_UNBIASED_PROGRAM_FREEZE_"
        + datetime.now().strftime("%Y%m%d_%H%M%S")
    )
)

OUTPUT_ROOT.mkdir(
    parents=True,
    exist_ok=False,
)

(OUTPUT_ROOT / "program_spectra").mkdir()
(OUTPUT_ROOT / "donor_effects").mkdir()
(OUTPUT_ROOT / "figures").mkdir()
(OUTPUT_ROOT / "audit").mkdir()
(OUTPUT_ROOT / "freeze").mkdir()


###############################################################################
# 1. Final frozen program plan
###############################################################################

PLAN = {
    "Astrocytes": {
        "k": 7,
        "threshold": 0.05,
        "primary_to_sensitivity": {
            1: 1,
            2: 2,
            3: 3,
            4: 4,
            5: 5,
            6: 6,
            7: 7,
        },
        "core_primary_programs": [1, 2, 3, 4, 5, 6, 7],
        "exploratory_primary_programs": [],
        "excluded_primary_programs": [],
    },

    "Microglia_macrophage": {
        "k": 6,
        "threshold": 0.01,
        "primary_to_sensitivity": {
            1: 1,
            2: 2,
            3: 3,
            4: 4,
            5: 5,
            6: 6,
        },
        "core_primary_programs": [1, 2, 3, 4, 5],
        "exploratory_primary_programs": [6],
        "excluded_primary_programs": [],
    },

    "OPC": {
        "k": 6,
        "threshold": 0.05,
        "primary_to_sensitivity": {
            1: 4,
            2: 2,
            3: 5,
            4: 1,
            5: 3,
            6: 6,
        },
        "core_primary_programs": [2, 3, 6],
        "exploratory_primary_programs": [1],
        "excluded_primary_programs": [4, 5],
    },

    "Oligodendrocytes": {
        "k": 6,
        "threshold": 0.01,
        "primary_to_sensitivity": {
            1: 1,
            2: 2,
            3: 3,
            4: 4,
            5: 5,
            6: 6,
        },
        "core_primary_programs": [1, 2, 3, 4, 5, 6],
        "exploratory_primary_programs": [],
        "excluded_primary_programs": [],
    },
}


H5AD_FILES = {
    "primary": {
        "Astrocytes":
            PRIMARY_H5AD_ROOT / "Astrocytes_primary.h5ad",

        "Microglia_macrophage":
            PRIMARY_H5AD_ROOT / "Microglia_macrophage_primary.h5ad",

        "OPC":
            PRIMARY_H5AD_ROOT / "OPC_primary.h5ad",

        "Oligodendrocytes":
            PRIMARY_H5AD_ROOT / "Oligodendrocytes_primary.h5ad",
    },

    "sensitivity": {
        "Astrocytes":
            SENSITIVITY_H5AD_ROOT / "Astrocytes_sensitivity.h5ad",

        "Microglia_macrophage":
            SENSITIVITY_H5AD_ROOT
            / "Microglia_macrophage_sensitivity.h5ad",

        "OPC":
            SENSITIVITY_H5AD_ROOT / "OPC_sensitivity.h5ad",

        "Oligodendrocytes":
            SENSITIVITY_H5AD_ROOT
            / "Oligodendrocytes_sensitivity.h5ad",
    },
}


###############################################################################
# 2. Utility functions
###############################################################################

def now_string() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def threshold_tag(value: float) -> str:
    return str(value).replace(".", "_")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def locate_unique(
    root: Path,
    pattern: str,
) -> Path:

    matches = [
        path
        for path in root.rglob(pattern)
        if path.is_file()
        and "cnmf_tmp" not in path.parts
    ]

    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one file matching {pattern} "
            f"in {root}; found {len(matches)}:\n"
            + "\n".join(str(path) for path in matches)
        )

    return matches[0]


def read_program_matrix(
    path: Path,
    k: int,
) -> pd.DataFrame:

    table = pd.read_csv(
        path,
        sep="\t",
        index_col=0,
    )

    if table.shape[0] == k:
        result = table.copy()

    elif table.shape[1] == k:
        result = table.T.copy()

    else:
        raise ValueError(
            f"Cannot orient program matrix {path}; "
            f"shape={table.shape}, expected one dimension={k}"
        )

    result.index = [
        f"Program_{index + 1}"
        for index in range(k)
    ]

    result.columns = result.columns.astype(str)

    return result


def read_usage_matrix(
    path: Path,
    k: int,
) -> pd.DataFrame:

    table = pd.read_csv(
        path,
        sep="\t",
        index_col=0,
    )

    if table.shape[1] == k:
        result = table.copy()

    elif table.shape[0] == k:
        result = table.T.copy()

    else:
        raise ValueError(
            f"Cannot orient usage matrix {path}; "
            f"shape={table.shape}, expected one dimension={k}"
        )

    result.columns = [
        f"Program_{index + 1}"
        for index in range(k)
    ]

    result.index = result.index.astype(str)

    return result


def normalize_usage_rows(
    usage: pd.DataFrame,
) -> pd.DataFrame:

    values = usage.to_numpy(dtype=float)

    if np.any(values < 0):
        raise ValueError("Negative cNMF usages detected.")

    row_sums = values.sum(axis=1)

    if np.any(row_sums <= 0):
        raise ValueError("One or more cells have zero total usage.")

    normalized = values / row_sums[:, None]

    return pd.DataFrame(
        normalized,
        index=usage.index,
        columns=usage.columns,
    )


def clr_transform(
    donor_usage: pd.DataFrame,
    pseudocount: float = 1e-6,
) -> pd.DataFrame:

    values = donor_usage.to_numpy(dtype=float)

    log_values = np.log(
        values + pseudocount
    )

    clr_values = (
        log_values
        - log_values.mean(
            axis=1,
            keepdims=True,
        )
    )

    return pd.DataFrame(
        clr_values,
        index=donor_usage.index,
        columns=donor_usage.columns,
    )


def exact_permutation_test(
    values: np.ndarray,
    groups: np.ndarray,
) -> tuple[float, float]:

    groups = np.asarray(groups)
    values = np.asarray(values, dtype=float)

    tle_mask = groups == "TLE"
    normal_mask = groups == "Normal"

    if tle_mask.sum() != 4 or normal_mask.sum() != 4:
        raise ValueError(
            "Exact test requires 4 TLE and 4 Normal donors."
        )

    observed = (
        values[tle_mask].mean()
        - values[normal_mask].mean()
    )

    effects = []

    for tle_indices in itertools.combinations(
        range(len(values)),
        4,
    ):
        permutation_mask = np.zeros(
            len(values),
            dtype=bool,
        )

        permutation_mask[
            list(tle_indices)
        ] = True

        effect = (
            values[permutation_mask].mean()
            - values[~permutation_mask].mean()
        )

        effects.append(effect)

    effects = np.asarray(effects)

    p_value = np.mean(
        np.abs(effects)
        >= abs(observed) - 1e-12
    )

    return float(observed), float(p_value)


def hedges_g(
    values: np.ndarray,
    groups: np.ndarray,
) -> float:

    tle = values[groups == "TLE"]
    normal = values[groups == "Normal"]

    n1 = len(tle)
    n0 = len(normal)

    pooled_variance = (
        (n1 - 1) * np.var(tle, ddof=1)
        + (n0 - 1) * np.var(normal, ddof=1)
    ) / (
        n1 + n0 - 2
    )

    if (
        not np.isfinite(pooled_variance)
        or pooled_variance <= 0
    ):
        return float("nan")

    d_value = (
        np.mean(tle)
        - np.mean(normal)
    ) / np.sqrt(pooled_variance)

    correction = (
        1
        - 3 / (
            4 * (n1 + n0) - 9
        )
    )

    return float(correction * d_value)


def lodo_direction_stability(
    values: np.ndarray,
    groups: np.ndarray,
) -> float:

    full_effect, _ = exact_permutation_test(
        values,
        groups,
    )

    full_sign = np.sign(full_effect)

    if full_sign == 0:
        return 0.0

    consistent = []

    for omitted in range(len(values)):
        keep = np.arange(len(values)) != omitted

        reduced_values = values[keep]
        reduced_groups = groups[keep]

        tle_mean = reduced_values[
            reduced_groups == "TLE"
        ].mean()

        normal_mean = reduced_values[
            reduced_groups == "Normal"
        ].mean()

        reduced_effect = (
            tle_mean - normal_mean
        )

        consistent.append(
            np.sign(reduced_effect)
            == full_sign
        )

    return float(
        np.mean(consistent)
    )


def benjamini_hochberg(
    p_values: pd.Series,
) -> pd.Series:

    values = p_values.to_numpy(dtype=float)
    number = len(values)

    order = np.argsort(values)
    ranked = values[order]

    adjusted_ranked = (
        ranked
        * number
        / np.arange(1, number + 1)
    )

    adjusted_ranked = np.minimum.accumulate(
        adjusted_ranked[::-1]
    )[::-1]

    adjusted_ranked = np.minimum(
        adjusted_ranked,
        1.0,
    )

    adjusted = np.empty(number)
    adjusted[order] = adjusted_ranked

    return pd.Series(
        adjusted,
        index=p_values.index,
    )


def standardize_gene_vector(
    values: pd.Series,
) -> pd.Series:

    values = values.astype(float)

    standard_deviation = values.std(ddof=0)

    if (
        not np.isfinite(standard_deviation)
        or standard_deviation <= 0
    ):
        return values * 0

    return (
        values - values.mean()
    ) / standard_deviation


###############################################################################
# 3. Validate inputs
###############################################################################

if not CNMF_ROOT.is_dir():
    raise FileNotFoundError(
        f"cNMF output root not found: {CNMF_ROOT}"
    )

for sampling, lineage_files in H5AD_FILES.items():
    for lineage, path in lineage_files.items():

        if not path.is_file():
            raise FileNotFoundError(
                f"Missing {sampling} H5AD for "
                f"{lineage}: {path}"
            )

print("All cNMF and H5AD inputs passed validation.")


###############################################################################
# 4. Freeze program matrices and donor-level usages
###############################################################################

program_lock_rows = []
top_gene_rows = []
donor_effect_rows = []
donor_usage_rows = []

overall_started = time.time()

for lineage, parameters in PLAN.items():

    k_value = parameters["k"]
    threshold = parameters["threshold"]
    tag = threshold_tag(threshold)

    spectra_by_sampling = {}
    usage_by_sampling = {}
    metadata_by_sampling = {}

    for sampling in (
        "primary",
        "sensitivity",
    ):

        run_directory = (
            CNMF_ROOT
            / lineage
            / sampling
            / (
                f"{lineage}_"
                f"{sampling}_"
                f"final_k{k_value}"
            )
        )

        spectra_path = locate_unique(
            run_directory,
            (
                f"*.gene_spectra_score."
                f"k_{k_value}.dt_{tag}*.txt"
            ),
        )

        usage_path = locate_unique(
            run_directory,
            (
                f"*.usages."
                f"k_{k_value}.dt_{tag}."
                f"consensus.txt"
            ),
        )

        spectra = read_program_matrix(
            spectra_path,
            k_value,
        )

        usage = normalize_usage_rows(
            read_usage_matrix(
                usage_path,
                k_value,
            )
        )

        h5ad_path = H5AD_FILES[
            sampling
        ][lineage]

        adata = ad.read_h5ad(
            h5ad_path,
            backed="r",
        )

        metadata = adata.obs[
            [
                "donor_id",
                "condition",
            ]
        ].copy()

        metadata.index = metadata.index.astype(str)
        metadata["donor_id"] = (
            metadata["donor_id"]
            .astype(str)
        )

        metadata["condition"] = (
            metadata["condition"]
            .astype(str)
        )

        adata.file.close()

        common_cells = (
            usage.index.intersection(
                metadata.index
            )
        )

        if len(common_cells) != len(usage):
            raise RuntimeError(
                f"{lineage}/{sampling}: only "
                f"{len(common_cells)}/{len(usage)} usage cells "
                "matched H5AD metadata."
            )

        usage = usage.loc[common_cells]
        metadata = metadata.loc[common_cells]

        spectra_by_sampling[sampling] = spectra
        usage_by_sampling[sampling] = usage
        metadata_by_sampling[sampling] = metadata

    ###########################################################################
    # Donor-level usage and CLR transformation
    ###########################################################################

    donor_clr_by_sampling = {}
    donor_condition_by_sampling = {}

    for sampling in (
        "primary",
        "sensitivity",
    ):

        usage = usage_by_sampling[sampling]
        metadata = metadata_by_sampling[sampling]

        combined = usage.copy()
        combined["donor_id"] = metadata[
            "donor_id"
        ]

        donor_usage = combined.groupby(
            "donor_id"
        )[usage.columns].mean()

        condition_map = (
            metadata[
                [
                    "donor_id",
                    "condition",
                ]
            ]
            .drop_duplicates()
            .set_index("donor_id")
        )

        if condition_map.index.duplicated().any():
            raise RuntimeError(
                f"{lineage}/{sampling}: donor-condition conflict."
            )

        condition_map = condition_map.loc[
            donor_usage.index
        ]

        donor_clr = clr_transform(
            donor_usage
        )

        donor_clr_by_sampling[sampling] = donor_clr
        donor_condition_by_sampling[sampling] = (
            condition_map["condition"]
        )

        donor_output = donor_clr.copy()
        donor_output.insert(
            0,
            "condition",
            condition_map["condition"],
        )
        donor_output.insert(
            0,
            "donor_id",
            donor_output.index,
        )

        donor_output.to_csv(
            OUTPUT_ROOT
            / "donor_effects"
            / (
                f"{lineage}__"
                f"{sampling}__"
                "donor_CLR_usage.csv"
            ),
            index=False,
        )

    ###########################################################################
    # Freeze programs and average matched spectra
    ###########################################################################

    mapping = parameters[
        "primary_to_sensitivity"
    ]

    core_programs = set(
        parameters[
            "core_primary_programs"
        ]
    )

    exploratory_programs = set(
        parameters[
            "exploratory_primary_programs"
        ]
    )

    excluded_programs = set(
        parameters[
            "excluded_primary_programs"
        ]
    )

    for primary_number in range(
        1,
        k_value + 1,
    ):

        sensitivity_number = mapping[
            primary_number
        ]

        if primary_number in core_programs:
            evidence_status = (
                "CORE_REPLICATED"
            )

        elif primary_number in exploratory_programs:
            evidence_status = (
                "EXPLORATORY"
            )

        elif primary_number in excluded_programs:
            evidence_status = (
                "EXCLUDED_UNSTABLE"
            )

        else:
            raise RuntimeError(
                f"{lineage} Program {primary_number} "
                "has no evidence status."
            )

        program_id = (
            f"{lineage}__P{primary_number}"
        )

        primary_vector = spectra_by_sampling[
            "primary"
        ].loc[
            f"Program_{primary_number}"
        ]

        sensitivity_vector = spectra_by_sampling[
            "sensitivity"
        ].loc[
            f"Program_{sensitivity_number}"
        ]

        common_genes = sorted(
            set(primary_vector.index)
            & set(sensitivity_vector.index)
        )

        primary_vector = primary_vector.loc[
            common_genes
        ]

        sensitivity_vector = (
            sensitivity_vector.loc[
                common_genes
            ]
        )

        correlation = float(
            np.corrcoef(
                primary_vector,
                sensitivity_vector,
            )[0, 1]
        )

        meta_spectrum = (
            standardize_gene_vector(
                primary_vector
            )
            + standardize_gene_vector(
                sensitivity_vector
            )
        ) / 2

        meta_spectrum = meta_spectrum.sort_values(
            ascending=False
        )

        spectrum_table = pd.DataFrame(
            {
                "gene": meta_spectrum.index,
                "meta_program_score":
                    meta_spectrum.values,

                "primary_score":
                    primary_vector.loc[
                        meta_spectrum.index
                    ].values,

                "sensitivity_score":
                    sensitivity_vector.loc[
                        meta_spectrum.index
                    ].values,
            }
        )

        spectrum_table.to_csv(
            OUTPUT_ROOT
            / "program_spectra"
            / f"{program_id}__meta_spectrum.csv",
            index=False,
        )

        for rank, row in enumerate(
            spectrum_table.head(100).itertuples(
                index=False
            ),
            start=1,
        ):
            top_gene_rows.append(
                {
                    "lineage": lineage,
                    "program_id": program_id,
                    "evidence_status":
                        evidence_status,
                    "rank": rank,
                    "gene": row.gene,
                    "meta_program_score":
                        row.meta_program_score,
                }
            )

        program_lock_rows.append(
            {
                "lineage": lineage,
                "k": k_value,
                "threshold": threshold,
                "program_id": program_id,
                "primary_program":
                    primary_number,
                "sensitivity_program":
                    sensitivity_number,
                "cross_sampling_correlation":
                    correlation,
                "evidence_status":
                    evidence_status,
                "allowed_for_formal_testing":
                    evidence_status
                    == "CORE_REPLICATED",
            }
        )

        #######################################################################
        # Donor-level disease effect in both cell samples
        #######################################################################

        for sampling in (
            "primary",
            "sensitivity",
        ):

            program_number = (
                primary_number
                if sampling == "primary"
                else sensitivity_number
            )

            program_column = (
                f"Program_{program_number}"
            )

            donor_clr = donor_clr_by_sampling[
                sampling
            ]

            conditions = (
                donor_condition_by_sampling[
                    sampling
                ]
                .loc[donor_clr.index]
                .to_numpy()
            )

            values = donor_clr[
                program_column
            ].to_numpy(dtype=float)

            effect, p_value = (
                exact_permutation_test(
                    values,
                    conditions,
                )
            )

            effect_size = hedges_g(
                values,
                conditions,
            )

            lodo_stability = (
                lodo_direction_stability(
                    values,
                    conditions,
                )
            )

            donor_effect_rows.append(
                {
                    "lineage": lineage,
                    "program_id": program_id,
                    "evidence_status":
                        evidence_status,
                    "sampling": sampling,
                    "tested_program_number":
                        program_number,
                    "effect_TLE_minus_Normal":
                        effect,
                    "hedges_g": effect_size,
                    "exact_p_value": p_value,
                    "lodo_direction_stability":
                        lodo_stability,
                    "direction": (
                        "TLE_higher"
                        if effect > 0
                        else "TLE_lower"
                        if effect < 0
                        else "no_difference"
                    ),
                }
            )

            for donor_id, condition, value in zip(
                donor_clr.index,
                conditions,
                values,
            ):
                donor_usage_rows.append(
                    {
                        "lineage": lineage,
                        "program_id": program_id,
                        "sampling": sampling,
                        "donor_id": donor_id,
                        "condition": condition,
                        "clr_usage": value,
                    }
                )


###############################################################################
# 5. Statistical tables
###############################################################################

program_lock = pd.DataFrame(
    program_lock_rows
)

top_genes = pd.DataFrame(
    top_gene_rows
)

donor_effects = pd.DataFrame(
    donor_effect_rows
)

donor_usage_long = pd.DataFrame(
    donor_usage_rows
)

formal_primary_mask = (
    (donor_effects["sampling"] == "primary")
    & (
        donor_effects["evidence_status"]
        == "CORE_REPLICATED"
    )
)

donor_effects[
    "formal_primary_fdr"
] = np.nan

donor_effects.loc[
    formal_primary_mask,
    "formal_primary_fdr",
] = benjamini_hochberg(
    donor_effects.loc[
        formal_primary_mask,
        "exact_p_value",
    ]
).to_numpy()

primary_effects = (
    donor_effects[
        donor_effects["sampling"]
        == "primary"
    ]
    .set_index(
        [
            "lineage",
            "program_id",
        ]
    )
)

sensitivity_effects = (
    donor_effects[
        donor_effects["sampling"]
        == "sensitivity"
    ]
    .set_index(
        [
            "lineage",
            "program_id",
        ]
    )
)

replication_rows = []

for index, primary_row in primary_effects.iterrows():

    sensitivity_row = sensitivity_effects.loc[
        index
    ]

    same_direction = (
        np.sign(
            primary_row[
                "effect_TLE_minus_Normal"
            ]
        )
        == np.sign(
            sensitivity_row[
                "effect_TLE_minus_Normal"
            ]
        )
    )

    formal_candidate = bool(
        primary_row["evidence_status"]
        == "CORE_REPLICATED"
        and np.isfinite(
            primary_row[
                "formal_primary_fdr"
            ]
        )
        and primary_row[
            "formal_primary_fdr"
        ] < 0.10
        and same_direction
        and primary_row[
            "lodo_direction_stability"
        ] >= 0.75
        and sensitivity_row[
            "lodo_direction_stability"
        ] >= 0.75
    )

    directional_candidate = bool(
        primary_row["evidence_status"]
        == "CORE_REPLICATED"
        and primary_row[
            "exact_p_value"
        ] < 0.10
        and same_direction
        and primary_row[
            "lodo_direction_stability"
        ] >= 0.75
        and sensitivity_row[
            "lodo_direction_stability"
        ] >= 0.75
    )

    replication_rows.append(
        {
            "lineage": index[0],
            "program_id": index[1],
            "evidence_status":
                primary_row[
                    "evidence_status"
                ],

            "primary_effect":
                primary_row[
                    "effect_TLE_minus_Normal"
                ],

            "primary_hedges_g":
                primary_row[
                    "hedges_g"
                ],

            "primary_exact_p":
                primary_row[
                    "exact_p_value"
                ],

            "primary_fdr":
                primary_row[
                    "formal_primary_fdr"
                ],

            "primary_lodo":
                primary_row[
                    "lodo_direction_stability"
                ],

            "sensitivity_effect":
                sensitivity_row[
                    "effect_TLE_minus_Normal"
                ],

            "sensitivity_exact_p":
                sensitivity_row[
                    "exact_p_value"
                ],

            "sensitivity_lodo":
                sensitivity_row[
                    "lodo_direction_stability"
                ],

            "same_direction":
                same_direction,

            "formal_TLE_program_candidate":
                formal_candidate,

            "directional_support_only":
                directional_candidate
                and not formal_candidate,
        }
    )

replication_summary = pd.DataFrame(
    replication_rows
)


program_lock.to_csv(
    OUTPUT_ROOT
    / "freeze"
    / "32B3A_FINAL_PROGRAM_LOCK.csv",
    index=False,
)

top_genes.to_csv(
    OUTPUT_ROOT
    / "program_spectra"
    / "32B3A_TOP100_GENES.csv",
    index=False,
)

donor_effects.to_csv(
    OUTPUT_ROOT
    / "donor_effects"
    / "32B3A_ALL_DONOR_EFFECTS.csv",
    index=False,
)

donor_usage_long.to_csv(
    OUTPUT_ROOT
    / "donor_effects"
    / "32B3A_DONOR_USAGE_LONG.csv",
    index=False,
)

replication_summary.to_csv(
    OUTPUT_ROOT
    / "donor_effects"
    / "32B3A_PROGRAM_TLE_EFFECT_SUMMARY.csv",
    index=False,
)


###############################################################################
# 6. Effect figure
###############################################################################

plot_data = replication_summary[
    replication_summary[
        "evidence_status"
    ] != "EXCLUDED_UNSTABLE"
].copy()

plot_data = plot_data.sort_values(
    [
        "lineage",
        "primary_effect",
    ]
)

figure_height = max(
    6,
    0.32 * len(plot_data),
)

figure, axis = plt.subplots(
    figsize=(9, figure_height)
)

positions = np.arange(
    len(plot_data)
)

axis.axvline(
    0,
    linestyle="--",
    linewidth=1,
)

axis.scatter(
    plot_data["primary_effect"],
    positions,
    label="Primary balanced sample",
    marker="o",
)

axis.scatter(
    plot_data["sensitivity_effect"],
    positions,
    label="Independent cell sample",
    marker="x",
)

axis.set_yticks(
    positions
)

axis.set_yticklabels(
    plot_data["program_id"]
)

axis.set_xlabel(
    "Donor-level CLR effect: TLE − Normal"
)

axis.set_title(
    "Unbiased cNMF program effects at donor level"
)

axis.legend()
figure.tight_layout()

figure.savefig(
    OUTPUT_ROOT
    / "figures"
    / "32B3A_DONOR_PROGRAM_EFFECTS.png",
    dpi=300,
)

figure.savefig(
    OUTPUT_ROOT
    / "figures"
    / "32B3A_DONOR_PROGRAM_EFFECTS.pdf",
)

plt.close(figure)


###############################################################################
# 7. Machine-readable lock
###############################################################################

final_lock = {
    "stage": "32B3A",
    "created_at": now_string(),

    "discovery_completed_before_disease_testing":
        True,

    "disease_labels_used_for_program_discovery":
        False,

    "senescence_gene_sets_used":
        False,

    "formal_program_count": int(
        (
            program_lock[
                "evidence_status"
            ]
            == "CORE_REPLICATED"
        ).sum()
    ),

    "exploratory_program_count": int(
        (
            program_lock[
                "evidence_status"
            ]
            == "EXPLORATORY"
        ).sum()
    ),

    "excluded_program_count": int(
        (
            program_lock[
                "evidence_status"
            ]
            == "EXCLUDED_UNSTABLE"
        ).sum()
    ),

    "formal_TLE_program_candidates": (
        replication_summary.loc[
            replication_summary[
                "formal_TLE_program_candidate"
            ],
            "program_id",
        ].tolist()
    ),

    "directional_support_programs": (
        replication_summary.loc[
            replication_summary[
                "directional_support_only"
            ],
            "program_id",
        ].tolist()
    ),

    "next_stage": (
        "32B3B_INDEPENDENT_BIOLOGICAL_AND_"
        "SENESCENCE_PROGRAM_INTERPRETATION"
    ),
}

with (
    OUTPUT_ROOT
    / "freeze"
    / "32B3A_ANALYSIS_LOCK.json"
).open(
    "w",
    encoding="utf-8",
) as handle:

    json.dump(
        final_lock,
        handle,
        indent=2,
        ensure_ascii=False,
    )


###############################################################################
# 8. Report
###############################################################################

formal_candidates = (
    replication_summary.loc[
        replication_summary[
            "formal_TLE_program_candidate"
        ],
        "program_id",
    ].tolist()
)

directional_candidates = (
    replication_summary.loc[
        replication_summary[
            "directional_support_only"
        ],
        "program_id",
    ].tolist()
)

report_lines = [
    "# 32B-3A Unbiased cNMF Program Freeze",
    "",
    "## Status: COMPLETED",
    "",
    f"- Core replicated programs: "
    f"{final_lock['formal_program_count']}",

    f"- Exploratory programs: "
    f"{final_lock['exploratory_program_count']}",

    f"- Excluded unstable programs: "
    f"{final_lock['excluded_program_count']}",

    "",
    "## Disease-effect analysis",
    "",
    "- Statistical unit: donor.",
    "- Design: 4 Normal versus 4 TLE donors.",
    "- Primary inference: primary balanced cell sample.",
    "- Independent cell sampling: robustness analysis only.",
    "- Exact donor-label permutation was used.",
    "- FDR was calculated across core replicated programs.",
    "- No senescence genes were used in this stage.",
    "",
    "## Formal TLE-associated program candidates",
    "",
    (
        ", ".join(formal_candidates)
        if formal_candidates
        else "None met the pre-specified formal gate."
    ),
    "",
    "## Directional support only",
    "",
    (
        ", ".join(directional_candidates)
        if directional_candidates
        else "None."
    ),
    "",
    "## Next stage",
    "",
    (
        "Annotate frozen unbiased programs using general biological "
        "pathways first, then independently test overlap with the "
        "frozen senescence and competing-response programs."
    ),
]

(
    OUTPUT_ROOT
    / "32B3A_EXECUTIVE_SUMMARY.md"
).write_text(
    "\n".join(report_lines),
    encoding="utf-8",
)


###############################################################################
# 9. Manifest and ZIP
###############################################################################

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
            "size_bytes":
                path.stat().st_size,
            "sha256":
                sha256_file(path),
        }
    )

manifest = pd.DataFrame(
    manifest_rows
)

manifest.to_csv(
    OUTPUT_ROOT
    / "32B3A_FILE_MANIFEST.csv",
    index=False,
)

review_zip = Path(
    str(OUTPUT_ROOT)
    + ".zip"
)

with zipfile.ZipFile(
    review_zip,
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

review_hash = sha256_file(
    review_zip
)

Path(
    str(review_zip)
    + ".sha256.txt"
).write_text(
    f"{review_hash}  {review_zip.name}\n",
    encoding="ascii",
)


###############################################################################
# 10. Final output
###############################################################################

elapsed = time.time() - overall_started

print()
print("=" * 72)
print("32B-3A completed successfully")
print("=" * 72)
print(f"Elapsed: {elapsed:.2f} seconds")
print(f"Core programs: {final_lock['formal_program_count']}")
print(f"Exploratory programs: {final_lock['exploratory_program_count']}")
print(f"Excluded programs: {final_lock['excluded_program_count']}")
print()
print("Formal TLE-associated candidates:")
print(formal_candidates)
print()
print("Directional-support programs:")
print(directional_candidates)
print()
print("Output:")
print(OUTPUT_ROOT)
print()
print("ZIP:")
print(review_zip)
print()
print("SHA256:")
print(review_hash)
print("=" * 72)