from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sys
import zipfile
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT = Path(
    os.environ.get(
        "TLE_PROJECT",
        (
            __import__("os").environ["P5_WORK_ROOT"]
        ),
    )
)

ROOT = (
    PROJECT
    / "99_frozen_releases"
)

CNMF_WORK = (
    PROJECT
    / "04_tables"
    / "32B_CNMf_WORK_20260728_215801"
)

STAMP = datetime.now().strftime(
    "%Y%m%d_%H%M%S"
)

OUT = (
    ROOT
    / f"34A2_DONOR_LEVEL_PROGRAM_MECHANISM_INPUT_AUDIT_{STAMP}"
)


H5AD = {
    "Astrocytes":
        "Astrocytes_primary.h5ad",

    "Microglia_macrophage":
        "Microglia_macrophage_primary.h5ad",

    "Oligodendrocytes":
        "Oligodendrocytes_primary.h5ad",
}


PACKAGES = [
    "anndata",
    "numpy",
    "pandas",
    "scipy",
    "scanpy",
    "statsmodels",
    "decoupler",
    "liana",
]


def sha(
    path: Path,
) -> str:

    digest = hashlib.sha256()

    with path.open(
        "rb"
    ) as handle:

        for block in iter(
            lambda: handle.read(
                1024 * 1024
            ),
            b"",
        ):

            digest.update(
                block
            )

    return digest.hexdigest()


def latest(
    prefix: str,
    required: list[str],
) -> Path:

    roots = sorted(
        [
            path
            for path in ROOT.glob(
                prefix + "*"
            )
            if path.is_dir()
        ],
        key=lambda path:
            path.stat().st_mtime,
        reverse=True,
    )

    for root in roots:

        complete = all(
            len(
                [
                    path
                    for path in root.rglob(
                        name
                    )
                    if path.is_file()
                ]
            ) == 1
            for name in required
        )

        if complete:

            return root.resolve()

    raise FileNotFoundError(
        f"没有找到完整release：{prefix}"
    )


def one(
    root: Path,
    name: str,
) -> Path:

    hits = [
        path
        for path in root.rglob(
            name
        )
        if path.is_file()
    ]

    if len(
        hits
    ) != 1:

        raise RuntimeError(
            f"文件未能唯一定位：{name}\n"
            + "\n".join(
                map(
                    str,
                    hits,
                )
            )
        )

    return hits[
        0
    ].resolve()


def find_h5ad(
    name: str,
) -> Path:

    hits = [
        path
        for path in CNMF_WORK.rglob(
            name
        )
        if path.is_file()
    ]

    if len(
        hits
    ) != 1:

        raise RuntimeError(
            f"{name}找到{len(hits)}个。"
        )

    return hits[
        0
    ].resolve()


def normalize_condition(
    value,
) -> str:

    value = str(
        value
    ).strip().lower()

    if value in {
        "normal",
        "control",
        "ctrl",
        "con",
    }:

        return "Normal"

    if value in {
        "tle",
        "epilepsy",
        "case",
        "disease",
    }:

        return "TLE"

    raise ValueError(
        f"无法识别condition：{value!r}"
    )


def audit_count_matrix(
    matrix,
) -> tuple[
    bool,
    float,
    float,
]:

    if sparse.issparse(
        matrix
    ):

        values = matrix.data

    else:

        values = np.asarray(
            matrix
        ).ravel()

    values = np.asarray(
        values,
        dtype=float,
    )

    values = values[
        np.isfinite(
            values
        )
    ]

    if values.size == 0:

        return (
            False,
            np.nan,
            np.nan,
        )

    integer_fraction = float(
        np.mean(
            np.abs(
                values
                - np.rint(
                    values
                )
            ) < 1e-6
        )
    )

    minimum_value = float(
        values.min()
    )

    count_like = bool(
        minimum_value >= -1e-8
        and integer_fraction >= 0.995
    )

    return (
        count_like,
        integer_fraction,
        minimum_value,
    )


def sum_rows(
    matrix,
    mask: np.ndarray,
) -> np.ndarray:

    subset = matrix[
        mask
    ]

    if sparse.issparse(
        subset
    ):

        return np.asarray(
            subset.sum(
                axis=0
            )
        ).ravel().astype(
            float
        )

    return np.asarray(
        subset,
        dtype=float,
    ).sum(
        axis=0
    )


def package_row(
    name: str,
) -> dict:

    available = (
        importlib.util.find_spec(
            name
        )
        is not None
    )

    version = ""
    import_error = ""

    if available:

        try:

            module = __import__(
                name
            )

            version = str(
                getattr(
                    module,
                    "__version__",
                    "UNKNOWN",
                )
            )

        except Exception as exc:

            import_error = (
                f"{type(exc).__name__}: {exc}"
            )

    return {
        "package":
            name,

        "available":
            bool(
                available
            ),

        "version":
            version,

        "import_error":
            import_error,
    }


###############################################################################
# Locate frozen releases
###############################################################################


B5 = latest(
    "34A1B5_FINAL_CONTRASTIVEVI_CONSENSUS_",
    [
        "34A1B5_FINAL_PROGRAM_CONSENSUS.csv",
        "34A1B5_FINAL_CONSENSUS_LOCK.json",
    ],
)


B2 = latest(
    "34A1B2_FROZEN_PROGRAM_SALIENT_ALIGNMENT_",
    [
        "34A1B2_CELL_FROZEN_PROGRAM_SCORES.csv",
        "34A1B2_PROGRAM_SCORE_COVERAGE.csv",
    ],
)


CONSENSUS_PATH = one(
    B5,
    "34A1B5_FINAL_PROGRAM_CONSENSUS.csv",
)


LOCK_PATH = one(
    B5,
    "34A1B5_FINAL_CONSENSUS_LOCK.json",
)


SCORES_PATH = one(
    B2,
    "34A1B2_CELL_FROZEN_PROGRAM_SCORES.csv",
)


COVERAGE_PATH = one(
    B2,
    "34A1B2_PROGRAM_SCORE_COVERAGE.csv",
)


lock = json.loads(
    LOCK_PATH.read_text(
        encoding="utf-8"
    )
)


if (
    lock.get(
        "status"
    )
    !=
    "34A1B5_FINAL_GLIAL_PROGRAM_CONSENSUS_FROZEN"
):

    raise RuntimeError(
        "34A1B5锁文件状态不正确。"
    )


if (
    lock.get(
        "next_stage"
    )
    !=
    "34A2_DONOR_LEVEL_PROGRAM_MECHANISM_INPUT_AUDIT"
):

    raise RuntimeError(
        "34A1B5未授权进入34A2。"
    )


consensus = pd.read_csv(
    CONSENSUS_PATH
)


scores = pd.read_csv(
    SCORES_PATH
)


coverage = pd.read_csv(
    COVERAGE_PATH
)


required_consensus = {
    "lineage",
    "program_id",
    "final_tier",
    "final_role",
    "exclude_from_core",
}


if (
    required_consensus
    - set(
        consensus.columns
    )
):

    raise RuntimeError(
        "34A1B5共识表缺少必要列。"
    )


required_score_columns = {
    "cell_id",
    "donor_id",
    "condition",
    "lineage",
}


if (
    required_score_columns
    - set(
        scores.columns
    )
):

    raise RuntimeError(
        "34A1B2评分表缺少必要列。"
    )


if scores[
    "cell_id"
].astype(
    str
).duplicated().any():

    raise RuntimeError(
        "评分表存在重复cell_id。"
    )


consensus[
    "program_id"
] = consensus[
    "program_id"
].astype(
    str
)


consensus[
    "lineage"
] = consensus[
    "lineage"
].astype(
    str
)


for program_id in consensus[
    "program_id"
]:

    for score_name in [
        "TOP50_WEIGHTED",
        "TOP100_WEIGHTED",
    ]:

        score_column = (
            f"{program_id}__"
            f"{score_name}"
        )

        if score_column not in scores.columns:

            raise RuntimeError(
                f"缺少评分列：{score_column}"
            )


###############################################################################
# Audit lineages and construct donor-level inputs
###############################################################################


lineage_rows = []
pseudobulk_rows = []
donor_score_rows = []
correlation_rows = []


for (
    lineage,
    filename,
) in H5AD.items():

    print(
        "=" * 78
    )

    print(
        "审计谱系：",
        lineage,
    )

    print(
        "=" * 78
    )

    h5ad_path = find_h5ad(
        filename
    )

    adata = ad.read_h5ad(
        h5ad_path
    )

    if (
        {
            "donor_id",
            "condition",
        }
        - set(
            adata.obs.columns
        )
    ):

        raise RuntimeError(
            f"{lineage} H5AD缺少"
            "donor_id或condition。"
        )

    adata.obs[
        "donor_id"
    ] = adata.obs[
        "donor_id"
    ].astype(
        str
    )

    adata.obs[
        "condition_model"
    ] = [
        normalize_condition(
            value
        )
        for value in adata.obs[
            "condition"
        ]
    ]

    donors = sorted(
        adata.obs[
            "donor_id"
        ].unique()
    )

    donor_map = adata.obs[
        [
            "donor_id",
            "condition_model",
        ]
    ].drop_duplicates()

    group_counts = (
        donor_map
        .groupby(
            "condition_model"
        )[
            "donor_id"
        ]
        .nunique()
        .to_dict()
    )

    if (
        len(
            donors
        ) != 8
        or group_counts
        != {
            "Normal": 4,
            "TLE": 4,
        }
    ):

        raise RuntimeError(
            f"{lineage}不是8名供者且4/4分组："
            f"{group_counts}"
        )

    (
        count_like,
        integer_fraction,
        minimum_value,
    ) = audit_count_matrix(
        adata.X
    )

    if not count_like:

        raise RuntimeError(
            f"{lineage}的X不是可靠原始count矩阵。"
        )

    lineage_scores = scores.loc[
        scores[
            "lineage"
        ].astype(
            str
        ) == lineage
    ].copy()

    score_index = pd.Index(
        lineage_scores[
            "cell_id"
        ].astype(
            str
        )
    )

    positions = score_index.get_indexer(
        adata.obs_names.astype(
            str
        )
    )

    if (
        np.any(
            positions < 0
        )
        or len(
            np.unique(
                positions
            )
        ) != adata.n_obs
    ):

        raise RuntimeError(
            f"{lineage} H5AD与评分表不能一一对应。"
        )

    ordered_scores = (
        lineage_scores
        .iloc[
            positions
        ]
        .reset_index(
            drop=True
        )
    )

    donor_array = adata.obs[
        "donor_id"
    ].to_numpy(
        dtype=str
    )

    condition_array = adata.obs[
        "condition_model"
    ].to_numpy(
        dtype=str
    )

    if not np.array_equal(
        ordered_scores[
            "donor_id"
        ].astype(
            str
        ).to_numpy(),
        donor_array,
    ):

        raise RuntimeError(
            f"{lineage} donor_id对齐失败。"
        )

    ordered_condition = np.asarray(
        [
            normalize_condition(
                value
            )
            for value in ordered_scores[
                "condition"
            ]
        ],
        dtype=str,
    )

    if not np.array_equal(
        ordered_condition,
        condition_array,
    ):

        raise RuntimeError(
            f"{lineage} condition对齐失败。"
        )

    donor_vectors = []

    for donor_id in donors:

        donor_mask = (
            donor_array
            == donor_id
        )

        donor_conditions = pd.unique(
            condition_array[
                donor_mask
            ]
        )

        if len(
            donor_conditions
        ) != 1:

            raise RuntimeError(
                f"{lineage}/{donor_id}"
                "对应多个condition。"
            )

        donor_vector = sum_rows(
            adata.X,
            donor_mask,
        )

        donor_vectors.append(
            donor_vector
        )

        pseudobulk_rows.append(
            {
                "lineage":
                    lineage,

                "donor_id":
                    donor_id,

                "condition":
                    donor_conditions[
                        0
                    ],

                "n_cells":
                    int(
                        donor_mask.sum()
                    ),

                "library_size":
                    float(
                        donor_vector.sum()
                    ),

                "detected_genes":
                    int(
                        np.sum(
                            donor_vector > 0
                        )
                    ),

                "genes_with_counts_ge_10":
                    int(
                        np.sum(
                            donor_vector >= 10
                        )
                    ),
            }
        )

    donor_matrix = np.vstack(
        donor_vectors
    )

    donors_detected_per_gene = np.sum(
        donor_matrix > 0,
        axis=0,
    )

    cells_per_donor = (
        adata.obs
        .groupby(
            "donor_id",
            observed=True,
        )
        .size()
    )

    lineage_rows.append(
        {
            "lineage":
                lineage,

            "h5ad_path":
                str(
                    h5ad_path
                ),

            "h5ad_sha256":
                sha(
                    h5ad_path
                ),

            "n_cells":
                int(
                    adata.n_obs
                ),

            "n_genes":
                int(
                    adata.n_vars
                ),

            "donors":
                8,

            "normal_donors":
                4,

            "tle_donors":
                4,

            "count_like":
                count_like,

            "integer_fraction":
                integer_fraction,

            "minimum_value":
                minimum_value,

            "genes_detected_in_all_8_donors":
                int(
                    np.sum(
                        donors_detected_per_gene
                        == 8
                    )
                ),

            "genes_detected_in_at_least_4_donors":
                int(
                    np.sum(
                        donors_detected_per_gene
                        >= 4
                    )
                ),

            "minimum_cells_per_donor":
                int(
                    cells_per_donor.min()
                ),

            "maximum_cells_per_donor":
                int(
                    cells_per_donor.max()
                ),
        }
    )

    lineage_programs = consensus.loc[
        consensus[
            "lineage"
        ] == lineage
    ].copy()

    donor_score_wide_rows = []

    for donor_id in donors:

        donor_mask = (
            donor_array
            == donor_id
        )

        donor_condition = pd.unique(
            condition_array[
                donor_mask
            ]
        )[
            0
        ]

        wide_row = {
            "donor_id":
                donor_id,

            "condition":
                donor_condition,
        }

        for program_id in lineage_programs[
            "program_id"
        ]:

            for score_name in [
                "TOP50_WEIGHTED",
                "TOP100_WEIGHTED",
            ]:

                score_column = (
                    f"{program_id}__"
                    f"{score_name}"
                )

                values = pd.to_numeric(
                    ordered_scores.loc[
                        donor_mask,
                        score_column,
                    ],
                    errors="coerce",
                ).to_numpy(
                    dtype=float
                )

                if not np.isfinite(
                    values
                ).all():

                    raise RuntimeError(
                        f"{lineage}/{score_column}"
                        "含非有限值。"
                    )

                mean_score = float(
                    values.mean()
                )

                donor_score_rows.append(
                    {
                        "lineage":
                            lineage,

                        "program_id":
                            program_id,

                        "score_definition":
                            score_name,

                        "donor_id":
                            donor_id,

                        "condition":
                            donor_condition,

                        "n_cells":
                            int(
                                donor_mask.sum()
                            ),

                        "mean_score":
                            mean_score,

                        "median_score":
                            float(
                                np.median(
                                    values
                                )
                            ),

                        "sd_score":
                            float(
                                np.std(
                                    values,
                                    ddof=1,
                                )
                            ),
                    }
                )

                wide_row[
                    score_column
                ] = mean_score

        donor_score_wide_rows.append(
            wide_row
        )

    donor_score_wide = pd.DataFrame(
        donor_score_wide_rows
    )

    score_columns = [
        column
        for column in donor_score_wide.columns
        if column
        not in {
            "donor_id",
            "condition",
        }
    ]

    for (
        first_index,
        first_score,
    ) in enumerate(
        score_columns
    ):

        for second_score in score_columns[
            first_index + 1:
        ]:

            correlation_rows.append(
                {
                    "lineage":
                        lineage,

                    "first_score":
                        first_score,

                    "second_score":
                        second_score,

                    "n_donors":
                        8,

                    "pearson_r":
                        float(
                            donor_score_wide[
                                first_score
                            ].corr(
                                donor_score_wide[
                                    second_score
                                ],
                                method="pearson",
                            )
                        ),

                    "spearman_rho":
                        float(
                            donor_score_wide[
                                first_score
                            ].corr(
                                donor_score_wide[
                                    second_score
                                ],
                                method="spearman",
                            )
                        ),

                    "descriptive_only":
                        True,
                }
            )

    del adata


###############################################################################
# Assemble audits
###############################################################################


lineage_audit = pd.DataFrame(
    lineage_rows
)


pseudobulk_audit = pd.DataFrame(
    pseudobulk_rows
)


donor_scores = pd.DataFrame(
    donor_score_rows
)


correlations = pd.DataFrame(
    correlation_rows
)


package_audit = pd.DataFrame(
    [
        package_row(
            package_name
        )
        for package_name in PACKAGES
    ]
)


score_audit_rows = []


for row in consensus.itertuples(
    index=False
):

    for score_name in [
        "TOP50_WEIGHTED",
        "TOP100_WEIGHTED",
    ]:

        subset = donor_scores.loc[
            (
                donor_scores[
                    "program_id"
                ] == row.program_id
            )
            &
            (
                donor_scores[
                    "score_definition"
                ] == score_name
            )
        ]

        if len(
            subset
        ) != 8:

            raise RuntimeError(
                f"{row.program_id}/{score_name}"
                "不是8名供者。"
            )

        donor_sd = float(
            subset[
                "mean_score"
            ].std(
                ddof=1
            )
        )

        score_audit_rows.append(
            {
                "lineage":
                    row.lineage,

                "program_id":
                    row.program_id,

                "final_tier":
                    row.final_tier,

                "final_role":
                    row.final_role,

                "score_definition":
                    score_name,

                "donors":
                    8,

                "normal_donors":
                    int(
                        (
                            subset[
                                "condition"
                            ] == "Normal"
                        ).sum()
                    ),

                "tle_donors":
                    int(
                        (
                            subset[
                                "condition"
                            ] == "TLE"
                        ).sum()
                    ),

                "donor_mean_score_sd":
                    donor_sd,

                "finite_scores":
                    bool(
                        np.isfinite(
                            subset[
                                "mean_score"
                            ]
                        ).all()
                    ),

                "nonzero_donor_variance":
                    bool(
                        donor_sd > 1e-12
                    ),
            }
        )


score_audit = pd.DataFrame(
    score_audit_rows
)


primary_coverage = coverage.loc[
    (
        coverage[
            "program_id"
        ].isin(
            [
                "Oligodendrocytes__P5",
                "Astrocytes__P6",
            ]
        )
    )
    &
    (
        coverage[
            "score_definition"
        ] == "TOP100_WEIGHTED"
    )
].copy()


if len(
    primary_coverage
) != 2:

    raise RuntimeError(
        "无法在34A1B2覆盖表中定位"
        "P5和Astro P6的Top100记录。"
    )


input_pass = bool(
    len(
        lineage_audit
    ) == 3
    and lineage_audit[
        "count_like"
    ].all()
    and lineage_audit[
        "genes_detected_in_at_least_4_donors"
    ].ge(
        5000
    ).all()
    and score_audit[
        "finite_scores"
    ].all()
    and score_audit[
        "nonzero_donor_variance"
    ].all()
    and primary_coverage[
        "coverage_fraction"
    ].ge(
        0.80
    ).all()
)


required_environment = package_audit.loc[
    package_audit[
        "package"
    ].isin(
        [
            "decoupler",
            "statsmodels",
        ]
    )
].copy()


environment_ready = bool(
    len(
        required_environment
    ) == 2
    and required_environment[
        "available"
    ].all()
    and required_environment[
        "import_error"
    ].fillna(
        ""
    ).eq(
        ""
    ).all()
)


if (
    input_pass
    and environment_ready
):

    STATUS = (
        "GO_FOR_34A3_DONOR_LEVEL_"
        "INTRINSIC_TF_ACTIVITY"
    )

    NEXT_STAGE = (
        "34A3_DONOR_LEVEL_"
        "INTRINSIC_TF_ACTIVITY"
    )

elif input_pass:

    STATUS = (
        "INPUTS_PASS_MECHANISM_ENVIRONMENT_"
        "PREPARATION_REQUIRED"
    )

    NEXT_STAGE = (
        "34A2B_MECHANISM_ENVIRONMENT_"
        "AND_RESOURCE_FREEZE"
    )

else:

    STATUS = (
        "HOLD_FOR_34A2_INPUT_REPAIR"
    )

    NEXT_STAGE = (
        "34A2_INPUT_REPAIR"
    )


plan = pd.DataFrame(
    [
        [
            "PRIMARY_INTRINSIC_TF_ACTIVITY",
            "Oligodendrocytes",
            "Oligodendrocytes__P5",
            "Oligodendrocytes__P6",
            "",
            "donor",
            "continuous TOP100 frozen program score",
            True,
        ],
        [
            "SECONDARY_INTRINSIC_TF_ACTIVITY",
            "Astrocytes",
            "Astrocytes__P6",
            "Astrocytes__P4",
            "Astrocytes__P2",
            "donor",
            "continuous TOP100 frozen program score",
            True,
        ],
        [
            "MYELOID_CONDITIONAL_SENSITIVITY",
            "Microglia_macrophage",
            "Microglia_macrophage__P4",
            "",
            "Microglia_macrophage__P2",
            "donor",
            "continuous TOP100 frozen program score",
            True,
        ],
        [
            "CELL_CELL_COMMUNICATION_LATER",
            "multi-lineage",
            "Oligodendrocytes__P5",
            (
                "Astrocytes__P6;"
                "Microglia_macrophage__P4"
            ),
            "Astrocytes__P2",
            "donor",
            "receiver-program association",
            True,
        ],
    ],
    columns=[
        "analysis_component",
        "lineage",
        "primary_program",
        "comparator",
        "negative_control",
        "inferential_unit",
        "state_representation",
        "authorized",
    ],
)


###############################################################################
# Freeze outputs
###############################################################################


OUT.mkdir(
    parents=True,
    exist_ok=False,
)


for folder in [
    "01_input_audit",
    "02_donor_scores",
    "03_pseudobulk_audit",
    "04_environment",
    "05_plan",
    "06_freeze",
    "07_logs",
]:

    (
        OUT
        / folder
    ).mkdir(
        parents=True,
        exist_ok=True,
    )


lineage_audit.to_csv(
    OUT
    / "01_input_audit"
    / "34A2_LINEAGE_INPUT_AUDIT.csv",
    index=False,
)


score_audit.to_csv(
    OUT
    / "01_input_audit"
    / "34A2_PROGRAM_SCORE_INPUT_AUDIT.csv",
    index=False,
)


primary_coverage.to_csv(
    OUT
    / "01_input_audit"
    / "34A2_PRIMARY_PROGRAM_COVERAGE_AUDIT.csv",
    index=False,
)


donor_scores.to_csv(
    OUT
    / "02_donor_scores"
    / "34A2_DONOR_PROGRAM_SCORE_TABLE.csv",
    index=False,
)


correlations.to_csv(
    OUT
    / "02_donor_scores"
    / "34A2_DONOR_PROGRAM_SCORE_CORRELATIONS.csv",
    index=False,
)


pseudobulk_audit.to_csv(
    OUT
    / "03_pseudobulk_audit"
    / "34A2_DONOR_PSEUDOBULK_INPUT_AUDIT.csv",
    index=False,
)


package_audit.to_csv(
    OUT
    / "04_environment"
    / "34A2_PACKAGE_CAPABILITY_AUDIT.csv",
    index=False,
)


plan.to_csv(
    OUT
    / "05_plan"
    / "34A2_MECHANISM_ANALYSIS_PLAN.csv",
    index=False,
)


freeze = {
    "stage":
        "34A2",

    "created_at":
        datetime.now().isoformat(),

    "status":
        STATUS,

    "source_34A1B5_release":
        str(
            B5
        ),

    "source_34A1B2_release":
        str(
            B2
        ),

    "primary_program":
        "Oligodendrocytes__P5",

    "secondary_program":
        "Astrocytes__P6",

    "inferential_unit":
        "donor",

    "state_representation":
        "continuous frozen program score",

    "input_pass":
        input_pass,

    "mechanism_environment_ready":
        environment_ready,

    "new_model_fitted":
        False,

    "new_disease_test_performed":
        False,

    "cell_level_p_values_authorized":
        False,

    "discrete_state_label_authorized":
        False,

    "senescent_cell_label_authorized":
        False,

    "causal_claim_authorized":
        False,

    "next_stage":
        NEXT_STAGE,
}


(
    OUT
    / "06_freeze"
    / "34A2_MECHANISM_INPUT_AUDIT_LOCK.json"
).write_text(
    json.dumps(
        freeze,
        ensure_ascii=False,
        indent=2,
    ),
    encoding="utf-8",
)


summary_lines = [
    "# 34A2 donor-level program mechanism input audit",
    "",
    f"## Status: {STATUS}",
    "",
    "- Donor remains the inferential unit.",
    "- Programs remain continuous frozen scores.",
    (
        "- No discrete-state, senescent-cell "
        "or causal label is authorized."
    ),
    "",
    "## Next stage",
    "",
    NEXT_STAGE,
]


(
    OUT
    / "34A2_EXECUTIVE_SUMMARY.md"
).write_text(
    "\n".join(
        summary_lines
    ),
    encoding="utf-8",
)


provenance = {
    "stage":
        "34A2",

    "python":
        sys.version,

    "inputs":
        {
            str(
                CONSENSUS_PATH
            ):
                sha(
                    CONSENSUS_PATH
                ),

            str(
                LOCK_PATH
            ):
                sha(
                    LOCK_PATH
                ),

            str(
                SCORES_PATH
            ):
                sha(
                    SCORES_PATH
                ),

            str(
                COVERAGE_PATH
            ):
                sha(
                    COVERAGE_PATH
                ),
        },
}


(
    OUT
    / "07_logs"
    / "34A2_INPUT_PROVENANCE.json"
).write_text(
    json.dumps(
        provenance,
        ensure_ascii=False,
        indent=2,
    ),
    encoding="utf-8",
)


output_files = sorted(
    path
    for path in OUT.rglob(
        "*"
    )
    if path.is_file()
)


pd.DataFrame(
    [
        {
            "relative_path":
                path.relative_to(
                    OUT
                ).as_posix(),

            "size_bytes":
                path.stat().st_size,

            "sha256":
                sha(
                    path
                ),
        }
        for path in output_files
    ]
).to_csv(
    OUT
    / "34A2_FILE_MANIFEST.csv",
    index=False,
)


ZIP_PATH = Path(
    str(
        OUT
    )
    + ".zip"
)


with zipfile.ZipFile(
    ZIP_PATH,
    "w",
    zipfile.ZIP_DEFLATED,
    compresslevel=9,
) as archive:

    for path in sorted(
        OUT.rglob(
            "*"
        )
    ):

        if path.is_file():

            archive.write(
                path,
                path.relative_to(
                    OUT
                ).as_posix(),
            )


ZIP_HASH = sha(
    ZIP_PATH
)


Path(
    str(
        ZIP_PATH
    )
    + ".sha256.txt"
).write_text(
    (
        f"{ZIP_HASH}  "
        f"{ZIP_PATH.name}\n"
    ),
    encoding="ascii",
)


print(
    "=" * 78
)

print(
    "34A2 completed"
)

print(
    "状态：",
    STATUS,
)

print(
    "谱系审计：",
    len(
        lineage_audit
    ),
    "/ 3",
)

print(
    "供者程序评分：",
    len(
        donor_scores
    ),
    "行",
)

print(
    "机制环境已就绪：",
    environment_ready,
)

print(
    "下一阶段：",
    NEXT_STAGE,
)

print(
    "ZIP：",
    ZIP_PATH,
)

print(
    "SHA256：",
    ZIP_HASH,
)

print(
    "=" * 78
)
