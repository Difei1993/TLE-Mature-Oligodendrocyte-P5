from __future__ import annotations

import hashlib
import json
import sys
import zipfile
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


###############################################################################
# 34A-1B2
# Frozen cNMF program alignment with contrastiveVI latent subspaces
#
# Primary principles:
#   1. cNMF programs are frozen before contrastiveVI interpretation.
#   2. Individual latent dimensions are not compared across seeds.
#   3. Alignment is assessed against the complete salient/background subspace.
#   4. Donor remains the biological inferential unit.
#   5. No cell-level disease P value is calculated.
#   6. No final senescence/inflammation/metabolism classification is made here.
###############################################################################

PROJECT = Path(
    __import__("os").environ["P5_WORK_ROOT"]
)

RELEASE_ROOT = (
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

OUT_ROOT = (
    RELEASE_ROOT
    / f"34A1B2_FROZEN_PROGRAM_SALIENT_ALIGNMENT_{STAMP}"
)


LINEAGES = {
    "Astrocytes": {
        "h5ad_name":
            "Astrocytes_primary.h5ad",
    },

    "Microglia_macrophage": {
        "h5ad_name":
            "Microglia_macrophage_primary.h5ad",
    },

    "Oligodendrocytes": {
        "h5ad_name":
            "Oligodendrocytes_primary.h5ad",
    },
}


SEEDS = [
    340111,
    340112,
    340113,
]


SCORE_DEFINITIONS = {
    "TOP50_WEIGHTED": 50,
    "TOP100_WEIGHTED": 100,
}


CENTERING_MODES = [
    "CELL_GLOBAL",
    "CELL_CONDITION_CENTERED",
    "CELL_DONOR_CENTERED",
    "DONOR_CENTROID_GLOBAL",
    "DONOR_CENTROID_CONDITION_CENTERED",
]


###############################################################################
# Utility functions
###############################################################################

def sha256_file(
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


def latest_complete_release(
    prefix: str,
    required_files: list[str],
) -> Path:

    candidates = sorted(
        [
            path
            for path in RELEASE_ROOT.glob(
                prefix + "*"
            )
            if path.is_dir()
        ],
        key=lambda path:
            path.stat().st_mtime,
        reverse=True,
    )

    for root in candidates:

        complete = all(
            len(
                [
                    path
                    for path in root.rglob(
                        filename
                    )
                    if path.is_file()
                ]
            ) == 1
            for filename in required_files
        )

        if complete:
            return root.resolve()

    raise FileNotFoundError(
        "没有找到完整release："
        f"{prefix}; required={required_files}"
    )


def one_file(
    root: Path,
    filename: str,
) -> Path:

    hits = [
        path
        for path in root.rglob(
            filename
        )
        if path.is_file()
    ]

    if len(hits) != 1:

        raise RuntimeError(
            f"文件未能唯一定位：{filename}\n"
            + "\n".join(
                map(
                    str,
                    hits,
                )
            )
        )

    return hits[0].resolve()


def find_unique_file(
    root: Path,
    filename: str,
) -> Path:

    hits = [
        path
        for path in root.rglob(
            filename
        )
        if path.is_file()
    ]

    if len(hits) != 1:

        raise RuntimeError(
            f"{filename}找到{len(hits)}个：\n"
            + "\n".join(
                map(
                    str,
                    hits,
                )
            )
        )

    return hits[0].resolve()


def as_boolean(
    values: pd.Series,
) -> pd.Series:

    if pd.api.types.is_bool_dtype(
        values
    ):
        return values.fillna(
            False
        )

    return (
        values.astype(str)
        .str.strip()
        .str.lower()
        .isin(
            [
                "true",
                "t",
                "1",
                "yes",
            ]
        )
    )


def normalize_gene(
    value,
) -> str:

    if value is None:
        return ""

    value = str(
        value
    ).strip().upper()

    if value in {
        "",
        "NA",
        "N/A",
        "NONE",
        "NULL",
        "---",
    }:
        return ""

    return value


def program_lineage(
    program_id: str,
) -> str | None:

    if program_id.startswith(
        "Astrocytes__"
    ):
        return "Astrocytes"

    if program_id.startswith(
        "Microglia_macrophage__"
    ):
        return "Microglia_macrophage"

    if program_id.startswith(
        "Oligodendrocytes__"
    ):
        return "Oligodendrocytes"

    return None


def find_first_column(
    table: pd.DataFrame,
    candidates: list[str],
) -> str:

    normalized_lookup = {
        str(column).strip().lower():
            str(column)
        for column in table.columns
    }

    for candidate in candidates:

        normalized = (
            candidate.strip().lower()
        )

        if normalized in normalized_lookup:

            return normalized_lookup[
                normalized
            ]

    raise RuntimeError(
        "未找到必要列。候选："
        + ", ".join(
            candidates
        )
        + "; 实际列："
        + ", ".join(
            map(
                str,
                table.columns,
            )
        )
    )


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


def group_center(
    matrix: np.ndarray,
    groups: np.ndarray,
) -> np.ndarray:

    matrix = np.asarray(
        matrix,
        dtype=float,
    )

    if matrix.ndim == 1:

        matrix = matrix.reshape(
            -1,
            1,
        )

    groups = np.asarray(
        groups
    )

    if len(groups) != matrix.shape[0]:

        raise ValueError(
            "分组标签数与矩阵行数不一致。"
        )

    centered = matrix.copy()

    for group in pd.unique(
        groups
    ):

        mask = groups == group

        centered[
            mask,
            :,
        ] -= centered[
            mask,
            :,
        ].mean(
            axis=0,
            keepdims=True,
        )

    return centered


def global_center(
    matrix: np.ndarray,
) -> np.ndarray:

    matrix = np.asarray(
        matrix,
        dtype=float,
    )

    if matrix.ndim == 1:

        matrix = matrix.reshape(
            -1,
            1,
        )

    return (
        matrix
        - matrix.mean(
            axis=0,
            keepdims=True,
        )
    )


def linear_cka(
    first: np.ndarray,
    second: np.ndarray,
) -> float:

    first = global_center(
        first
    )

    second = global_center(
        second
    )

    if first.shape[0] != second.shape[0]:

        raise ValueError(
            "CKA两个输入的行数不一致。"
        )

    first_cross_second = (
        first.T
        @ second
    )

    first_self = (
        first.T
        @ first
    )

    second_self = (
        second.T
        @ second
    )

    numerator = float(
        np.linalg.norm(
            first_cross_second,
            ord="fro",
        ) ** 2
    )

    denominator = float(
        np.linalg.norm(
            first_self,
            ord="fro",
        )
        *
        np.linalg.norm(
            second_self,
            ord="fro",
        )
    )

    if (
        not np.isfinite(
            denominator
        )
        or denominator <= 0
    ):
        return np.nan

    return numerator / denominator


def donor_centroids(
    latent: np.ndarray,
    program_score: np.ndarray,
    donor_id: np.ndarray,
    condition: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:

    latent = np.asarray(
        latent,
        dtype=float,
    )

    program_score = np.asarray(
        program_score,
        dtype=float,
    )

    donor_id = np.asarray(
        donor_id
    )

    condition = np.asarray(
        condition
    )

    donor_order = sorted(
        pd.unique(
            donor_id
        )
    )

    latent_rows = []
    score_rows = []
    condition_rows = []

    for donor in donor_order:

        mask = donor_id == donor

        donor_conditions = pd.unique(
            condition[
                mask
            ]
        )

        if len(
            donor_conditions
        ) != 1:

            raise RuntimeError(
                f"供者{donor}对应多个condition。"
            )

        latent_rows.append(
            latent[
                mask,
                :,
            ].mean(
                axis=0
            )
        )

        score_rows.append(
            float(
                program_score[
                    mask
                ].mean()
            )
        )

        condition_rows.append(
            donor_conditions[
                0
            ]
        )

    return (
        np.vstack(
            latent_rows
        ),
        np.asarray(
            score_rows,
            dtype=float,
        ),
        np.asarray(
            donor_order,
            dtype=str,
        ),
        np.asarray(
            condition_rows,
            dtype=str,
        ),
    )


def prepare_for_mode(
    latent: np.ndarray,
    program_score: np.ndarray,
    donor_id: np.ndarray,
    condition: np.ndarray,
    mode: str,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:

    latent = np.asarray(
        latent,
        dtype=float,
    )

    program_score = np.asarray(
        program_score,
        dtype=float,
    )

    if mode == "CELL_GLOBAL":

        return (
            latent,
            program_score,
        )

    if mode == "CELL_CONDITION_CENTERED":

        return (
            group_center(
                latent,
                condition,
            ),
            group_center(
                program_score,
                condition,
            ).ravel(),
        )

    if mode == "CELL_DONOR_CENTERED":

        return (
            group_center(
                latent,
                donor_id,
            ),
            group_center(
                program_score,
                donor_id,
            ).ravel(),
        )

    (
        donor_latent,
        donor_score,
        _,
        donor_condition,
    ) = donor_centroids(
        latent=latent,
        program_score=program_score,
        donor_id=donor_id,
        condition=condition,
    )

    if mode == "DONOR_CENTROID_GLOBAL":

        return (
            donor_latent,
            donor_score,
        )

    if (
        mode
        ==
        "DONOR_CENTROID_CONDITION_CENTERED"
    ):

        return (
            group_center(
                donor_latent,
                donor_condition,
            ),
            group_center(
                donor_score,
                donor_condition,
            ).ravel(),
        )

    raise ValueError(
        f"未知centering mode：{mode}"
    )


def read_program_spectrum(
    path: Path,
) -> pd.DataFrame:

    table = pd.read_csv(
        path
    )

    gene_column = find_first_column(
        table,
        [
            "gene",
            "human_gene",
            "gene_symbol",
            "symbol",
        ],
    )

    score_column = find_first_column(
        table,
        [
            "meta_program_score",
            "program_score",
            "score",
            "weight",
        ],
    )

    table = table[
        [
            gene_column,
            score_column,
        ]
    ].copy()

    table.columns = [
        "gene",
        "meta_program_score",
    ]

    table[
        "gene"
    ] = table[
        "gene"
    ].map(
        normalize_gene
    )

    table[
        "meta_program_score"
    ] = pd.to_numeric(
        table[
            "meta_program_score"
        ],
        errors="coerce",
    )

    table = table.loc[
        table[
            "gene"
        ].ne("")
        &
        np.isfinite(
            table[
                "meta_program_score"
            ]
        )
    ].copy()

    table = (
        table
        .sort_values(
            "meta_program_score",
            ascending=False,
        )
        .drop_duplicates(
            "gene",
            keep="first",
        )
        .reset_index(
            drop=True
        )
    )

    if table.empty:

        raise RuntimeError(
            f"程序谱为空：{path}"
        )

    return table


def calculate_program_score(
    adata,
    normalized_var_names: list[str],
    gene_to_index: dict[str, int],
    spectrum: pd.DataFrame,
    top_n: int,
) -> tuple[
    np.ndarray,
    pd.DataFrame,
]:

    selected = spectrum.head(
        top_n
    ).copy()

    selected[
        "detected"
    ] = selected[
        "gene"
    ].isin(
        gene_to_index
    )

    detected = selected.loc[
        selected[
            "detected"
        ]
    ].copy()

    detected[
        "matrix_index"
    ] = [
        gene_to_index[
            gene
        ]
        for gene in detected[
            "gene"
        ]
    ]

    if len(
        detected
    ) < 10:

        raise RuntimeError(
            f"Top{top_n}检测到的基因不足10："
            f"{len(detected)}"
        )

    indices = detected[
        "matrix_index"
    ].to_numpy(
        dtype=int
    )

    if sparse.issparse(
        adata.X
    ):

        total_counts = np.asarray(
            adata.X.sum(
                axis=1
            )
        ).ravel()

        matrix = (
            adata.X[
                :,
                indices,
            ]
            .toarray()
            .astype(
                np.float64
            )
        )

    else:

        full_matrix = np.asarray(
            adata.X,
            dtype=np.float64,
        )

        total_counts = full_matrix.sum(
            axis=1
        )

        matrix = full_matrix[
            :,
            indices,
        ]

    if np.any(
        total_counts <= 0
    ):

        raise RuntimeError(
            "存在总count为0的细胞。"
        )

    matrix = np.log1p(
        matrix
        /
        total_counts[
            :,
            None,
        ]
        *
        10000.0
    )

    gene_means = matrix.mean(
        axis=0
    )

    gene_sds = matrix.std(
        axis=0,
        ddof=0,
    )

    variable_mask = (
        np.isfinite(
            gene_sds
        )
        &
        (
            gene_sds
            > 1e-12
        )
    )

    matrix = matrix[
        :,
        variable_mask,
    ]

    gene_means = gene_means[
        variable_mask
    ]

    gene_sds = gene_sds[
        variable_mask
    ]

    detected = detected.iloc[
        np.flatnonzero(
            variable_mask
        )
    ].copy()

    if len(
        detected
    ) < 10:

        raise RuntimeError(
            f"Top{top_n}去除零方差基因后"
            f"不足10：{len(detected)}"
        )

    z_matrix = (
        matrix
        - gene_means[
            None,
            :,
        ]
    ) / gene_sds[
        None,
        :,
    ]

    weights = np.maximum(
        detected[
            "meta_program_score"
        ].to_numpy(
            dtype=float
        ),
        0.0,
    )

    if (
        not np.isfinite(
            weights
        ).all()
        or weights.sum() <= 0
    ):

        weights = np.ones(
            len(
                detected
            ),
            dtype=float,
        )

    weights = (
        weights
        /
        weights.sum()
    )

    score = (
        z_matrix
        @ weights
    )

    detected[
        "normalized_weight"
    ] = weights

    detected[
        "top_n_definition"
    ] = top_n

    detected[
        "genes_in_frozen_top_set"
    ] = len(
        selected
    )

    detected[
        "genes_detected_and_variable"
    ] = len(
        detected
    )

    return (
        score.astype(
            np.float64
        ),
        detected,
    )


###############################################################################
# Locate upstream releases
###############################################################################

A0_ROOT = latest_complete_release(
    "34A0_GLIAL_STATE_CONSENSUS_INPUT_AUDIT_",
    [
        "34A0_MASTER_PROGRAM_EVIDENCE_MATRIX.csv",
        "34A0_CONTRASTIVEVI_TARGET_REGISTRY.csv",
        "34A0_STATE_CONSENSUS_LOCK.json",
    ],
)


B3A_ROOT = latest_complete_release(
    "32B3A_UNBIASED_PROGRAM_FREEZE_",
    [
        "32B3A_FINAL_PROGRAM_LOCK.csv",
        "32B3A_ANALYSIS_LOCK.json",
    ],
)


B1_ROOT = latest_complete_release(
    "34A1B1_CONTRASTIVEVI_MULTISEED_STABILITY_",
    [
        "34A1B1_MULTISEED_STABILITY_LOCK.json",
        "34A1B1_LINEAGE_STABILITY_CLASSIFICATION.csv",
        "34A1B1_MODEL_RUN_SUMMARY.csv",
    ],
)


MASTER_PATH = one_file(
    A0_ROOT,
    "34A0_MASTER_PROGRAM_EVIDENCE_MATRIX.csv",
)


TARGET_PATH = one_file(
    A0_ROOT,
    "34A0_CONTRASTIVEVI_TARGET_REGISTRY.csv",
)


A0_LOCK_PATH = one_file(
    A0_ROOT,
    "34A0_STATE_CONSENSUS_LOCK.json",
)


B1_LOCK_PATH = one_file(
    B1_ROOT,
    "34A1B1_MULTISEED_STABILITY_LOCK.json",
)


STABILITY_PATH = one_file(
    B1_ROOT,
    "34A1B1_LINEAGE_STABILITY_CLASSIFICATION.csv",
)


with A0_LOCK_PATH.open(
    "r",
    encoding="utf-8",
) as handle:

    a0_lock = json.load(
        handle
    )


with B1_LOCK_PATH.open(
    "r",
    encoding="utf-8",
) as handle:

    b1_lock = json.load(
        handle
    )


if (
    a0_lock.get(
        "status"
    )
    !=
    "GO_FOR_34A1_CONTRASTIVEVI_TARGETED_THREE_LINEAGE_ANALYSIS"
):

    raise RuntimeError(
        "34A0没有授权contrastiveVI分析。"
    )


if not str(
    b1_lock.get(
        "status",
        "",
    )
).startswith(
    (
        "GO_FOR_34A1B2",
        "CONDITIONAL_FOR_34A1B2",
    )
):

    raise RuntimeError(
        "34A1B1没有授权进入程序对齐。"
    )


if int(
    b1_lock.get(
        "models_completed",
        0,
    )
) != 9:

    raise RuntimeError(
        "34A1B1完成模型数不是9。"
    )


###############################################################################
# Read frozen program registries
###############################################################################

master = pd.read_csv(
    MASTER_PATH
)


targets = pd.read_csv(
    TARGET_PATH
)


stability = pd.read_csv(
    STABILITY_PATH
)


required_master_columns = {
    "program_id",
    "passes_content_qc_for_interpretation",
}


missing_master_columns = (
    required_master_columns
    - set(
        master.columns
    )
)


if missing_master_columns:

    raise RuntimeError(
        "34A0 master缺少必要列："
        + ", ".join(
            sorted(
                missing_master_columns
            )
        )
    )


master[
    "passes_content_qc_for_interpretation"
] = as_boolean(
    master[
        "passes_content_qc_for_interpretation"
    ]
)


master[
    "analysis_lineage"
] = master[
    "program_id"
].astype(
    str
).map(
    program_lineage
)


eligible_programs = master.loc[
    master[
        "passes_content_qc_for_interpretation"
    ]
    &
    master[
        "analysis_lineage"
    ].isin(
        LINEAGES
    )
].copy()


###############################################################################
# Audit eligible programs in the three prespecified lineages
###############################################################################

if eligible_programs.empty:

    raise RuntimeError(
        "三个预设谱系中没有内容QC可解释程序。"
    )


if eligible_programs[
    "program_id"
].duplicated().any():

    duplicated_programs = (
        eligible_programs.loc[
            eligible_programs[
                "program_id"
            ].duplicated(
                keep=False
            ),
            "program_id",
        ]
        .astype(str)
        .tolist()
    )

    raise RuntimeError(
        "内容QC程序存在重复program_id："
        + ", ".join(
            duplicated_programs
        )
    )


eligible_count_by_lineage = (
    eligible_programs
    .groupby(
        "analysis_lineage"
    )[
        "program_id"
    ]
    .nunique()
    .to_dict()
)


missing_eligible_lineages = sorted(
    set(
        LINEAGES
    )
    - set(
        eligible_count_by_lineage
    )
)


if missing_eligible_lineages:

    raise RuntimeError(
        "以下预设谱系没有内容QC可解释程序："
        + ", ".join(
            missing_eligible_lineages
        )
    )


print(
    "三个预设谱系的内容QC可解释程序数：",
    len(
        eligible_programs
    ),
)


print(
    "各谱系程序数：",
    eligible_count_by_lineage,
)


target_columns = [
    column
    for column in [
        "program_id",
        "priority",
        "target_role",
        "contrastivevi_question",
    ]
    if column in targets.columns
]


target_metadata = targets[
    target_columns
].drop_duplicates(
    "program_id"
)


eligible_programs = (
    eligible_programs
    .merge(
        target_metadata,
        on="program_id",
        how="left",
        validate="one_to_one",
    )
)


eligible_programs[
    "pre_registered_target"
] = eligible_programs[
    "program_id"
].isin(
    set(
        targets[
            "program_id"
        ].astype(
            str
        )
    )
)


###############################################################################
# Prepare outputs
###############################################################################

OUT_ROOT.mkdir(
    parents=True,
    exist_ok=False,
)


for subdirectory in [
    "01_program_coverage",
    "02_cell_program_scores",
    "03_seed_level_alignment",
    "04_cross_seed_summary",
    "05_program_ranking",
    "06_freeze",
    "07_logs",
]:

    (
        OUT_ROOT
        / subdirectory
    ).mkdir(
        parents=True,
        exist_ok=True,
    )


coverage_rows = []
mapping_rows = []
cell_score_frames = []
alignment_rows = []


###############################################################################
# Score frozen programs and align with latent subspaces
###############################################################################

for lineage, configuration in LINEAGES.items():

    print(
        "\n"
        + "=" * 78
    )

    print(
        "处理谱系：",
        lineage,
    )

    print(
        "=" * 78
    )

    h5ad_path = find_unique_file(
        CNMF_WORK,
        configuration[
            "h5ad_name"
        ],
    )


    adata = ad.read_h5ad(
        h5ad_path
    )


    required_obs_columns = {
        "donor_id",
        "condition",
    }


    missing_obs_columns = (
        required_obs_columns
        - set(
            adata.obs.columns
        )
    )


    if missing_obs_columns:

        raise RuntimeError(
            f"{lineage}缺少元数据列："
            + ", ".join(
                sorted(
                    missing_obs_columns
                )
            )
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


    normalized_var_names = [
        normalize_gene(
            value
        )
        for value in adata.var_names
    ]


    gene_to_index = {}


    for index, gene in enumerate(
        normalized_var_names
    ):

        if (
            gene
            and gene not in gene_to_index
        ):

            gene_to_index[
                gene
            ] = index


    lineage_programs = eligible_programs.loc[
        eligible_programs[
            "analysis_lineage"
        ] == lineage
    ].copy()


    if lineage_programs.empty:

        raise RuntimeError(
            f"{lineage}没有内容QC可解释程序。"
        )


    score_lookup = {}


    cell_score_table = pd.DataFrame(
        {
            "cell_id":
                adata.obs_names.astype(
                    str
                ),

            "donor_id":
                adata.obs[
                    "donor_id"
                ].to_numpy(),

            "condition":
                adata.obs[
                    "condition_model"
                ].to_numpy(),

            "lineage":
                lineage,
        }
    )


    for program_row in lineage_programs.itertuples(
        index=False
    ):

        program_id = str(
            program_row.program_id
        )


        spectrum_path = one_file(
            B3A_ROOT,
            f"{program_id}__meta_spectrum.csv",
        )


        spectrum = read_program_spectrum(
            spectrum_path
        )


        for (
            score_definition,
            top_n,
        ) in SCORE_DEFINITIONS.items():

            (
                program_score,
                detected_mapping,
            ) = calculate_program_score(
                adata=adata,
                normalized_var_names=
                    normalized_var_names,
                gene_to_index=
                    gene_to_index,
                spectrum=spectrum,
                top_n=top_n,
            )


            score_key = (
                program_id,
                score_definition,
            )


            score_lookup[
                score_key
            ] = program_score


            output_column = (
                f"{program_id}__"
                f"{score_definition}"
            )


            cell_score_table[
                output_column
            ] = program_score


            coverage_rows.append(
                {
                    "lineage":
                        lineage,

                    "program_id":
                        program_id,

                    "score_definition":
                        score_definition,

                    "top_n":
                        top_n,

                    "frozen_top_genes":
                        min(
                            top_n,
                            len(
                                spectrum
                            ),
                        ),

                    "detected_variable_genes":
                        len(
                            detected_mapping
                        ),

                    "coverage_fraction":
                        (
                            len(
                                detected_mapping
                            )
                            /
                            min(
                                top_n,
                                len(
                                    spectrum
                                ),
                            )
                        ),

                    "spectrum_path":
                        str(
                            spectrum_path
                        ),

                    "spectrum_sha256":
                        sha256_file(
                            spectrum_path
                        ),
                }
            )


            detected_mapping = (
                detected_mapping
                .copy()
            )


            detected_mapping.insert(
                0,
                "lineage",
                lineage,
            )


            detected_mapping.insert(
                1,
                "program_id",
                program_id,
            )


            detected_mapping.insert(
                2,
                "score_definition",
                score_definition,
            )


            mapping_rows.append(
                detected_mapping
            )


    cell_score_frames.append(
        cell_score_table
    )


    h5ad_index = pd.Index(
        adata.obs_names.astype(
            str
        )
    )


    for seed in SEEDS:

        latent_path = one_file(
            B1_ROOT,
            (
                f"34A1B1_{lineage}_"
                f"seed_{seed}_latents.npz"
            ),
        )


        with np.load(
            latent_path,
            allow_pickle=True,
        ) as latent_file:

            required_arrays = {
                "cell_id",
                "donor_id",
                "condition",
                "salient",
                "background",
            }

            missing_arrays = (
                required_arrays
                - set(
                    latent_file.files
                )
            )

            if missing_arrays:

                raise RuntimeError(
                    f"{lineage} seed={seed}潜变量文件缺少数组："
                    + ", ".join(
                        sorted(
                            missing_arrays
                        )
                    )
                )

            latent_cell_ids = np.asarray(
                latent_file[
                    "cell_id"
                ]
            ).astype(
                str
            )

            latent_donor = np.asarray(
                latent_file[
                    "donor_id"
                ]
            ).astype(
                str
            )

            latent_condition = np.asarray(
                latent_file[
                    "condition"
                ]
            ).astype(
                str
            )

            salient = np.asarray(
                latent_file[
                    "salient"
                ],
                dtype=float,
            )

            background = np.asarray(
                latent_file[
                    "background"
                ],
                dtype=float,
            )


        positions = h5ad_index.get_indexer(
            latent_cell_ids
        )


        if np.any(
            positions < 0
        ):

            raise RuntimeError(
                f"{lineage} seed={seed} latent"
                "细胞无法与H5AD匹配。"
            )


        if len(
            np.unique(
                positions
            )
        ) != adata.n_obs:

            raise RuntimeError(
                f"{lineage} seed={seed}不是"
                "H5AD细胞的一一对应。"
            )


        h5ad_donor_reordered = (
            adata.obs[
                "donor_id"
            ]
            .to_numpy()[
                positions
            ]
            .astype(
                str
            )
        )


        h5ad_condition_reordered = (
            adata.obs[
                "condition_model"
            ]
            .to_numpy()[
                positions
            ]
            .astype(
                str
            )
        )


        if not np.array_equal(
            latent_donor,
            h5ad_donor_reordered,
        ):

            raise RuntimeError(
                f"{lineage} seed={seed}"
                " donor_id不一致。"
            )


        if not np.array_equal(
            latent_condition,
            h5ad_condition_reordered,
        ):

            raise RuntimeError(
                f"{lineage} seed={seed}"
                " condition不一致。"
            )


        for (
            program_id,
            score_definition,
        ), program_score in score_lookup.items():

            reordered_score = (
                program_score[
                    positions
                ]
            )


            program_metadata = (
                lineage_programs.loc[
                    lineage_programs[
                        "program_id"
                    ] == program_id
                ]
            )


            if len(
                program_metadata
            ) != 1:

                raise RuntimeError(
                    "程序元数据无法唯一定位："
                    f"{program_id}"
                )


            program_metadata = (
                program_metadata.iloc[
                    0
                ]
            )


            for mode in CENTERING_MODES:

                (
                    salient_input,
                    score_input,
                ) = prepare_for_mode(
                    latent=salient,
                    program_score=
                        reordered_score,
                    donor_id=
                        latent_donor,
                    condition=
                        latent_condition,
                    mode=mode,
                )


                (
                    background_input,
                    background_score_input,
                ) = prepare_for_mode(
                    latent=background,
                    program_score=
                        reordered_score,
                    donor_id=
                        latent_donor,
                    condition=
                        latent_condition,
                    mode=mode,
                )


                if not np.allclose(
                    score_input,
                    background_score_input,
                    atol=1e-12,
                ):

                    raise RuntimeError(
                        "salient/background使用的"
                        "程序分数不一致。"
                    )


                salient_cka = linear_cka(
                    salient_input,
                    score_input,
                )


                background_cka = linear_cka(
                    background_input,
                    score_input,
                )


                ratio = (
                    salient_cka
                    / background_cka

                    if (
                        np.isfinite(
                            background_cka
                        )
                        and background_cka
                        > 1e-12
                    )

                    else np.nan
                )


                alignment_rows.append(
                    {
                        "lineage":
                            lineage,

                        "program_id":
                            program_id,

                        "score_definition":
                            score_definition,

                        "seed":
                            seed,

                        "centering_mode":
                            mode,

                        "observations_used":
                            int(
                                len(
                                    score_input
                                )
                            ),

                        "salient_linear_CKA":
                            salient_cka,

                        "background_linear_CKA":
                            background_cka,

                        "salient_minus_background_CKA":
                            (
                                salient_cka
                                - background_cka
                            ),

                        "salient_to_background_ratio":
                            ratio,

                        "pre_registered_target":
                            bool(
                                program_metadata[
                                    "pre_registered_target"
                                ]
                            ),

                        "target_role":
                            (
                                ""
                                if pd.isna(
                                    program_metadata.get(
                                        "target_role",
                                        np.nan,
                                    )
                                )
                                else str(
                                    program_metadata.get(
                                        "target_role",
                                        "",
                                    )
                                )
                            ),

                        "technical_stability_class":
                            str(
                                stability.loc[
                                    stability[
                                        "lineage"
                                    ] == lineage,
                                    "technical_stability_class",
                                ].iloc[
                                    0
                                ]
                            ),

                        "latent_path":
                            str(
                                latent_path
                            ),
                    }
                )


    del adata


###############################################################################
# Assemble outputs
###############################################################################

coverage_table = pd.DataFrame(
    coverage_rows
)


mapping_table = pd.concat(
    mapping_rows,
    ignore_index=True,
)


cell_score_table = pd.concat(
    cell_score_frames,
    ignore_index=True,
)


alignment_table = pd.DataFrame(
    alignment_rows
)


expected_alignment_rows = (
    len(
        eligible_programs
    )
    *
    len(
        SCORE_DEFINITIONS
    )
    *
    len(
        SEEDS
    )
    *
    len(
        CENTERING_MODES
    )
)


if len(
    alignment_table
) != expected_alignment_rows:

    raise RuntimeError(
        "对齐结果行数异常："
        f"expected={expected_alignment_rows}, "
        f"actual={len(alignment_table)}"
    )


summary_table = (
    alignment_table
    .groupby(
        [
            "lineage",
            "program_id",
            "score_definition",
            "centering_mode",
            "pre_registered_target",
            "target_role",
            "technical_stability_class",
        ],
        as_index=False,
        dropna=False,
    )
    .agg(
        seeds_completed=(
            "seed",
            "nunique",
        ),

        median_salient_CKA=(
            "salient_linear_CKA",
            "median",
        ),

        minimum_salient_CKA=(
            "salient_linear_CKA",
            "min",
        ),

        maximum_salient_CKA=(
            "salient_linear_CKA",
            "max",
        ),

        median_background_CKA=(
            "background_linear_CKA",
            "median",
        ),

        median_salient_minus_background_CKA=(
            "salient_minus_background_CKA",
            "median",
        ),

        minimum_salient_minus_background_CKA=(
            "salient_minus_background_CKA",
            "min",
        ),

        maximum_salient_minus_background_CKA=(
            "salient_minus_background_CKA",
            "max",
        ),

        median_salient_to_background_ratio=(
            "salient_to_background_ratio",
            "median",
        ),

        seeds_with_salient_above_background=(
            "salient_minus_background_CKA",
            lambda values:
                int(
                    (
                        np.asarray(
                            values
                        )
                        > 0
                    ).sum()
                ),
        ),
    )
)


primary_ranking = summary_table.loc[
    (
        summary_table[
            "score_definition"
        ]
        == "TOP100_WEIGHTED"
    )
    &
    (
        summary_table[
            "centering_mode"
        ]
        ==
        "DONOR_CENTROID_CONDITION_CENTERED"
    )
].copy()


primary_ranking[
    "rank_by_salient_CKA_within_lineage"
] = (
    primary_ranking
    .groupby(
        "lineage"
    )[
        "median_salient_CKA"
    ]
    .rank(
        method="min",
        ascending=False,
    )
    .astype(
        int
    )
)


primary_ranking[
    "rank_by_salient_specific_delta_within_lineage"
] = (
    primary_ranking
    .groupby(
        "lineage"
    )[
        "median_salient_minus_background_CKA"
    ]
    .rank(
        method="min",
        ascending=False,
    )
    .astype(
        int
    )
)


supportive_ranking = summary_table.loc[
    (
        summary_table[
            "score_definition"
        ]
        == "TOP100_WEIGHTED"
    )
    &
    (
        summary_table[
            "centering_mode"
        ]
        ==
        "CELL_DONOR_CENTERED"
    )
].copy()


supportive_ranking[
    "rank_by_salient_CKA_within_lineage"
] = (
    supportive_ranking
    .groupby(
        "lineage"
    )[
        "median_salient_CKA"
    ]
    .rank(
        method="min",
        ascending=False,
    )
    .astype(
        int
    )
)


supportive_ranking[
    "rank_by_salient_specific_delta_within_lineage"
] = (
    supportive_ranking
    .groupby(
        "lineage"
    )[
        "median_salient_minus_background_CKA"
    ]
    .rank(
        method="min",
        ascending=False,
    )
    .astype(
        int
    )
)


###############################################################################
# Save result tables
###############################################################################

coverage_table.to_csv(
    OUT_ROOT
    / "01_program_coverage"
    / "34A1B2_PROGRAM_SCORE_COVERAGE.csv",
    index=False,
)


mapping_table.to_csv(
    OUT_ROOT
    / "01_program_coverage"
    / "34A1B2_PROGRAM_SCORE_GENES_AND_WEIGHTS.csv",
    index=False,
)


cell_score_table.to_csv(
    OUT_ROOT
    / "02_cell_program_scores"
    / "34A1B2_CELL_FROZEN_PROGRAM_SCORES.csv",
    index=False,
)


alignment_table.to_csv(
    OUT_ROOT
    / "03_seed_level_alignment"
    / "34A1B2_SEED_LEVEL_PROGRAM_ALIGNMENT.csv",
    index=False,
)


summary_table.to_csv(
    OUT_ROOT
    / "04_cross_seed_summary"
    / "34A1B2_CROSS_SEED_ALIGNMENT_SUMMARY.csv",
    index=False,
)


primary_ranking.sort_values(
    [
        "lineage",
        "rank_by_salient_specific_delta_within_lineage",
    ]
).to_csv(
    OUT_ROOT
    / "05_program_ranking"
    / "34A1B2_PRIMARY_DONOR_LEVEL_PROGRAM_RANKING.csv",
    index=False,
)


supportive_ranking.sort_values(
    [
        "lineage",
        "rank_by_salient_specific_delta_within_lineage",
    ]
).to_csv(
    OUT_ROOT
    / "05_program_ranking"
    / "34A1B2_SUPPORTIVE_WITHIN_DONOR_PROGRAM_RANKING.csv",
    index=False,
)


###############################################################################
# Freeze stage
###############################################################################

STATUS = (
    "READY_FOR_34A1B3_"
    "PROGRAM_ALIGNMENT_REVIEW_AND_LODO_TARGET_SELECTION"
)


lock = {
    "stage":
        "34A1B2",

    "created_at":
        datetime.now().isoformat(),

    "status":
        STATUS,

    "programs_aligned":
        int(
            len(
                eligible_programs
            )
        ),

    "score_definitions":
        list(
            SCORE_DEFINITIONS
        ),

    "seeds":
        SEEDS,

    "centering_modes":
        CENTERING_MODES,

    "primary_alignment_metric":
        (
            "TOP100_WEIGHTED with "
            "DONOR_CENTROID_CONDITION_CENTERED linear CKA"
        ),

    "supportive_alignment_metric":
        (
            "TOP100_WEIGHTED with "
            "CELL_DONOR_CENTERED linear CKA"
        ),

    "individual_latent_dimensions_interpreted":
        False,

    "cell_level_disease_p_values_calculated":
        False,

    "alignment_metrics_are_descriptive":
        True,

    "final_biological_state_classification_performed":
        False,

    "program_refitted":
        False,

    "program_gene_reselection_performed":
        False,

    "donor_is_inferential_unit":
        True,

    "senescent_cell_label_authorized":
        False,

    "causal_claim_authorized":
        False,

    "next_stage":
        (
            "34A1B3_PROGRAM_ALIGNMENT_REVIEW_"
            "AND_LODO_TARGET_SELECTION"
        ),
}


(
    OUT_ROOT
    / "06_freeze"
    / "34A1B2_PROGRAM_ALIGNMENT_LOCK.json"
).write_text(
    json.dumps(
        lock,
        ensure_ascii=False,
        indent=2,
    ),
    encoding="utf-8",
)


###############################################################################
# Summary
###############################################################################

summary_lines = [
    "# 34A-1B2 Frozen-program Alignment",
    "",
    f"## Status: {STATUS}",
    "",
    "## Design",
    "",
    (
        f"- Content-QC-eligible frozen programs: "
        f"{len(eligible_programs)}."
    ),
    (
        "- Program scores were calculated from frozen "
        "top-50 and top-100 meta-spectrum genes."
    ),
    (
        "- Alignment was assessed against the complete "
        "salient and background latent subspaces."
    ),
    (
        "- Three contrastiveVI seeds were summarized."
    ),
    (
        "- Individual latent axes were not interpreted."
    ),
    "",
    "## Primary review table",
    "",
]


for row in primary_ranking.sort_values(
    [
        "lineage",
        "rank_by_salient_specific_delta_within_lineage",
    ]
).itertuples(
    index=False
):

    summary_lines.append(
        (
            f"- {row.lineage} / {row.program_id}: "
            f"salient CKA={row.median_salient_CKA:.4f}; "
            f"background CKA={row.median_background_CKA:.4f}; "
            f"delta="
            f"{row.median_salient_minus_background_CKA:.4f}; "
            f"salient>background seeds="
            f"{row.seeds_with_salient_above_background}/3; "
            f"pre-registered="
            f"{bool(row.pre_registered_target)}."
        )
    )


summary_lines.extend(
    [
        "",
        "## Interpretation boundary",
        "",
        (
            "- Alignment does not establish that a program "
            "is TLE-associated; disease effects were tested "
            "in the upstream donor-level analysis."
        ),
        (
            "- Alignment does not prove cellular senescence, "
            "inflammation, metabolism or causality."
        ),
        (
            "- Donor-centroid analyses contain only eight "
            "biological units and are descriptive."
        ),
        (
            "- Final target selection requires joint review "
            "of donor-level disease effect, content QC, "
            "seed stability and program alignment."
        ),
        "",
        "## Next stage",
        "",
        lock[
            "next_stage"
        ],
    ]
)


(
    OUT_ROOT
    / "34A1B2_EXECUTIVE_SUMMARY.md"
).write_text(
    "\n".join(
        summary_lines
    ),
    encoding="utf-8",
)


###############################################################################
# Provenance and manifest
###############################################################################

provenance = {
    "stage":
        "34A1B2",

    "created_at":
        datetime.now().isoformat(),

    "34A0_release":
        str(
            A0_ROOT
        ),

    "32B3A_release":
        str(
            B3A_ROOT
        ),

    "34A1B1_release":
        str(
            B1_ROOT
        ),

    "master_evidence_matrix_sha256":
        sha256_file(
            MASTER_PATH
        ),

    "target_registry_sha256":
        sha256_file(
            TARGET_PATH
        ),

    "stability_table_sha256":
        sha256_file(
            STABILITY_PATH
        ),

    "new_contrastivevi_model_fitted":
        False,

    "new_disease_test_performed":
        False,

    "program_refitted":
        False,
}


(
    OUT_ROOT
    / "07_logs"
    / "34A1B2_INPUT_PROVENANCE.json"
).write_text(
    json.dumps(
        provenance,
        ensure_ascii=False,
        indent=2,
    ),
    encoding="utf-8",
)


(
    OUT_ROOT
    / "07_logs"
    / "34A1B2_python_session.txt"
).write_text(
    (
        f"python={sys.version}\n"
        f"numpy={np.__version__}\n"
        f"pandas={pd.__version__}\n"
        f"anndata={ad.__version__}\n"
    ),
    encoding="utf-8",
)


output_files = sorted(
    path
    for path in OUT_ROOT.rglob(
        "*"
    )
    if path.is_file()
)


manifest = pd.DataFrame(
    [
        {
            "relative_path":
                path.relative_to(
                    OUT_ROOT
                ).as_posix(),

            "size_bytes":
                path.stat().st_size,

            "sha256":
                sha256_file(
                    path
                ),
        }
        for path in output_files
    ]
)


manifest.to_csv(
    OUT_ROOT
    / "34A1B2_FILE_MANIFEST.csv",
    index=False,
)


ZIP_PATH = Path(
    str(
        OUT_ROOT
    )
    + ".zip"
)


with zipfile.ZipFile(
    ZIP_PATH,
    mode="w",
    compression=zipfile.ZIP_DEFLATED,
    compresslevel=9,
) as archive:

    for path in sorted(
        OUT_ROOT.rglob(
            "*"
        )
    ):

        if path.is_file():

            archive.write(
                path,
                arcname=path.relative_to(
                    OUT_ROOT
                ).as_posix(),
            )


ZIP_HASH = sha256_file(
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


###############################################################################
# Console output
###############################################################################

print(
    "\n"
    + "=" * 78
)


print(
    "34A-1B2 completed"
)


print(
    "状态：",
    STATUS,
)


print(
    "对齐程序数：",
    len(
        eligible_programs
    ),
)


print(
    "对齐结果行数：",
    len(
        alignment_table
    ),
)


print(
    "结果目录：",
    OUT_ROOT,
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