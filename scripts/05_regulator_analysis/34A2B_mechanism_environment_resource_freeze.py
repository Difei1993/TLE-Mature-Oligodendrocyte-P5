from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import zipfile
from datetime import datetime
from pathlib import Path

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd


PROJECT = Path(
    os.environ.get(
        "TLE_PROJECT",
        __import__("os").environ["P5_WORK_ROOT"],
    )
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

OUT = (
    RELEASE_ROOT
    / f"34A2B_MECHANISM_ENVIRONMENT_RESOURCE_FREEZE_{STAMP}"
)

EXPECTED_DECOUPLER = "2.2.0"
TMIN = 5

H5AD_NAMES = {
    "Astrocytes":
        "Astrocytes_primary.h5ad",

    "Microglia_macrophage":
        "Microglia_macrophage_primary.h5ad",

    "Oligodendrocytes":
        "Oligodendrocytes_primary.h5ad",
}


def sha256(
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
            for filename in required
        )

        if complete:
            return root.resolve()

    raise FileNotFoundError(
        f"没有找到完整release：{prefix}"
    )


def one(
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

    if len(
        hits
    ) != 1:

        raise RuntimeError(
            f"文件未能唯一定位：{filename}\n"
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
    filename: str,
) -> Path:

    hits = [
        path
        for path in CNMF_WORK.rglob(
            filename
        )
        if path.is_file()
    ]

    if len(
        hits
    ) != 1:

        raise RuntimeError(
            f"{filename}找到{len(hits)}个：\n"
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


def norm_gene(
    value,
) -> str:

    value = (
        ""
        if value is None
        else str(
            value
        ).strip().upper()
    )

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


def fetch_collectri(
    remove_complexes: bool,
) -> pd.DataFrame:

    last_error = None

    for attempt in range(
        1,
        4,
    ):

        try:

            result = dc.op.collectri(
                organism="human",
                remove_complexes=
                    remove_complexes,
                license="academic",
                verbose=True,
            )

            if (
                not isinstance(
                    result,
                    pd.DataFrame,
                )
                or result.empty
            ):

                raise RuntimeError(
                    "CollecTRI返回空表。"
                )

            return result.copy()

        except Exception as exc:

            last_error = exc

            print(
                (
                    f"下载失败，第{attempt}/3次："
                    f"{type(exc).__name__}: {exc}"
                ),
                flush=True,
            )

            if attempt < 3:
                time.sleep(
                    10
                )

    raise RuntimeError(
        f"CollecTRI下载失败：{last_error}"
    )


def prepare_network(
    raw: pd.DataFrame,
) -> tuple[
    pd.DataFrame,
    dict,
]:

    required = {
        "source",
        "target",
        "weight",
    }

    missing = (
        required
        - set(
            raw.columns
        )
    )

    if missing:

        raise RuntimeError(
            "CollecTRI缺少必要列："
            f"{sorted(missing)}"
        )

    network = raw.copy()

    network[
        "source"
    ] = network[
        "source"
    ].astype(
        str
    ).str.strip()

    network[
        "target"
    ] = network[
        "target"
    ].astype(
        str
    ).str.strip()

    network[
        "target_norm"
    ] = network[
        "target"
    ].map(
        norm_gene
    )

    network[
        "weight"
    ] = pd.to_numeric(
        network[
            "weight"
        ],
        errors="coerce",
    )

    network = network.loc[
        network[
            "source"
        ].ne("")
        &
        network[
            "target_norm"
        ].ne("")
        &
        np.isfinite(
            network[
                "weight"
            ]
        )
    ].copy()

    conflicts = (
        network
        .groupby(
            [
                "source",
                "target_norm",
            ],
            observed=True,
        )[
            "weight"
        ]
        .nunique()
        .gt(
            1
        )
    )

    if conflicts.any():

        raise RuntimeError(
            (
                "CollecTRI存在同一source-target"
                "的冲突权重："
                f"{int(conflicts.sum())}对。"
            )
        )

    duplicate_rows = int(
        network.duplicated(
            [
                "source",
                "target_norm",
            ],
            keep=False,
        ).sum()
    )

    network = (
        network
        .sort_values(
            [
                "source",
                "target_norm",
            ]
        )
        .drop_duplicates(
            [
                "source",
                "target_norm",
            ],
            keep="first",
        )
        .reset_index(
            drop=True
        )
    )

    summary = {
        "raw_rows":
            int(
                len(
                    raw
                )
            ),

        "valid_unique_interactions":
            int(
                len(
                    network
                )
            ),

        "unique_sources":
            int(
                network[
                    "source"
                ].nunique()
            ),

        "unique_targets":
            int(
                network[
                    "target_norm"
                ].nunique()
            ),

        "duplicate_rows_before_collapse":
            duplicate_rows,
    }

    return (
        network,
        summary,
    )


###############################################################################
# Locate frozen inputs
###############################################################################


A2 = latest(
    "34A2_DONOR_LEVEL_PROGRAM_MECHANISM_INPUT_AUDIT_",
    [
        "34A2_MECHANISM_INPUT_AUDIT_LOCK.json",
    ],
)


B5 = latest(
    "34A1B5_FINAL_CONTRASTIVEVI_CONSENSUS_",
    [
        "34A1B5_FINAL_PROGRAM_CONSENSUS.csv",
    ],
)


B2 = latest(
    "34A1B2_FROZEN_PROGRAM_SALIENT_ALIGNMENT_",
    [
        "34A1B2_PROGRAM_SCORE_GENES_AND_WEIGHTS.csv",
    ],
)


A2_LOCK_PATH = one(
    A2,
    "34A2_MECHANISM_INPUT_AUDIT_LOCK.json",
)


CONSENSUS_PATH = one(
    B5,
    "34A1B5_FINAL_PROGRAM_CONSENSUS.csv",
)


PROGRAM_GENES_PATH = one(
    B2,
    "34A1B2_PROGRAM_SCORE_GENES_AND_WEIGHTS.csv",
)


a2_lock = json.loads(
    A2_LOCK_PATH.read_text(
        encoding="utf-8"
    )
)


if (
    a2_lock.get(
        "status"
    )
    !=
    "INPUTS_PASS_MECHANISM_ENVIRONMENT_PREPARATION_REQUIRED"
):

    raise RuntimeError(
        "34A2状态不正确："
        f"{a2_lock.get('status')}"
    )


if not a2_lock.get(
    "input_pass",
    False,
):

    raise RuntimeError(
        "34A2生物学输入没有通过。"
    )


if (
    a2_lock.get(
        "next_stage"
    )
    !=
    "34A2B_MECHANISM_ENVIRONMENT_AND_RESOURCE_FREEZE"
):

    raise RuntimeError(
        "34A2未授权进入34A2B。"
    )


if str(
    dc.__version__
) != EXPECTED_DECOUPLER:

    raise RuntimeError(
        (
            "decoupler版本错误："
            f"expected={EXPECTED_DECOUPLER}, "
            f"actual={dc.__version__}"
        )
    )


consensus = pd.read_csv(
    CONSENSUS_PATH
)


program_genes = pd.read_csv(
    PROGRAM_GENES_PATH
)


for column in [
    "lineage",
    "program_id",
    "final_tier",
    "final_role",
]:

    if column not in consensus.columns:

        raise RuntimeError(
            f"共识表缺少列：{column}"
        )


for column in [
    "lineage",
    "program_id",
    "score_definition",
    "gene",
]:

    if column not in program_genes.columns:

        raise RuntimeError(
            f"程序基因表缺少列：{column}"
        )


consensus[
    "lineage"
] = consensus[
    "lineage"
].astype(
    str
)


consensus[
    "program_id"
] = consensus[
    "program_id"
].astype(
    str
)


top100 = program_genes.loc[
    program_genes[
        "score_definition"
    ].astype(
        str
    ).eq(
        "TOP100_WEIGHTED"
    )
].copy()


top100[
    "gene_norm"
] = top100[
    "gene"
].map(
    norm_gene
)


top100 = top100.loc[
    top100[
        "gene_norm"
    ].ne("")
].copy()


for program_id in consensus[
    "program_id"
]:

    n_genes = top100.loc[
        top100[
            "program_id"
        ].astype(
            str
        ).eq(
            program_id
        ),
        "gene_norm",
    ].nunique()

    if n_genes != 100:

        raise RuntimeError(
            (
                f"{program_id}的Top100"
                f"基因数为{n_genes}。"
            )
        )


###############################################################################
# Output directories and environment
###############################################################################


OUT.mkdir(
    parents=True,
    exist_ok=False,
)


for folder in [
    "01_environment",
    "02_raw_resources",
    "03_core_resources",
    "04_lineage_networks",
    "05_program_orthogonal_networks",
    "06_audit",
    "07_freeze",
    "08_logs",
]:

    (
        OUT
        / folder
    ).mkdir(
        parents=True,
        exist_ok=True,
    )


pip_freeze = subprocess.check_output(
    [
        sys.executable,
        "-m",
        "pip",
        "freeze",
    ],
    text=True,
)


(
    OUT
    / "01_environment"
    / "34A2B_PIP_FREEZE.txt"
).write_text(
    pip_freeze,
    encoding="utf-8",
)


environment = {
    "python":
        sys.version,

    "python_executable":
        sys.executable,

    "platform":
        platform.platform(),

    "decoupler":
        str(
            dc.__version__
        ),

    "anndata":
        str(
            ad.__version__
        ),

    "numpy":
        str(
            np.__version__
        ),

    "pandas":
        str(
            pd.__version__
        ),
}


(
    OUT
    / "01_environment"
    / "34A2B_ENVIRONMENT.json"
).write_text(
    json.dumps(
        environment,
        ensure_ascii=False,
        indent=2,
    ),
    encoding="utf-8",
)


###############################################################################
# Download and freeze CollecTRI
###############################################################################


print(
    "下载CollecTRI（保留复合物）……",
    flush=True,
)


raw_with = fetch_collectri(
    remove_complexes=False
)


print(
    "下载CollecTRI（移除复合物）……",
    flush=True,
)


raw_without = fetch_collectri(
    remove_complexes=True
)


RAW_WITH_PATH = (
    OUT
    / "02_raw_resources"
    / "34A2B_COLLECTRI_HUMAN_WITH_COMPLEXES_RAW.csv"
)


RAW_WITHOUT_PATH = (
    OUT
    / "02_raw_resources"
    / "34A2B_COLLECTRI_HUMAN_NO_COMPLEXES_RAW.csv"
)


raw_with.to_csv(
    RAW_WITH_PATH,
    index=False,
)


raw_without.to_csv(
    RAW_WITHOUT_PATH,
    index=False,
)


(
    core_with,
    summary_with,
) = prepare_network(
    raw_with
)


(
    core_without,
    summary_without,
) = prepare_network(
    raw_without
)


CORE_WITH_PATH = (
    OUT
    / "03_core_resources"
    / "34A2B_COLLECTRI_HUMAN_WITH_COMPLEXES_CORE.csv"
)


CORE_WITHOUT_PATH = (
    OUT
    / "03_core_resources"
    / "34A2B_COLLECTRI_HUMAN_NO_COMPLEXES_CORE.csv"
)


core_with[
    [
        "source",
        "target_norm",
        "weight",
    ]
].rename(
    columns={
        "target_norm":
            "target",
    }
).to_csv(
    CORE_WITH_PATH,
    index=False,
)


core_without[
    [
        "source",
        "target_norm",
        "weight",
    ]
].rename(
    columns={
        "target_norm":
            "target",
    }
).to_csv(
    CORE_WITHOUT_PATH,
    index=False,
)


###############################################################################
# Lineage coverage and program-orthogonal networks
###############################################################################


lineage_rows = []
program_rows = []


for (
    lineage,
    h5ad_name,
) in H5AD_NAMES.items():

    print(
        f"审计谱系：{lineage}",
        flush=True,
    )

    h5ad_path = find_h5ad(
        h5ad_name
    )

    adata = ad.read_h5ad(
        h5ad_path,
        backed="r",
    )

    gene_universe = {
        norm_gene(
            value
        )
        for value in adata.var_names
    }

    gene_universe.discard(
        ""
    )

    lineage_network = core_with.loc[
        core_with[
            "target_norm"
        ].isin(
            gene_universe
        )
    ].copy()

    target_counts = (
        lineage_network
        .groupby(
            "source",
            observed=True,
        )[
            "target_norm"
        ]
        .nunique()
    )

    eligible_sources = set(
        target_counts.loc[
            target_counts
            >= TMIN
        ].index.astype(
            str
        )
    )

    lineage_network = lineage_network.loc[
        lineage_network[
            "source"
        ].astype(
            str
        ).isin(
            eligible_sources
        )
    ].copy()

    lineage_path = (
        OUT
        / "04_lineage_networks"
        / f"34A2B_{lineage}_COLLECTRI_TMIN5.csv"
    )

    lineage_network[
        [
            "source",
            "target_norm",
            "weight",
        ]
    ].rename(
        columns={
            "target_norm":
                "target",
        }
    ).to_csv(
        lineage_path,
        index=False,
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
                sha256(
                    h5ad_path
                ),

            "genes_in_h5ad":
                len(
                    gene_universe
                ),

            "interactions_tmin5":
                len(
                    lineage_network
                ),

            "sources_tmin5":
                lineage_network[
                    "source"
                ].nunique(),

            "targets":
                lineage_network[
                    "target_norm"
                ].nunique(),

            "resource_path":
                str(
                    lineage_path
                ),

            "resource_sha256":
                sha256(
                    lineage_path
                ),
        }
    )

    lineage_programs = consensus.loc[
        consensus[
            "lineage"
        ].eq(
            lineage
        )
    ]

    for row in lineage_programs.itertuples(
        index=False
    ):

        program_id = str(
            row.program_id
        )

        program_gene_set = set(
            top100.loc[
                top100[
                    "program_id"
                ].astype(
                    str
                ).eq(
                    program_id
                ),
                "gene_norm",
            ].astype(
                str
            )
        )

        overlap = sorted(
            program_gene_set
            &
            set(
                lineage_network[
                    "target_norm"
                ].astype(
                    str
                )
            )
        )

        orthogonal = lineage_network.loc[
            ~lineage_network[
                "target_norm"
            ].isin(
                program_gene_set
            )
        ].copy()

        orthogonal_counts = (
            orthogonal
            .groupby(
                "source",
                observed=True,
            )[
                "target_norm"
            ]
            .nunique()
        )

        orthogonal_sources = set(
            orthogonal_counts.loc[
                orthogonal_counts
                >= TMIN
            ].index.astype(
                str
            )
        )

        orthogonal = orthogonal.loc[
            orthogonal[
                "source"
            ].astype(
                str
            ).isin(
                orthogonal_sources
            )
        ].copy()

        orthogonal_path = (
            OUT
            / "05_program_orthogonal_networks"
            / (
                f"34A2B_{program_id}_"
                "EXCLUDE_TOP100_TMIN5.csv"
            )
        )

        orthogonal[
            [
                "source",
                "target_norm",
                "weight",
            ]
        ].rename(
            columns={
                "target_norm":
                    "target",
            }
        ).to_csv(
            orthogonal_path,
            index=False,
        )

        program_rows.append(
            {
                "lineage":
                    lineage,

                "program_id":
                    program_id,

                "final_tier":
                    str(
                        row.final_tier
                    ),

                "final_role":
                    str(
                        row.final_role
                    ),

                "top100_collectri_target_overlap":
                    len(
                        overlap
                    ),

                "top100_overlap_fraction":
                    len(
                        overlap
                    )
                    / 100.0,

                "baseline_interactions_tmin5":
                    len(
                        lineage_network
                    ),

                "baseline_sources_tmin5":
                    lineage_network[
                        "source"
                    ].nunique(),

                "orthogonal_interactions_tmin5":
                    len(
                        orthogonal
                    ),

                "orthogonal_sources_tmin5":
                    orthogonal[
                        "source"
                    ].nunique(),

                "sources_lost_after_exclusion":
                    (
                        lineage_network[
                            "source"
                        ].nunique()
                        -
                        orthogonal[
                            "source"
                        ].nunique()
                    ),

                "overlap_genes":
                    ";".join(
                        overlap
                    ),

                "orthogonal_resource_path":
                    str(
                        orthogonal_path
                    ),

                "orthogonal_resource_sha256":
                    sha256(
                        orthogonal_path
                    ),
            }
        )

    try:
        adata.file.close()
    except Exception:
        pass

    del adata


lineage_audit = pd.DataFrame(
    lineage_rows
)


program_audit = pd.DataFrame(
    program_rows
)


resource_summary = pd.DataFrame(
    [
        {
            "resource":
                "CollecTRI_human_with_complexes",

            **summary_with,

            "raw_sha256":
                sha256(
                    RAW_WITH_PATH
                ),

            "core_sha256":
                sha256(
                    CORE_WITH_PATH
                ),
        },

        {
            "resource":
                "CollecTRI_human_no_complexes",

            **summary_without,

            "raw_sha256":
                sha256(
                    RAW_WITHOUT_PATH
                ),

            "core_sha256":
                sha256(
                    CORE_WITHOUT_PATH
                ),
        },
    ]
)


lineage_audit.to_csv(
    OUT
    / "06_audit"
    / "34A2B_LINEAGE_NETWORK_COVERAGE.csv",
    index=False,
)


program_audit.to_csv(
    OUT
    / "06_audit"
    / "34A2B_PROGRAM_ORTHOGONALIZATION_AUDIT.csv",
    index=False,
)


resource_summary.to_csv(
    OUT
    / "06_audit"
    / "34A2B_COLLECTRI_RESOURCE_SUMMARY.csv",
    index=False,
)


###############################################################################
# Final authorization
###############################################################################


primary = program_audit.loc[
    program_audit[
        "program_id"
    ].isin(
        {
            "Oligodendrocytes__P5",
            "Astrocytes__P6",
        }
    )
]


resource_pass = bool(
    len(
        primary
    ) == 2
    and primary[
        "orthogonal_sources_tmin5"
    ].ge(
        100
    ).all()
    and primary[
        "orthogonal_interactions_tmin5"
    ].ge(
        1000
    ).all()
    and lineage_audit[
        "sources_tmin5"
    ].ge(
        100
    ).all()
)


if resource_pass:

    STATUS = (
        "GO_FOR_34A3_DONOR_LEVEL_INTRINSIC_TF_ACTIVITY"
    )

    NEXT_STAGE = (
        "34A3_DONOR_LEVEL_INTRINSIC_TF_ACTIVITY"
    )

else:

    STATUS = (
        "HOLD_FOR_34A2B_RESOURCE_COVERAGE_REVIEW"
    )

    NEXT_STAGE = (
        "34A2B_RESOURCE_COVERAGE_REVIEW"
    )


lock = {
    "stage":
        "34A2B",

    "created_at":
        datetime.now().isoformat(),

    "status":
        STATUS,

    "source_34A2_release":
        str(
            A2
        ),

    "source_34A1B5_release":
        str(
            B5
        ),

    "source_34A1B2_release":
        str(
            B2
        ),

    "decoupler_version":
        str(
            dc.__version__
        ),

    "collectri_organism":
        "human",

    "collectri_license_argument":
        "academic",

    "primary_network_keeps_complexes":
        True,

    "no_complexes_sensitivity_frozen":
        True,

    "minimum_targets_per_tf":
        TMIN,

    "primary_tf_network":
        (
            "program-specific CollecTRI "
            "excluding its frozen Top100 genes"
        ),

    "supportive_tf_network":
        (
            "lineage-matched full CollecTRI network"
        ),

    "planned_normalization":
        (
            "library-size normalization "
            "followed by log1p"
        ),

    "inferential_unit":
        "donor",

    "planned_permutation":
        (
            "restricted permutations within "
            "Normal and TLE groups"
        ),

    "new_tf_activity_calculated":
        False,

    "new_disease_test_performed":
        False,

    "cell_level_p_values_authorized":
        False,

    "causal_claim_authorized":
        False,

    "liana_deferred":
        True,

    "resource_pass":
        resource_pass,

    "next_stage":
        NEXT_STAGE,
}


(
    OUT
    / "07_freeze"
    / "34A2B_RESOURCE_FREEZE_LOCK.json"
).write_text(
    json.dumps(
        lock,
        ensure_ascii=False,
        indent=2,
    ),
    encoding="utf-8",
)


summary_lines = [
    "# 34A2B mechanism environment and resource freeze",
    "",
    f"## Status: {STATUS}",
    "",
    f"- decoupler: {dc.__version__}",
    (
        "- Human CollecTRI with complexes "
        "is the primary network."
    ),
    (
        "- A no-complex network is frozen "
        "for sensitivity analysis."
    ),
    (
        "- Program-specific primary networks "
        "exclude each program's frozen Top100 "
        "genes to reduce circularity."
    ),
    "- Donor remains the inferential unit.",
    "- LIANA is deferred to the communication stage.",
    "",
    "## Next stage",
    "",
    NEXT_STAGE,
]


(
    OUT
    / "34A2B_EXECUTIVE_SUMMARY.md"
).write_text(
    "\n".join(
        summary_lines
    ),
    encoding="utf-8",
)


provenance = {
    "stage":
        "34A2B",

    "python":
        sys.version,

    "inputs":
        {
            str(
                A2_LOCK_PATH
            ):
                sha256(
                    A2_LOCK_PATH
                ),

            str(
                CONSENSUS_PATH
            ):
                sha256(
                    CONSENSUS_PATH
                ),

            str(
                PROGRAM_GENES_PATH
            ):
                sha256(
                    PROGRAM_GENES_PATH
                ),
        },

    "resources":
        {
            str(
                RAW_WITH_PATH
            ):
                sha256(
                    RAW_WITH_PATH
                ),

            str(
                RAW_WITHOUT_PATH
            ):
                sha256(
                    RAW_WITHOUT_PATH
                ),

            str(
                CORE_WITH_PATH
            ):
                sha256(
                    CORE_WITH_PATH
                ),

            str(
                CORE_WITHOUT_PATH
            ):
                sha256(
                    CORE_WITHOUT_PATH
                ),
        },
}


(
    OUT
    / "08_logs"
    / "34A2B_PROVENANCE.json"
).write_text(
    json.dumps(
        provenance,
        ensure_ascii=False,
        indent=2,
    ),
    encoding="utf-8",
)


###############################################################################
# Manifest and ZIP
###############################################################################


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
                sha256(
                    path
                ),
        }
        for path in output_files
    ]
).to_csv(
    OUT
    / "34A2B_FILE_MANIFEST.csv",
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


ZIP_HASH = sha256(
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
    "34A2B completed"
)

print(
    "状态：",
    STATUS,
)

print(
    "decoupler：",
    dc.__version__,
)

print(
    "谱系网络：",
    len(
        lineage_audit
    ),
    "/ 3",
)

print(
    "程序正交网络：",
    len(
        program_audit
    ),
    "/ 7",
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
