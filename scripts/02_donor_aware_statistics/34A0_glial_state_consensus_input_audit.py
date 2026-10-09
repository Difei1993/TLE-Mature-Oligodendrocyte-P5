from __future__ import annotations

import hashlib
import json
import sys
import zipfile
from datetime import datetime
from pathlib import Path

import pandas as pd

PROJECT = Path(
    __import__("os").environ["P5_WORK_ROOT"]
)

RELEASE_ROOT = PROJECT / "99_frozen_releases"

STAMP = datetime.now().strftime(
    "%Y%m%d_%H%M%S"
)

OUT = (
    RELEASE_ROOT
    / f"34A0_GLIAL_STATE_CONSENSUS_INPUT_AUDIT_{STAMP}"
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()

    with path.open("rb") as f:
        for block in iter(
            lambda: f.read(
                1024 * 1024
            ),
            b"",
        ):
            h.update(block)

    return h.hexdigest()


def latest_release(
    prefix: str,
    required: list[str],
) -> Path:

    candidates = sorted(
        [
            p
            for p in RELEASE_ROOT.glob(
                prefix + "*"
            )
            if p.is_dir()
        ],
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )

    for root in candidates:
        complete = all(
            len(
                [
                    p
                    for p in root.rglob(name)
                    if p.is_file()
                ]
            ) == 1
            for name in required
        )

        if complete:
            return root.resolve()

    raise FileNotFoundError(
        f"没有找到完整release：{prefix}\n"
        f"必要文件：{required}"
    )


def one_file(
    root: Path,
    name: str,
) -> Path:

    hits = [
        p
        for p in root.rglob(name)
        if p.is_file()
    ]

    if len(hits) != 1:
        raise RuntimeError(
            f"文件未能唯一定位：{name}\n"
            + "\n".join(
                map(
                    str,
                    hits,
                )
            )
        )

    return hits[0].resolve()


def bool_series(
    x: pd.Series,
) -> pd.Series:

    if pd.api.types.is_bool_dtype(x):
        return x.fillna(False)

    return (
        x.astype(str)
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


def check_unique(
    df: pd.DataFrame,
    name: str,
) -> None:

    if "program_id" not in df.columns:
        raise RuntimeError(
            f"{name}缺少program_id。"
        )

    if (
        df.program_id.isna().any()
        or df.program_id.duplicated().any()
    ):
        raise RuntimeError(
            f"{name}的program_id缺失或重复。"
        )


def text(
    row: pd.Series,
    column: str,
) -> str:

    value = row.get(
        column,
        "",
    )

    return (
        ""
        if pd.isna(value)
        else str(value).strip()
    )


def axis_tags(
    row: pd.Series,
) -> str:

    combined = " ".join(
        text(
            row,
            column,
        )
        for column in [
            "safe_label",
            "upstream_interpretation",
            "interpretation",
            "dominant_general_module",
            "content_flag_audit",
            "content_flag_interpretation",
            "formal_hit_names",
            "exploratory_hit_names",
            "top_general_pathways",
            "reviewed_anchor_pathways",
        ]
    ).lower()

    tags = []

    if any(
        term in combined
        for term in [
            "senescence",
            "stable_cell_cycle_exit",
            "stable cell-cycle",
            "fridman",
            "sasp",
            "ddr",
            "genome_instability",
        ]
    ):
        tags.append(
            "PARTIAL_SENESCENCE_ASSOCIATED"
        )

    if any(
        term in combined
        for term in [
            "inflamm",
            "ifn",
            "nfkb",
            "jak_stat",
            "jak-stat",
            "immediate",
            "stimulus-response",
            "stress",
            "oxidative",
            "chaperone",
            "immune-receptor",
        ]
    ):
        tags.append(
            "INFLAMMATION_OR_STRESS"
        )

    if any(
        term in combined
        for term in [
            "phagoc",
            "phagolys",
            "lipid",
            "lysosom",
            "autophagy",
            "myelin",
            "cholesterol",
            "metabol",
            "proteostasis",
            "iron",
        ]
    ):
        tags.append(
            "PHAGOCYTIC_LIPID_METABOLIC_REMODELING"
        )

    return ";".join(
        tags
        or [
            "LINEAGE_OR_HOMEOSTATIC_CONTEXT"
        ]
    )


def use_status(
    row: pd.Series,
) -> str:

    flag = text(
        row,
        "content_flag_audit",
    )

    review = text(
        row,
        "review_priority_audit",
    )

    status = text(
        row,
        "program_evidence_status",
    )

    tier = text(
        row,
        "evidence_tier",
    )

    interpretation = text(
        row,
        "interpretation",
    )

    program_id = text(
        row,
        "program_id",
    )

    if (
        flag
        == "LIKELY_LINEAGE_INCOMPATIBLE_RNA_PROGRAM"
    ):
        return (
            "EXCLUDE_LINEAGE_INCOMPATIBLE"
        )

    if (
        flag
        == "ACTIVE_CELL_CYCLE_COMPETING_STATE"
    ):
        return (
            "EXCLUDE_ACTIVE_PROLIFERATION"
        )

    if (
        "EXCLUDED" in status
        or "EXCLUDED_UNSTABLE" in review
    ):
        return "EXCLUDE_UNSTABLE"

    if (
        flag
        == "LOW_ANNOTATABILITY_LNCRNA_DOMINANT"
    ):
        return (
            "MANUAL_REVIEW_LOW_ANNOTATABILITY"
        )

    if (
        program_id
        == "Oligodendrocytes__P5"
    ):
        return (
            "PRIMARY_HIGH_PRIORITY_EXPLORATORY"
        )

    if (
        tier
        == "CONTEXTUAL_SENESCENCE_FEATURES"
        or interpretation
        == "SENESCENCE_FEATURES_WITHOUT_DISEASE_ASSOCIATION"
    ):
        return (
            "CONTEXTUAL_SENESCENCE_FEATURES_"
            "NO_DISEASE_ASSOCIATION"
        )

    if (
        tier
        == "COMPETING_RESPONSE"
        or interpretation
        == "COMPETING_RESPONSE_DOMINANT"
    ):
        return (
            "COMPETING_RESPONSE_CONTEXT"
        )

    if tier == "LINEAGE_CONTEXT":
        return "LINEAGE_CONTEXT"

    return (
        "SECONDARY_CONTEXT_OR_EXPLORATORY"
    )


###############################################################################
# 1. Locate frozen upstream releases
###############################################################################

B3A = latest_release(
    "32B3A_UNBIASED_PROGRAM_FREEZE_",
    [
        "32B3A_PROGRAM_TLE_EFFECT_SUMMARY.csv",
        "32B3A_FINAL_PROGRAM_LOCK.csv",
        "32B3A_ANALYSIS_LOCK.json",
    ],
)

B3B0 = latest_release(
    "32B3B0_PROGRAM_CONTENT_AUDIT_",
    [
        "32B3B0_PROGRAM_CONTENT_AUDIT.csv",
        "32B3B0_CONTENT_AUDIT_LOCK.json",
    ],
)

B3B1B = latest_release(
    "32B3B1B_PATHWAY_AND_26E_INTERPRETATION_",
    [
        "32B3B1B_PROGRAM_INTERPRETATION.csv",
        "32B3B1B_26E_OVERLAP_WITH_DISEASE_DIRECTION.csv",
        "32B3B1B_ANALYSIS_LOCK.json",
    ],
)

B3C = latest_release(
    "32B3C_MANUSCRIPT_SAFE_PROGRAM_FREEZE_",
    [
        "32B3C_MANUSCRIPT_SAFE_PROGRAM_LABELS.csv",
        "32B3C_CLAIM_MATRIX.csv",
        "32B3C_VALIDATION_PRIORITY.csv",
        "32B3C_MANUSCRIPT_SAFE_LOCK.json",
    ],
)

paths = {
    "effects":
        one_file(
            B3A,
            "32B3A_PROGRAM_TLE_EFFECT_SUMMARY.csv",
        ),

    "program_lock":
        one_file(
            B3A,
            "32B3A_FINAL_PROGRAM_LOCK.csv",
        ),

    "B3A_lock":
        one_file(
            B3A,
            "32B3A_ANALYSIS_LOCK.json",
        ),

    "content":
        one_file(
            B3B0,
            "32B3B0_PROGRAM_CONTENT_AUDIT.csv",
        ),

    "B3B0_lock":
        one_file(
            B3B0,
            "32B3B0_CONTENT_AUDIT_LOCK.json",
        ),

    "interpretation":
        one_file(
            B3B1B,
            "32B3B1B_PROGRAM_INTERPRETATION.csv",
        ),

    "overlap":
        one_file(
            B3B1B,
            "32B3B1B_26E_OVERLAP_WITH_DISEASE_DIRECTION.csv",
        ),

    "B3B1B_lock":
        one_file(
            B3B1B,
            "32B3B1B_ANALYSIS_LOCK.json",
        ),

    "labels":
        one_file(
            B3C,
            "32B3C_MANUSCRIPT_SAFE_PROGRAM_LABELS.csv",
        ),

    "claims":
        one_file(
            B3C,
            "32B3C_CLAIM_MATRIX.csv",
        ),

    "priority":
        one_file(
            B3C,
            "32B3C_VALIDATION_PRIORITY.csv",
        ),

    "B3C_lock":
        one_file(
            B3C,
            "32B3C_MANUSCRIPT_SAFE_LOCK.json",
        ),
}


###############################################################################
# 2. Read and audit
###############################################################################

effects = pd.read_csv(
    paths["effects"]
).rename(
    columns={
        "evidence_status":
            "program_evidence_status"
    }
)

program_lock = pd.read_csv(
    paths["program_lock"]
).rename(
    columns={
        "evidence_status":
            "lock_evidence_status"
    }
)

content = pd.read_csv(
    paths["content"]
).rename(
    columns={
        "content_flag":
            "content_flag_audit",

        "review_priority":
            "review_priority_audit",
    }
)

interpretation = pd.read_csv(
    paths["interpretation"]
).rename(
    columns={
        "content_flag":
            "content_flag_interpretation",

        "review_priority":
            "review_priority_interpretation",
    }
)

overlap = pd.read_csv(
    paths["overlap"]
)

labels = pd.read_csv(
    paths["labels"]
)

claims = pd.read_csv(
    paths["claims"]
)

priority = pd.read_csv(
    paths["priority"]
)

for table, name in [
    (
        effects,
        "B3A effects",
    ),
    (
        program_lock,
        "B3A lock",
    ),
    (
        content,
        "B3B0 content",
    ),
    (
        interpretation,
        "B3B1B interpretation",
    ),
    (
        labels,
        "B3C labels",
    ),
]:
    check_unique(
        table,
        name,
    )

for column in [
    "formal_TLE_program_candidate",
    "directional_support_only",
    "same_direction",
]:
    effects[column] = bool_series(
        effects[column]
    )

for column in [
    "eligible_for_senescence_overlap",
    "eligible_for_competing_analysis",
]:
    content[column] = bool_series(
        content[column]
    )

base = set(
    effects.program_id
)

if (
    len(base) != 25
    or set(program_lock.program_id) != base
    or set(content.program_id) != base
):
    raise RuntimeError(
        "冻结25程序在B3A/B3B0之间不一致。"
    )

if (
    set(interpretation.program_id)
    != set(labels.program_id)
    or len(interpretation) != 11
):
    raise RuntimeError(
        "B3B1B解释表与B3C安全标签表"
        "不一致，或程序数不是11。"
    )

if not set(
    interpretation.program_id
).issubset(base):
    raise RuntimeError(
        "解释程序不是冻结25程序的子集。"
    )

if int(
    effects[
        "formal_TLE_program_candidate"
    ].sum()
) != 0:
    raise RuntimeError(
        "冻结结果不应存在正式TLE相关cNMF程序。"
    )

if (
    "Oligodendrocytes__P5"
    not in base
):
    raise RuntimeError(
        "缺少Oligodendrocytes__P5。"
    )

lineage_incompatible_count = int(
    (
        content[
            "content_flag_audit"
        ]
        == "LIKELY_LINEAGE_INCOMPATIBLE_RNA_PROGRAM"
    ).sum()
)

if lineage_incompatible_count != 8:
    raise RuntimeError(
        "谱系不兼容程序数不是8。"
    )

active_cycle_count = int(
    (
        content[
            "content_flag_audit"
        ]
        == "ACTIVE_CELL_CYCLE_COMPETING_STATE"
    ).sum()
)

if active_cycle_count != 1:
    raise RuntimeError(
        "活动性细胞周期竞争程序数不是1。"
    )


###############################################################################
# 3. Merge one program-level evidence matrix
###############################################################################

lock_keep = program_lock[
    [
        "program_id",
        "k",
        "threshold",
        "cross_sampling_correlation",
        "lock_evidence_status",
        "allowed_for_formal_testing",
    ]
]

content_keep = content[
    [
        "program_id",
        "dominant_general_module",
        "dominant_overlap_top50",
        "expected_lineage_overlap_top50",
        "incompatible_overlap_top50",
        "lnc_fraction_top20",
        "cell_cycle_overlap_top50",
        "stress_overlap_top50",
        "content_flag_audit",
        "eligible_for_senescence_overlap",
        "eligible_for_competing_analysis",
        "review_priority_audit",
    ]
]

interpret_keep = interpretation.drop(
    columns=[
        column
        for column in [
            "lineage",
            "review_priority_interpretation",
        ]
        if column in interpretation.columns
    ]
)

labels_keep = labels.drop(
    columns=[
        column
        for column in [
            "lineage",
            "primary_effect",
            "primary_hedges_g",
            "primary_exact_p",
            "primary_fdr",
            "sensitivity_effect",
            "same_direction",
            "upstream_interpretation",
        ]
        if column in labels.columns
    ]
)

master = effects.merge(
    lock_keep,
    on="program_id",
    validate="one_to_one",
)

master = master.merge(
    content_keep,
    on="program_id",
    validate="one_to_one",
)

master = master.merge(
    interpret_keep,
    on="program_id",
    how="left",
    validate="one_to_one",
)

master = master.merge(
    labels_keep,
    on="program_id",
    how="left",
    validate="one_to_one",
)

if (
    len(master) != 25
    or master.program_id.duplicated().any()
):
    raise RuntimeError(
        "主证据矩阵不是25行唯一程序。"
    )

master[
    "axis_tags_preliminary"
] = master.apply(
    axis_tags,
    axis=1,
)

master[
    "program_use_status"
] = master.apply(
    use_status,
    axis=1,
)

master[
    "has_multiple_axis_tags"
] = master[
    "axis_tags_preliminary"
].str.contains(";")

master[
    "passes_content_qc_for_interpretation"
] = master.program_id.isin(
    set(
        interpretation.program_id
    )
)

master[
    "is_primary_exploratory_program"
] = master.program_id.eq(
    "Oligodendrocytes__P5"
)


###############################################################################
# 4. Pre-register targeted contrastiveVI questions
###############################################################################

targets = pd.DataFrame(
    [
        [
            1,
            "Oligodendrocytes",
            "Oligodendrocytes__P5",
            "PRIMARY_DISEASE_LINKED_EXPLORATORY_TARGET",
            (
                "Test alignment of TLE-salient variation "
                "with frozen P5 usage and its "
                "myelin-lipid-lysosomal/stress spectrum."
            ),
        ],
        [
            2,
            "Astrocytes",
            "Astrocytes__P4",
            "CONTEXTUAL_SENESCENCE_FEATURE_TARGET",
            (
                "Test alignment with the "
                "chaperone-autophagy/proteostasis program "
                "without assuming disease association."
            ),
        ],
        [
            3,
            "Astrocytes",
            "Astrocytes__P6",
            "CONTEXTUAL_SENESCENCE_FEATURE_TARGET",
            (
                "Test alignment with the transport/projection "
                "program carrying senescence-associated features."
            ),
        ],
        [
            4,
            "Microglia_macrophage",
            "Microglia_macrophage__P3",
            "CONTEXTUAL_LIPID_SENESCENCE_FEATURE_TARGET",
            (
                "Test alignment with the myeloid lipid, "
                "oxygen-stress and phagolysosomal program."
            ),
        ],
        [
            5,
            "Astrocytes",
            "Astrocytes__P5",
            "INFLAMMATORY_COMPETING_RESPONSE",
            (
                "Test whether TLE-salient variation is better "
                "explained by inflammatory/immediate-response "
                "signaling."
            ),
        ],
        [
            6,
            "Microglia_macrophage",
            "Microglia_macrophage__P2",
            "PHAGOCYTIC_COMPETING_RESPONSE",
            (
                "Test alignment with motility and "
                "phagocytic signaling."
            ),
        ],
        [
            7,
            "Microglia_macrophage",
            "Microglia_macrophage__P4",
            "INFLAMMATORY_STRESS_COMPETING_RESPONSE",
            (
                "Test alignment with immune-receptor/stress "
                "responses rather than stable cell-cycle exit."
            ),
        ],
    ],
    columns=[
        "priority",
        "lineage",
        "program_id",
        "target_role",
        "contrastivevi_question",
    ],
)

if not set(
    targets.program_id
).issubset(base):
    raise RuntimeError(
        "预注册contrastiveVI目标不存在。"
    )

master[
    "pre_registered_for_contrastivevi"
] = master.program_id.isin(
    set(
        targets.program_id
    )
)

excluded = master[
    master[
        "program_use_status"
    ].str.startswith(
        "EXCLUDE_"
    )
    |
    master[
        "program_use_status"
    ].eq(
        "MANUAL_REVIEW_LOW_ANNOTATABILITY"
    )
].copy()

lineage_summary = (
    master.groupby(
        "lineage",
        dropna=False,
    )
    .agg(
        frozen_programs=(
            "program_id",
            "nunique",
        ),

        content_qc_eligible=(
            "passes_content_qc_for_interpretation",
            "sum",
        ),

        contrastivevi_targets=(
            "pre_registered_for_contrastivevi",
            "sum",
        ),

        directional_support_programs=(
            "directional_support_only",
            "sum",
        ),

        median_cross_sampling_correlation=(
            "cross_sampling_correlation",
            "median",
        ),
    )
    .reset_index()
)

axis_long = pd.DataFrame(
    [
        {
            "lineage":
                row.lineage,

            "program_id":
                row.program_id,

            "axis_tag_preliminary":
                tag,

            "program_use_status":
                row.program_use_status,

            "primary_effect":
                row.primary_effect,

            "primary_exact_p":
                row.primary_exact_p,

            "primary_fdr":
                row.primary_fdr,

            "primary_lodo":
                row.primary_lodo,

            "safe_label":
                getattr(
                    row,
                    "safe_label",
                    "",
                ),
        }
        for row in master.itertuples(
            index=False
        )
        for tag in str(
            row.axis_tags_preliminary
        ).split(";")
    ]
)


###############################################################################
# 5. Freeze outputs
###############################################################################

STATUS = (
    "GO_FOR_34A1_CONTRASTIVEVI_"
    "TARGETED_THREE_LINEAGE_ANALYSIS"
)

OUT.mkdir(
    parents=True,
    exist_ok=False,
)

master.sort_values(
    [
        "lineage",
        "program_id",
    ]
).to_csv(
    OUT
    / "34A0_MASTER_PROGRAM_EVIDENCE_MATRIX.csv",
    index=False,
)

axis_long.sort_values(
    [
        "lineage",
        "program_id",
        "axis_tag_preliminary",
    ]
).to_csv(
    OUT
    / "34A0_PROGRAM_AXIS_TAGS_LONG.csv",
    index=False,
)

targets.sort_values(
    "priority"
).to_csv(
    OUT
    / "34A0_CONTRASTIVEVI_TARGET_REGISTRY.csv",
    index=False,
)

excluded.sort_values(
    [
        "lineage",
        "program_id",
    ]
).to_csv(
    OUT
    / "34A0_EXCLUDED_OR_MANUAL_REVIEW_PROGRAMS.csv",
    index=False,
)

lineage_summary.to_csv(
    OUT
    / "34A0_LINEAGE_SUMMARY.csv",
    index=False,
)

claims.to_csv(
    OUT
    / "34A0_UPSTREAM_CLAIM_MATRIX.csv",
    index=False,
)

priority.to_csv(
    OUT
    / "34A0_UPSTREAM_VALIDATION_PRIORITY.csv",
    index=False,
)

overlap.to_csv(
    OUT
    / "34A0_FROZEN_26E_OVERLAP_WITH_DISEASE_DIRECTION.csv",
    index=False,
)

provenance = {
    "stage":
        "34A0",

    "created_at":
        datetime.now().isoformat(),

    "B3A_release":
        str(B3A),

    "B3B0_release":
        str(B3B0),

    "B3B1B_release":
        str(B3B1B),

    "B3C_release":
        str(B3C),

    "input_files": {
        name: {
            "path":
                str(path),

            "sha256":
                sha256(path),
        }
        for name, path in paths.items()
    },

    "new_biological_test_performed":
        False,

    "program_refitted":
        False,

    "program_gene_membership_changed":
        False,

    "contrastivevi_run":
        False,
}

(
    OUT
    / "34A0_INPUT_PROVENANCE.json"
).write_text(
    json.dumps(
        provenance,
        ensure_ascii=False,
        indent=2,
    ),
    encoding="utf-8",
)

lock_json = {
    "stage":
        "34A0",

    "created_at":
        datetime.now().isoformat(),

    "status":
        STATUS,

    "frozen_programs":
        int(
            len(master)
        ),

    "content_qc_eligible_programs":
        int(
            master[
                "passes_content_qc_for_interpretation"
            ].sum()
        ),

    "formal_TLE_programs":
        int(
            master[
                "formal_TLE_program_candidate"
            ].sum()
        ),

    "primary_exploratory_program":
        "Oligodendrocytes__P5",

    "contrastivevi_lineages": [
        "Astrocytes",
        "Microglia_macrophage",
        "Oligodendrocytes",
    ],

    "pre_registered_contrastivevi_programs":
        targets.program_id.tolist(),

    "axis_tags_are_preliminary_not_final_classification":
        True,

    "no_program_reselection_after_contrastivevi":
        True,

    "donor_is_inferential_unit":
        True,

    "cell_level_pseudoreplication_authorized":
        False,

    "senescent_cell_label_authorized":
        False,

    "causal_claim_authorized":
        False,

    "next_stage":
        "34A1_CONTRASTIVEVI_TARGETED_THREE_LINEAGE_ANALYSIS",
}

(
    OUT
    / "34A0_STATE_CONSENSUS_LOCK.json"
).write_text(
    json.dumps(
        lock_json,
        ensure_ascii=False,
        indent=2,
    ),
    encoding="utf-8",
)

summary = f"""# 34A-0 Glial-state Consensus Input Audit

## Status: {STATUS}

## Completed

- Integrated frozen cNMF disease effects, content QC, 26E interpretation and manuscript-safe labels.
- Created one 25-program evidence matrix.
- Added preliminary parallel tags for partial senescence, inflammation/stress and phagocytic-lipid-metabolic remodeling.
- Pre-registered seven targeted contrastiveVI questions before fitting the model.

## Not performed

- No new disease test.
- No program refitting or gene reselection.
- No final biological-state classification.
- No contrastiveVI analysis.

## Frozen counts

- Frozen cNMF programs: {len(master)}
- Programs passing content QC for 26E interpretation: {int(master.passes_content_qc_for_interpretation.sum())}
- Formal TLE-associated cNMF programs: {int(master.formal_TLE_program_candidate.sum())}
- Pre-registered contrastiveVI targets: {len(targets)}

## Primary exploratory program

- Oligodendrocytes__P5 remains exploratory and mixed; it is not a senescent-cell label.

## Next stage

{lock_json["next_stage"]}
"""

(
    OUT
    / "34A0_EXECUTIVE_SUMMARY.md"
).write_text(
    summary,
    encoding="utf-8",
)

(
    OUT
    / "34A0_python_session.txt"
).write_text(
    (
        f"python={sys.version}\n"
        f"pandas={pd.__version__}\n"
    ),
    encoding="utf-8",
)

output_files = sorted(
    path
    for path in OUT.rglob("*")
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
                sha256(path),
        }
        for path in output_files
    ]
).to_csv(
    OUT
    / "34A0_FILE_MANIFEST.csv",
    index=False,
)

zip_path = Path(
    str(OUT)
    + ".zip"
)

with zipfile.ZipFile(
    zip_path,
    "w",
    zipfile.ZIP_DEFLATED,
    compresslevel=9,
) as archive:

    for path in sorted(
        OUT.rglob("*")
    ):
        if path.is_file():
            archive.write(
                path,
                path.relative_to(
                    OUT
                ).as_posix(),
            )

zip_hash = sha256(
    zip_path
)

Path(
    str(zip_path)
    + ".sha256.txt"
).write_text(
    (
        f"{zip_hash}  "
        f"{zip_path.name}\n"
    ),
    encoding="ascii",
)

print(
    "=" * 78
)

print(
    "34A-0 completed"
)

print(
    "=" * 78
)

print(
    "状态：",
    STATUS,
)

print(
    "冻结程序数：",
    len(master),
)

print(
    "内容QC后可解释程序数：",
    int(
        master[
            "passes_content_qc_for_interpretation"
        ].sum()
    ),
)

print(
    "正式TLE程序数：",
    int(
        master[
            "formal_TLE_program_candidate"
        ].sum()
    ),
)

print(
    "contrastiveVI预注册目标数：",
    len(targets),
)

print(
    "结果目录：",
    OUT,
)

print(
    "ZIP：",
    zip_path,
)

print(
    "SHA256：",
    zip_hash,
)

print(
    "=" * 78
)