from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


TARGET_DONORS = [
    "H18.30.001",
    "H18.30.002",
    "H19.30.001",
    "H19.30.002",
    "H200.1023",
]

MIN_FORMAL_CELLS = 30


def clean_series(series: pd.Series) -> pd.Series:
    return (
        series.astype("string")
        .fillna("")
        .str.strip()
    )


def write_csv(
    path: Path,
    rows: list[dict[str, Any]],
    fieldnames: list[str],
) -> None:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open(
        "w",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def detect_donor_column(
    obs: pd.DataFrame,
) -> tuple[str, list[dict[str, Any]]]:
    audit_rows: list[dict[str, Any]] = []

    for column in obs.columns:
        values = clean_series(obs[column])

        target_mask = values.isin(TARGET_DONORS)

        matched_cells = int(
            target_mask.sum()
        )

        matched_donors = sorted(
            set(values[target_mask].tolist())
        )

        name_lower = str(column).lower()

        name_score = int(
            any(
                token in name_lower
                for token in (
                    "donor",
                    "subject",
                    "individual",
                    "participant",
                    "patient",
                )
            )
        )

        audit_rows.append(
            {
                "column": str(column),
                "matched_cells": matched_cells,
                "matched_target_donors":
                    len(matched_donors),
                "target_donor_values":
                    "|".join(matched_donors),
                "donor_like_column_name":
                    bool(name_score),
            }
        )

    ranked = sorted(
        audit_rows,
        key=lambda row: (
            -int(row["matched_target_donors"]),
            -int(row["matched_cells"]),
            -int(row["donor_like_column_name"]),
            str(row["column"]),
        ),
    )

    if not ranked:
        raise RuntimeError(
            "obs中没有任何元数据列。"
        )

    best = ranked[0]

    if int(best["matched_target_donors"]) == 0:
        raise RuntimeError(
            "没有任何obs字段包含目标SEA-AD供者ID。"
        )

    return str(best["column"]), ranked


def get_sample_values(
    matrix: Any,
    n_rows: int,
    n_cols: int,
) -> np.ndarray:
    row_count = min(200, n_rows)
    col_count = min(1000, n_cols)

    if row_count == 0 or col_count == 0:
        return np.asarray([], dtype=float)

    block = matrix[
        :row_count,
        :col_count,
    ]

    if hasattr(block, "to_memory"):
        block = block.to_memory()

    if sp.issparse(block):
        values = np.asarray(
            block.data,
            dtype=float,
        )
    else:
        values = np.asarray(
            block,
            dtype=float,
        ).ravel()

        values = values[
            values != 0
        ]

    values = values[
        np.isfinite(values)
    ]

    return values


def inspect_count_matrix(
    matrix: Any,
    n_rows: int,
    n_cols: int,
) -> dict[str, Any]:
    values = get_sample_values(
        matrix=matrix,
        n_rows=n_rows,
        n_cols=n_cols,
    )

    if values.size == 0:
        return {
            "nonzero_values_sampled": 0,
            "minimum": "",
            "maximum": "",
            "negative_values": 0,
            "integer_like_fraction": "",
            "raw_count_eligible": False,
        }

    integer_like_fraction = float(
        np.mean(
            np.isclose(
                values,
                np.round(values),
                atol=1e-8,
            )
        )
    )

    negative_values = int(
        np.sum(values < 0)
    )

    eligible = bool(
        negative_values == 0
        and integer_like_fraction >= 0.99
    )

    return {
        "nonzero_values_sampled":
            int(values.size),
        "minimum":
            float(np.min(values)),
        "maximum":
            float(np.max(values)),
        "negative_values":
            negative_values,
        "integer_like_fraction":
            integer_like_fraction,
        "raw_count_eligible":
            eligible,
    }


def get_lineage_label(
    file_path: Path,
) -> str:
    lower = file_path.name.lower()

    if "opc" in lower:
        return "OPC"

    if "oligodendrocyte" in lower:
        return "Oligodendrocyte"

    raise RuntimeError(
        f"无法从文件名判断细胞谱系：{file_path.name}"
    )


def materialize_matrix(
    matrix: Any,
) -> Any:
    if hasattr(matrix, "to_memory"):
        matrix = matrix.to_memory()

    if sp.issparse(matrix):
        return matrix.tocsr().copy()

    return np.asarray(matrix).copy()


def create_counts_adata(
    subset: ad.AnnData,
    matrix_source: str,
    donor_column: str,
    lineage: str,
) -> ad.AnnData:
    if matrix_source == "raw.X":
        if subset.raw is None:
            raise RuntimeError(
                "预期raw存在，但子集raw为空。"
            )

        counts = materialize_matrix(
            subset.raw.X
        )

        var = subset.raw.var.copy()

    elif matrix_source == "X":
        counts = materialize_matrix(
            subset.X
        )

        var = subset.var.copy()

    else:
        raise RuntimeError(
            f"不支持的矩阵来源：{matrix_source}"
        )

    obs = subset.obs.copy()

    donor_values = clean_series(
        obs[donor_column]
    )

    obs["donor_id_35c"] = donor_values
    obs["lineage_35c"] = lineage

    output = ad.AnnData(
        X=counts,
        obs=obs,
        var=var,
    )

    output.uns["35c_matrix_source"] = (
        matrix_source
    )

    output.uns["35c_donor_column"] = (
        donor_column
    )

    output.uns["35c_lineage"] = lineage
    output.uns["35c_target_donors"] = (
        TARGET_DONORS
    )

    output.uns["35c_min_formal_cells"] = (
        MIN_FORMAL_CELLS
    )

    return output


def process_file(
    file_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    lineage = get_lineage_label(
        file_path
    )

    print("")
    print(
        f"处理{lineage}：{file_path.name}",
        flush=True,
    )

    adata = ad.read_h5ad(
        file_path,
        backed="r",
    )

    try:
        obs = adata.obs.copy()

        donor_column, donor_audit = (
            detect_donor_column(obs)
        )

        print(
            f"  donor字段：{donor_column}",
            flush=True,
        )

        donor_values = clean_series(
            obs[donor_column]
        )

        counts_by_donor = (
            donor_values[
                donor_values.isin(
                    TARGET_DONORS
                )
            ]
            .value_counts()
            .reindex(
                TARGET_DONORS,
                fill_value=0,
            )
        )

        count_rows = []

        for donor in TARGET_DONORS:
            cell_count = int(
                counts_by_donor.loc[donor]
            )

            count_rows.append(
                {
                    "lineage": lineage,
                    "donor_id": donor,
                    "cell_count": cell_count,
                    "formal_eligible":
                        cell_count
                        >= MIN_FORMAL_CELLS,
                }
            )

        raw_audit = None

        if adata.raw is not None:
            raw_audit = inspect_count_matrix(
                matrix=adata.raw.X,
                n_rows=adata.raw.n_obs,
                n_cols=adata.raw.n_vars,
            )

        x_audit = inspect_count_matrix(
            matrix=adata.X,
            n_rows=adata.n_obs,
            n_cols=adata.n_vars,
        )

        if (
            raw_audit is not None
            and raw_audit[
                "raw_count_eligible"
            ]
        ):
            matrix_source = "raw.X"
            selected_matrix_audit = raw_audit

        elif x_audit["raw_count_eligible"]:
            matrix_source = "X"
            selected_matrix_audit = x_audit

        else:
            raise RuntimeError(
                f"{file_path.name}中没有通过整数计数检查的矩阵。"
            )

        print(
            f"  原始计数来源：{matrix_source}",
            flush=True,
        )

        target_mask = (
            donor_values.isin(
                TARGET_DONORS
            )
            .to_numpy()
        )

        if int(target_mask.sum()) == 0:
            raise RuntimeError(
                f"{file_path.name}没有目标供者细胞。"
            )

        target_subset = (
            adata[target_mask]
            .to_memory()
        )

        target_counts = create_counts_adata(
            subset=target_subset,
            matrix_source=matrix_source,
            donor_column=donor_column,
            lineage=lineage,
        )

        formal_donors = [
            donor
            for donor in TARGET_DONORS
            if int(counts_by_donor.loc[donor])
            >= MIN_FORMAL_CELLS
        ]

        formal_mask = (
            target_counts.obs[
                "donor_id_35c"
            ]
            .astype(str)
            .isin(formal_donors)
            .to_numpy()
        )

        formal_counts = (
            target_counts[formal_mask]
            .copy()
        )

        target_output = (
            output_dir
            / f"SEAAD_{lineage.upper()}_TARGET5_RAWCOUNTS.h5ad"
        )

        formal_output = (
            output_dir
            / f"SEAAD_{lineage.upper()}_FORMAL_RAWCOUNTS.h5ad"
        )

        target_counts.write_h5ad(
            target_output,
            compression="gzip",
        )

        formal_counts.write_h5ad(
            formal_output,
            compression="gzip",
        )

        donor_audit_rows = [
            {
                "lineage": lineage,
                **row,
            }
            for row in donor_audit
        ]

        matrix_rows = [
            {
                "lineage": lineage,
                "file_name": file_path.name,
                "matrix": "X",
                **x_audit,
                "selected_for_counts":
                    matrix_source == "X",
            }
        ]

        if raw_audit is not None:
            matrix_rows.append(
                {
                    "lineage": lineage,
                    "file_name":
                        file_path.name,
                    "matrix": "raw.X",
                    **raw_audit,
                    "selected_for_counts":
                        matrix_source
                        == "raw.X",
                }
            )

        return {
            "lineage": lineage,
            "donor_column": donor_column,
            "donor_counts": count_rows,
            "donor_audit":
                donor_audit_rows,
            "matrix_audit": matrix_rows,
            "formal_donors":
                formal_donors,
            "target_output":
                str(target_output),
            "formal_output":
                str(formal_output),
            "target_cells":
                int(target_counts.n_obs),
            "formal_cells":
                int(formal_counts.n_obs),
            "genes":
                int(target_counts.n_vars),
            "matrix_source":
                matrix_source,
            "selected_matrix_audit":
                selected_matrix_audit,
        }

    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input-dir",
        required=True,
    )

    parser.add_argument(
        "--output-dir",
        required=True,
    )

    args = parser.parse_args()

    input_dir = Path(
        args.input_dir
    ).resolve()

    output_dir = Path(
        args.output_dir
    ).resolve()

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    files = sorted(
        input_dir.glob("*.h5ad")
    )

    if len(files) != 2:
        raise RuntimeError(
            f"预期2个H5AD，实际发现{len(files)}个。"
        )

    results = [
        process_file(
            file_path=file_path,
            output_dir=output_dir,
        )
        for file_path in files
    ]

    donor_count_rows = [
        row
        for result in results
        for row in result["donor_counts"]
    ]

    donor_audit_rows = [
        row
        for result in results
        for row in result["donor_audit"]
    ]

    matrix_audit_rows = [
        row
        for result in results
        for row in result["matrix_audit"]
    ]

    write_csv(
        output_dir /
        "SEAAD_TARGET_DONOR_CELL_COUNTS_CONFIRMED.csv",
        donor_count_rows,
        [
            "lineage",
            "donor_id",
            "cell_count",
            "formal_eligible",
        ],
    )

    write_csv(
        output_dir /
        "SEAAD_DONOR_COLUMN_AUDIT.csv",
        donor_audit_rows,
        [
            "lineage",
            "column",
            "matched_cells",
            "matched_target_donors",
            "target_donor_values",
            "donor_like_column_name",
        ],
    )

    write_csv(
        output_dir /
        "SEAAD_MATRIX_SOURCE_AUDIT.csv",
        matrix_audit_rows,
        [
            "lineage",
            "file_name",
            "matrix",
            "nonzero_values_sampled",
            "minimum",
            "maximum",
            "negative_values",
            "integer_like_fraction",
            "raw_count_eligible",
            "selected_for_counts",
        ],
    )

    formal_sets = {
        result["lineage"]:
            set(result["formal_donors"])
        for result in results
    }

    common_formal_donors = sorted(
        set.intersection(
            *formal_sets.values()
        )
    )

    formal_rows = []

    for donor in TARGET_DONORS:
        row = {
            "donor_id": donor,
            "formal_in_OPC":
                donor
                in formal_sets.get(
                    "OPC",
                    set(),
                ),
            "formal_in_Oligodendrocyte":
                donor
                in formal_sets.get(
                    "Oligodendrocyte",
                    set(),
                ),
            "formal_in_both":
                donor
                in common_formal_donors,
        }

        formal_rows.append(row)

    write_csv(
        output_dir /
        "SEAAD_FORMAL_REFERENCE_DONORS.csv",
        formal_rows,
        [
            "donor_id",
            "formal_in_OPC",
            "formal_in_Oligodendrocyte",
            "formal_in_both",
        ],
    )

    summary = {
        "target_donors":
            TARGET_DONORS,
        "minimum_formal_cells":
            MIN_FORMAL_CELLS,
        "common_formal_donors":
            common_formal_donors,
        "results": [
            {
                "lineage":
                    result["lineage"],
                "donor_column":
                    result["donor_column"],
                "matrix_source":
                    result["matrix_source"],
                "target_cells":
                    result["target_cells"],
                "formal_cells":
                    result["formal_cells"],
                "genes":
                    result["genes"],
                "formal_donors":
                    result["formal_donors"],
                "target_output":
                    result["target_output"],
                "formal_output":
                    result["formal_output"],
            }
            for result in results
        ],
    }

    (
        output_dir /
        "SEAAD_SUBSET_EXTRACTION_SUMMARY.json"
    ).write_text(
        json.dumps(
            summary,
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print("")
    print("SEA-AD目标供者提取完成")
    print("")

    for result in results:
        print(
            f"{result['lineage']}:"
        )
        print(
            f"  counts来源="
            f"{result['matrix_source']}"
        )
        print(
            f"  5名目标供者细胞="
            f"{result['target_cells']}"
        )
        print(
            f"  正式合格供者细胞="
            f"{result['formal_cells']}"
        )
        print(
            f"  基因数="
            f"{result['genes']}"
        )
        print(
            "  正式供者="
            + ", ".join(
                result["formal_donors"]
            )
        )
        print("")

    print(
        "OPC与成熟少突共同合格供者："
        + ", ".join(
            common_formal_donors
        )
    )

    print("")
    print(f"输出目录：{output_dir}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(
            f"FATAL: {type(error).__name__}: {error}",
            file=sys.stderr,
        )
        raise