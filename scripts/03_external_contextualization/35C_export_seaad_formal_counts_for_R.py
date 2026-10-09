from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import shutil
import tempfile
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


EXPECTED_FORMAL_DONORS = {
    "H18.30.002",
    "H19.30.001",
    "H19.30.002",
}

GENE_SYMBOL_COLUMNS = (
    "feature_name",
    "gene_symbol",
    "gene_symbols",
    "gene_name",
    "symbol",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def clean_text(value: Any) -> str:
    if value is None:
        return ""

    return str(value).strip()


def choose_gene_symbols(
    var: pd.DataFrame,
) -> tuple[pd.Index, str]:
    for column in GENE_SYMBOL_COLUMNS:
        if column not in var.columns:
            continue

        values = (
            var[column]
            .astype("string")
            .fillna("")
            .str.strip()
        )

        nonempty_fraction = float(
            np.mean(values != "")
        )

        if nonempty_fraction >= 0.90:
            return pd.Index(values), column

    return pd.Index(
        var.index.astype(str)
    ), "var_names"


def materialize_sparse_matrix(
    matrix: Any,
) -> sp.csr_matrix:
    if hasattr(matrix, "to_memory"):
        matrix = matrix.to_memory()

    if sp.issparse(matrix):
        output = matrix.tocsr()
    else:
        output = sp.csr_matrix(
            np.asarray(matrix)
        )

    output.eliminate_zeros()

    return output


def validate_counts(
    matrix: sp.csr_matrix,
) -> dict[str, Any]:
    values = matrix.data

    if values.size == 0:
        raise RuntimeError(
            "计数矩阵没有任何非零值。"
        )

    if not np.all(np.isfinite(values)):
        raise RuntimeError(
            "计数矩阵包含NaN或无限值。"
        )

    negative_count = int(
        np.sum(values < 0)
    )

    integer_like_fraction = float(
        np.mean(
            np.isclose(
                values,
                np.round(values),
                atol=1e-8,
            )
        )
    )

    if negative_count > 0:
        raise RuntimeError(
            "计数矩阵包含负值。"
        )

    if integer_like_fraction < 0.99:
        raise RuntimeError(
            "矩阵不像原始整数计数；"
            f"integer-like fraction={integer_like_fraction:.6f}"
        )

    return {
        "nonzero_count": int(values.size),
        "minimum_nonzero": float(values.min()),
        "maximum": float(values.max()),
        "integer_like_fraction": integer_like_fraction,
        "negative_values": negative_count,
    }


def write_gzip_lines(
    path: Path,
    rows: list[str],
) -> None:
    with gzip.open(
        path,
        "wt",
        encoding="utf-8",
        newline="\n",
    ) as handle:
        for row in rows:
            handle.write(row)
            handle.write("\n")


def write_matrix_market_gzip(
    matrix_gene_by_cell: sp.spmatrix,
    output_path: Path,
) -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        temporary_mtx = (
            Path(temp_dir) / "matrix.mtx"
        )

        scipy.io.mmwrite(
            temporary_mtx,
            matrix_gene_by_cell,
            field="integer",
            symmetry="general",
        )

        with temporary_mtx.open("rb") as source:
            with gzip.open(
                output_path,
                "wb",
                compresslevel=6,
            ) as destination:
                shutil.copyfileobj(
                    source,
                    destination,
                    length=1024 * 1024,
                )


def detect_lineage(file_path: Path) -> str:
    name = file_path.name.lower()

    if "oligodendrocyte" in name:
        return "Mature_oligodendrocyte"

    if "opc" in name:
        return "OPC"

    raise RuntimeError(
        f"无法从文件名判断谱系：{file_path.name}"
    )


def export_one_h5ad(
    file_path: Path,
    output_root: Path,
) -> dict[str, Any]:
    lineage = detect_lineage(file_path)

    output_dir = output_root / lineage
    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("")
    print(f"读取：{file_path.name}")
    print(f"谱系：{lineage}")

    adata = ad.read_h5ad(file_path)

    if "donor_id_35c" not in adata.obs.columns:
        raise RuntimeError(
            f"{file_path.name}缺少donor_id_35c。"
        )

    donor_values = (
        adata.obs["donor_id_35c"]
        .astype("string")
        .fillna("")
        .str.strip()
    )

    observed_donors = set(
        donor_values.unique()
    )

    if observed_donors != EXPECTED_FORMAL_DONORS:
        raise RuntimeError(
            f"{lineage}正式供者不符合预期。"
            f"观察到：{sorted(observed_donors)}"
        )

    counts_cell_by_gene = (
        materialize_sparse_matrix(
            adata.X
        )
    )

    if counts_cell_by_gene.shape != (
        adata.n_obs,
        adata.n_vars,
    ):
        raise RuntimeError(
            "表达矩阵维度与AnnData不一致。"
        )

    count_audit = validate_counts(
        counts_cell_by_gene
    )

    gene_symbols, symbol_source = (
        choose_gene_symbols(
            adata.var
        )
    )

    if len(gene_symbols) != adata.n_vars:
        raise RuntimeError(
            "基因符号数量与基因数不一致。"
        )

    cell_barcodes = pd.Index(
        adata.obs_names.astype(str)
    )

    if cell_barcodes.has_duplicates:
        raise RuntimeError(
            "细胞barcode存在重复。"
        )

    matrix_path = (
        output_dir / "matrix.mtx.gz"
    )

    feature_path = (
        output_dir / "features.tsv.gz"
    )

    barcode_path = (
        output_dir / "barcodes.tsv.gz"
    )

    metadata_path = (
        output_dir / "cell_metadata.csv.gz"
    )

    var_metadata_path = (
        output_dir / "gene_metadata.csv.gz"
    )

    print(
        "写出稀疏矩阵："
        f"{adata.n_vars} genes × "
        f"{adata.n_obs} cells"
    )

    counts_gene_by_cell = (
        counts_cell_by_gene.transpose()
        .tocsc()
    )

    write_matrix_market_gzip(
        counts_gene_by_cell,
        matrix_path,
    )

    feature_rows = [
        f"{gene_id}\t{gene_symbol}"
        for gene_id, gene_symbol in zip(
            adata.var_names.astype(str),
            gene_symbols.astype(str),
        )
    ]

    write_gzip_lines(
        feature_path,
        feature_rows,
    )

    write_gzip_lines(
        barcode_path,
        cell_barcodes.astype(str).tolist(),
    )

    cell_metadata = (
        adata.obs.copy()
        .reset_index()
        .rename(
            columns={
                adata.obs.index.name or "index":
                    "cell_barcode"
            }
        )
    )

    if "cell_barcode" not in cell_metadata.columns:
        cell_metadata.insert(
            0,
            "cell_barcode",
            cell_barcodes.astype(str),
        )

    cell_metadata["lineage_35c"] = lineage

    cell_metadata.to_csv(
        metadata_path,
        index=False,
        compression="gzip",
    )

    gene_metadata = (
        adata.var.copy()
        .reset_index()
    )

    first_column = gene_metadata.columns[0]

    gene_metadata = gene_metadata.rename(
        columns={
            first_column: "gene_id"
        }
    )

    gene_metadata["gene_symbol_35c"] = (
        gene_symbols.astype(str).to_numpy()
    )

    gene_metadata[
        "gene_symbol_source_35c"
    ] = symbol_source

    gene_metadata.to_csv(
        var_metadata_path,
        index=False,
        compression="gzip",
    )

    donor_counts = (
        donor_values
        .value_counts()
        .sort_index()
    )

    donor_count_path = (
        output_dir / "donor_cell_counts.csv"
    )

    donor_count_table = pd.DataFrame(
        {
            "donor_id": donor_counts.index,
            "cell_count": donor_counts.values,
            "lineage": lineage,
        }
    )

    donor_count_table.to_csv(
        donor_count_path,
        index=False,
    )

    files = [
        matrix_path,
        feature_path,
        barcode_path,
        metadata_path,
        var_metadata_path,
        donor_count_path,
    ]

    manifest_rows = []

    for output_file in files:
        manifest_rows.append(
            {
                "lineage": lineage,
                "file_name": output_file.name,
                "full_path": str(output_file),
                "size_bytes":
                    output_file.stat().st_size,
                "sha256":
                    sha256_file(output_file),
            }
        )

    manifest_path = (
        output_dir / "export_manifest.csv"
    )

    pd.DataFrame(
        manifest_rows
    ).to_csv(
        manifest_path,
        index=False,
    )

    print(
        "供者细胞数："
        + ", ".join(
            f"{donor}={count}"
            for donor, count in donor_counts.items()
        )
    )

    return {
        "lineage": lineage,
        "input_file": str(file_path),
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "gene_symbol_source": symbol_source,
        **count_audit,
        "matrix_path": str(matrix_path),
        "metadata_path": str(metadata_path),
    }


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

    formal_files = sorted(
        input_dir.glob(
            "SEAAD_*_FORMAL_RAWCOUNTS.h5ad"
        )
    )

    if len(formal_files) != 2:
        raise RuntimeError(
            "预期2个正式H5AD文件，"
            f"实际找到{len(formal_files)}个。"
        )

    results = [
        export_one_h5ad(
            file_path=file_path,
            output_root=output_dir,
        )
        for file_path in formal_files
    ]

    summary_path = (
        output_dir /
        "SEAAD_R_EXPORT_SUMMARY.csv"
    )

    with summary_path.open(
        "w",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(results[0].keys()),
        )
        writer.writeheader()
        writer.writerows(results)

    print("")
    print("========================================")
    print("SEA-AD正式供者R矩阵导出完成")
    print("========================================")

    for result in results:
        print("")
        print(result["lineage"])
        print(f"  cells={result['cells']}")
        print(f"  genes={result['genes']}")
        print(
            "  gene_symbol_source="
            f"{result['gene_symbol_source']}"
        )
        print(
            "  integer_like_fraction="
            f"{result['integer_like_fraction']}"
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
            flush=True,
        )
        raise