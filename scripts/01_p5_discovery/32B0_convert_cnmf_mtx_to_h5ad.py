from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


PROJECT = Path(
    __import__("os").environ["P5_WORK_ROOT"]
)

INPUT_ROOT = (
    PROJECT
    / "99_frozen_releases"
    / "32B0_BALANCED_CNMf_INPUTS_20260728_215801"
    / "cnmf_input"
)

OUTPUT_ROOT = (
    PROJECT
    / "04_tables"
    / "32B_CNMf_WORK_20260728_215801"
    / "_h5ad_inputs"
)

LINEAGES = {
    "Astrocytes": "Astrocytes",
    "Microglia_macrophage": "Microglia_macrophage",
    "OPC": "OPC",
    "Oligodendrocytes": "Oligodendrocytes",
}


def read_one_column(path: Path) -> list[str]:
    table = pd.read_csv(
        path,
        sep="\t",
        header=None,
        dtype=str,
    )
    return table.iloc[:, 0].tolist()


def convert_lineage(folder_name: str, output_name: str) -> None:
    input_dir = INPUT_ROOT / folder_name / "balanced_primary"

    matrix_path = input_dir / "counts.mtx"
    barcodes_path = input_dir / "barcodes.tsv"
    genes_path = input_dir / "genes.tsv"
    metadata_path = input_dir / "cell_metadata.csv"

    required = [
        matrix_path,
        barcodes_path,
        genes_path,
        metadata_path,
    ]

    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing input files:\n" + "\n".join(missing)
        )

    print(f"\nReading {folder_name}...")

    # R写出的Matrix Market矩阵为 genes × cells；
    # AnnData要求 cells × genes，所以需要转置。
    matrix = mmread(matrix_path)

    if not sparse.issparse(matrix):
        matrix = sparse.coo_matrix(matrix)

    matrix = matrix.tocsr().transpose().tocsr()

    barcodes = read_one_column(barcodes_path)

    genes_table = pd.read_csv(
        genes_path,
        sep="\t",
        header=None,
        dtype=str,
    )

    gene_ids = genes_table.iloc[:, 0].tolist()

    if genes_table.shape[1] >= 2:
        gene_names = genes_table.iloc[:, 1].tolist()
    else:
        gene_names = gene_ids

    if matrix.shape != (len(barcodes), len(gene_ids)):
        raise ValueError(
            f"{folder_name}: matrix shape {matrix.shape} "
            f"does not match cells={len(barcodes)}, "
            f"genes={len(gene_ids)}"
        )

    metadata = pd.read_csv(
        metadata_path,
        dtype=str,
        index_col="cell_id",
    )

    missing_metadata = set(barcodes) - set(metadata.index)
    if missing_metadata:
        raise ValueError(
            f"{folder_name}: "
            f"{len(missing_metadata)} barcodes missing from metadata"
        )

    metadata = metadata.loc[barcodes].copy()

    var = pd.DataFrame(
        {
            "gene_id": gene_ids,
            "gene_name": gene_names,
        },
        index=pd.Index(gene_names, name="gene"),
    )

    adata = ad.AnnData(
        X=matrix,
        obs=metadata,
        var=var,
    )

    adata.obs_names = pd.Index(barcodes, name="cell_id")
    adata.var_names_make_unique()

    if np.any(adata.X.data < 0):
        raise ValueError(f"{folder_name}: negative counts detected")

    if np.any(np.abs(adata.X.data - np.rint(adata.X.data)) > 1e-8):
        raise ValueError(f"{folder_name}: non-integer counts detected")

    output_path = OUTPUT_ROOT / f"{output_name}_primary.h5ad"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    adata.write_h5ad(
        output_path,
        compression="gzip",
    )

    print(
        f"Saved: {output_path}\n"
        f"Shape: {adata.n_obs} cells × {adata.n_vars} genes"
    )


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    for folder_name, output_name in LINEAGES.items():
        convert_lineage(folder_name, output_name)

    print("\nAll four H5AD inputs were created successfully.")


if __name__ == "__main__":
    main()