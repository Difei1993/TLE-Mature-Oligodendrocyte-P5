###############################################################################
# 32B-2C
# Filtered cNMF consensus and threshold-sensitivity audit
#
# Reuses completed 100-replicate NMF results from 32B-2B.
# Does NOT rerun factorization.
# Does NOT use disease labels or senescence signatures.
###############################################################################

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import subprocess
import time
import zipfile
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


###############################################################################
# 0. Paths
###############################################################################

PROJECT = Path(
    __import__("os").environ["P5_WORK_ROOT"]
)

ENV_ROOT = (
    PROJECT
    / "04_tables"
    / "python_envs"
    / "cnmf310"
)

CNMF_EXE = ENV_ROOT / "Scripts" / "cnmf.exe"

SOURCE_ROOT = (
    PROJECT
    / "04_tables"
    / "32B2B_FINAL_CNMf_20260728_215801"
)

OUTPUT_ROOT = (
    PROJECT
    / "04_tables"
    / "32B2C_FILTERED_CONSENSUS_AUDIT_20260728_215801"
)

LOG_ROOT = OUTPUT_ROOT / "_logs"
REVIEW_ROOT = OUTPUT_ROOT / "_FILTERED_CONSENSUS_REVIEW"


###############################################################################
# 1. Frozen plan
###############################################################################

PLAN = {
    "Astrocytes": {
        "k": 7,
        "candidate_threshold": 0.02,
        "alternative_thresholds": [0.05],
    },
    "Microglia_macrophage": {
        "k": 6,
        "candidate_threshold": 0.01,
        "alternative_thresholds": [],
    },
    "OPC": {
        "k": 6,
        "candidate_threshold": 0.01,
        "alternative_thresholds": [0.05],
    },
    "Oligodendrocytes": {
        "k": 6,
        "candidate_threshold": 0.01,
        "alternative_thresholds": [],
    },
}

SAMPLINGS = ("primary", "sensitivity")


###############################################################################
# 2. Utilities
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


def threshold_tag(threshold: float) -> str:
    return str(threshold).replace(".", "_")


def run_command(
    label: str,
    arguments: list[str],
    log_path: Path,
) -> None:

    log_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    command = [
        str(CNMF_EXE),
        *arguments,
    ]

    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "1"
    environment["MKL_NUM_THREADS"] = "1"
    environment["OPENBLAS_NUM_THREADS"] = "1"
    environment["NUMEXPR_NUM_THREADS"] = "1"

    started = time.time()

    print()
    print("=" * 72)
    print(label)
    print(f"Started: {now_string()}")
    print("=" * 72)

    with log_path.open(
        "w",
        encoding="utf-8",
        errors="replace",
    ) as log_handle:

        log_handle.write(
            "Command:\n"
            + subprocess.list2cmdline(command)
            + "\n\n"
        )

        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            env=environment,
        )

        assert process.stdout is not None

        for line in process.stdout:
            print(line, end="")
            log_handle.write(line)
            log_handle.flush()

        return_code = process.wait()

    elapsed = time.time() - started

    if return_code != 0:
        raise RuntimeError(
            f"{label} failed with exit code {return_code}. "
            f"See: {log_path}"
        )

    print(
        f"Completed: {label}\n"
        f"Elapsed: {elapsed:.2f} seconds"
    )


def locate_one(
    run_directory: Path,
    pattern: str,
) -> Path:

    matches = sorted(
        path
        for path in run_directory.rglob(pattern)
        if path.is_file()
    )

    if not matches:
        raise FileNotFoundError(
            f"No file matching {pattern} in {run_directory}"
        )

    return matches[0]


def load_local_density(
    run_directory: Path,
    k_value: int,
) -> np.ndarray:

    cache_path = locate_one(
        run_directory,
        f"*.local_density_cache.k_{k_value}.merged.df.npz",
    )

    cache = np.load(
        cache_path,
        allow_pickle=True,
    )

    if "data" not in cache.files:
        raise ValueError(
            f"Missing data array in {cache_path}"
        )

    values = np.asarray(
        cache["data"]
    ).reshape(-1)

    return values.astype(float)


def read_spectra(
    path: Path,
    k_value: int,
) -> pd.DataFrame:

    table = pd.read_csv(
        path,
        sep="\t",
        index_col=0,
    )

    if table.shape[0] == k_value:
        result = table.copy()

    elif table.shape[1] == k_value:
        result = table.T.copy()

    else:
        raise ValueError(
            f"Cannot orient spectra file {path}: "
            f"shape={table.shape}, expected one dimension={k_value}"
        )

    result.index = [
        f"Program_{index + 1}"
        for index in range(result.shape[0])
    ]

    result.columns = result.columns.astype(str)

    return result


def match_programs(
    first: pd.DataFrame,
    second: pd.DataFrame,
    comparison_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:

    common_genes = sorted(
        set(first.columns)
        & set(second.columns)
    )

    if len(common_genes) < 500:
        raise ValueError(
            f"{comparison_name}: fewer than 500 shared genes."
        )

    first_matrix = first[
        common_genes
    ].to_numpy(dtype=float)

    second_matrix = second[
        common_genes
    ].to_numpy(dtype=float)

    correlation = np.corrcoef(
        np.vstack(
            [first_matrix, second_matrix]
        )
    )

    n_first = first_matrix.shape[0]

    correlation = correlation[
        :n_first,
        n_first:
    ]

    row_index, column_index = (
        linear_sum_assignment(
            -correlation
        )
    )

    match_table = pd.DataFrame(
        {
            "comparison": comparison_name,
            "first_program": [
                first.index[index]
                for index in row_index
            ],
            "second_program": [
                second.index[index]
                for index in column_index
            ],
            "pearson_correlation": [
                correlation[row, column]
                for row, column in zip(
                    row_index,
                    column_index,
                )
            ],
            "shared_genes": len(common_genes),
        }
    ).sort_values(
        "first_program"
    )

    correlation_table = pd.DataFrame(
        correlation,
        index=first.index,
        columns=second.index,
    )

    return match_table, correlation_table


def save_correlation_heatmap(
    correlation: pd.DataFrame,
    title: str,
    output_path: Path,
) -> None:

    figure, axis = plt.subplots(
        figsize=(6.5, 5.5)
    )

    image = axis.imshow(
        correlation.to_numpy(),
        vmin=-1,
        vmax=1,
        aspect="auto",
    )

    axis.set_xticks(
        range(correlation.shape[1])
    )

    axis.set_xticklabels(
        correlation.columns,
        rotation=45,
        ha="right",
    )

    axis.set_yticks(
        range(correlation.shape[0])
    )

    axis.set_yticklabels(
        correlation.index
    )

    axis.set_title(title)
    axis.set_xlabel("Second solution")
    axis.set_ylabel("First solution")

    figure.colorbar(
        image,
        ax=axis,
        label="Pearson correlation",
    )

    figure.tight_layout()

    figure.savefig(
        output_path,
        dpi=300,
    )

    plt.close(figure)


###############################################################################
# 3. Validate paths
###############################################################################

if not CNMF_EXE.is_file():
    raise FileNotFoundError(
        f"cNMF executable not found: {CNMF_EXE}"
    )

if not SOURCE_ROOT.is_dir():
    raise FileNotFoundError(
        f"32B2B source directory not found: {SOURCE_ROOT}"
    )

OUTPUT_ROOT.mkdir(
    parents=True,
    exist_ok=True,
)

LOG_ROOT.mkdir(
    parents=True,
    exist_ok=True,
)

REVIEW_ROOT.mkdir(
    parents=True,
    exist_ok=True,
)


###############################################################################
# 4. Run candidate and alternative filtered consensus
###############################################################################

removal_rows: list[dict[str, object]] = []
output_rows: list[dict[str, object]] = []

spectra_lookup: dict[
    tuple[str, str, float],
    pd.DataFrame
] = {}

overall_started = time.time()

for lineage, parameters in PLAN.items():

    k_value = parameters["k"]

    thresholds = [
        parameters["candidate_threshold"],
        *parameters["alternative_thresholds"],
    ]

    for sampling in SAMPLINGS:

        output_directory = (
            SOURCE_ROOT
            / lineage
            / sampling
        )

        run_name = (
            f"{lineage}_"
            f"{sampling}_"
            f"final_k{k_value}"
        )

        run_directory = (
            output_directory
            / run_name
        )

        if not run_directory.is_dir():
            raise FileNotFoundError(
                f"Run directory missing: {run_directory}"
            )

        density_values = load_local_density(
            run_directory,
            k_value,
        )

        for threshold in thresholds:

            tag = threshold_tag(
                threshold
            )

            log_path = (
                LOG_ROOT
                / (
                    f"{lineage}__"
                    f"{sampling}__"
                    f"consensus_dt_{tag}.log"
                )
            )

            run_command(
                label=(
                    f"{lineage}/{sampling}: "
                    f"consensus threshold={threshold}"
                ),
                arguments=[
                    "consensus",
                    "--output-dir",
                    str(output_directory),
                    "--name",
                    run_name,
                    "--components",
                    str(k_value),
                    "--local-density-threshold",
                    str(threshold),
                    "--show-clustering",
                ],
                log_path=log_path,
            )

            removed = int(
                np.sum(
                    density_values > threshold
                )
            )

            removal_rows.append(
                {
                    "lineage": lineage,
                    "sampling": sampling,
                    "k": k_value,
                    "threshold": threshold,
                    "total_spectra": len(
                        density_values
                    ),
                    "removed_spectra": removed,
                    "retained_spectra": (
                        len(density_values)
                        - removed
                    ),
                    "removed_fraction": (
                        removed
                        / len(density_values)
                    ),
                    "density_min": float(
                        np.min(density_values)
                    ),
                    "density_median": float(
                        np.median(density_values)
                    ),
                    "density_max": float(
                        np.max(density_values)
                    ),
                }
            )

            spectra_path = locate_one(
                run_directory,
                (
                    f"*.gene_spectra_score."
                    f"k_{k_value}.dt_{tag}*.txt"
                ),
            )

            spectra = read_spectra(
                spectra_path,
                k_value,
            )

            spectra_lookup[
                (
                    lineage,
                    sampling,
                    threshold,
                )
            ] = spectra

            review_directory = (
                REVIEW_ROOT
                / lineage
                / sampling
                / f"dt_{tag}"
            )

            review_directory.mkdir(
                parents=True,
                exist_ok=True,
            )

            review_patterns = [
                (
                    f"*.clustering."
                    f"k_{k_value}.dt_{tag}*"
                ),
                (
                    f"*.gene_spectra_score."
                    f"k_{k_value}.dt_{tag}*"
                ),
                (
                    f"*.top_genes."
                    f"k_{k_value}.dt_{tag}*"
                ),
            ]

            for pattern in review_patterns:

                for source_path in sorted(
                    run_directory.glob(
                        pattern
                    )
                ):

                    if not source_path.is_file():
                        continue

                    destination = (
                        review_directory
                        / source_path.name
                    )

                    shutil.copy2(
                        source_path,
                        destination,
                    )

                    output_rows.append(
                        {
                            "lineage": lineage,
                            "sampling": sampling,
                            "k": k_value,
                            "threshold": threshold,
                            "file_name": (
                                source_path.name
                            ),
                            "source_path": str(
                                source_path
                            ),
                            "size_bytes": (
                                source_path.stat().st_size
                            ),
                            "sha256": sha256_file(
                                source_path
                            ),
                        }
                    )

            shutil.copy2(
                log_path,
                review_directory
                / log_path.name,
            )


###############################################################################
# 5. Threshold-removal audit
###############################################################################

removal_table = pd.DataFrame(
    removal_rows
).sort_values(
    [
        "lineage",
        "sampling",
        "threshold",
    ]
)

removal_path = (
    OUTPUT_ROOT
    / "32B2C_THRESHOLD_REMOVAL_AUDIT.csv"
)

removal_table.to_csv(
    removal_path,
    index=False,
)

shutil.copy2(
    removal_path,
    REVIEW_ROOT / removal_path.name,
)


###############################################################################
# 6. Cross-sampling reproducibility at candidate threshold
###############################################################################

program_match_tables: list[pd.DataFrame] = []

for lineage, parameters in PLAN.items():

    candidate = parameters[
        "candidate_threshold"
    ]

    primary = spectra_lookup[
        (
            lineage,
            "primary",
            candidate,
        )
    ]

    sensitivity = spectra_lookup[
        (
            lineage,
            "sensitivity",
            candidate,
        )
    ]

    comparison_name = (
        f"{lineage}: "
        f"primary_vs_sensitivity_"
        f"dt_{threshold_tag(candidate)}"
    )

    matches, correlation = match_programs(
        primary,
        sensitivity,
        comparison_name,
    )

    matches["comparison_type"] = (
        "cross_sampling_candidate_threshold"
    )

    program_match_tables.append(
        matches
    )

    correlation_path = (
        REVIEW_ROOT
        / lineage
        / (
            "candidate_cross_sampling_"
            "correlation.csv"
        )
    )

    correlation.to_csv(
        correlation_path
    )

    save_correlation_heatmap(
        correlation,
        comparison_name,
        REVIEW_ROOT
        / lineage
        / (
            "candidate_cross_sampling_"
            "correlation.png"
        ),
    )


###############################################################################
# 7. Within-sampling threshold sensitivity
###############################################################################

for lineage, parameters in PLAN.items():

    alternatives = parameters[
        "alternative_thresholds"
    ]

    if not alternatives:
        continue

    candidate = parameters[
        "candidate_threshold"
    ]

    for alternative in alternatives:

        for sampling in SAMPLINGS:

            candidate_spectra = (
                spectra_lookup[
                    (
                        lineage,
                        sampling,
                        candidate,
                    )
                ]
            )

            alternative_spectra = (
                spectra_lookup[
                    (
                        lineage,
                        sampling,
                        alternative,
                    )
                ]
            )

            comparison_name = (
                f"{lineage}/{sampling}: "
                f"dt_{threshold_tag(candidate)}"
                f"_vs_"
                f"dt_{threshold_tag(alternative)}"
            )

            matches, correlation = (
                match_programs(
                    candidate_spectra,
                    alternative_spectra,
                    comparison_name,
                )
            )

            matches["comparison_type"] = (
                "within_sampling_threshold_sensitivity"
            )

            program_match_tables.append(
                matches
            )

            output_prefix = (
                f"{sampling}_"
                f"dt_{threshold_tag(candidate)}"
                f"_vs_"
                f"dt_{threshold_tag(alternative)}"
            )

            correlation.to_csv(
                REVIEW_ROOT
                / lineage
                / f"{output_prefix}.csv"
            )

            save_correlation_heatmap(
                correlation,
                comparison_name,
                REVIEW_ROOT
                / lineage
                / f"{output_prefix}.png",
            )


###############################################################################
# 8. Program-match table
###############################################################################

program_match_table = pd.concat(
    program_match_tables,
    ignore_index=True,
)

program_match_path = (
    OUTPUT_ROOT
    / "32B2C_PROGRAM_MATCH_AUDIT.csv"
)

program_match_table.to_csv(
    program_match_path,
    index=False,
)

shutil.copy2(
    program_match_path,
    REVIEW_ROOT / program_match_path.name,
)


###############################################################################
# 9. Summary metrics
###############################################################################

summary_rows = []

for (
    comparison,
    comparison_type,
), subset in program_match_table.groupby(
    [
        "comparison",
        "comparison_type",
    ]
):

    summary_rows.append(
        {
            "comparison": comparison,
            "comparison_type": comparison_type,
            "programs": len(subset),
            "minimum_matched_correlation": (
                subset[
                    "pearson_correlation"
                ].min()
            ),
            "median_matched_correlation": (
                subset[
                    "pearson_correlation"
                ].median()
            ),
            "mean_matched_correlation": (
                subset[
                    "pearson_correlation"
                ].mean()
            ),
            "programs_correlation_ge_0_70": int(
                np.sum(
                    subset[
                        "pearson_correlation"
                    ] >= 0.70
                )
            ),
            "programs_correlation_ge_0_50": int(
                np.sum(
                    subset[
                        "pearson_correlation"
                    ] >= 0.50
                )
            ),
        }
    )

summary_table = pd.DataFrame(
    summary_rows
)

summary_path = (
    OUTPUT_ROOT
    / "32B2C_PROGRAM_MATCH_SUMMARY.csv"
)

summary_table.to_csv(
    summary_path,
    index=False,
)

shutil.copy2(
    summary_path,
    REVIEW_ROOT / summary_path.name,
)


###############################################################################
# 10. Output manifest
###############################################################################

output_manifest = pd.DataFrame(
    output_rows
)

manifest_path = (
    OUTPUT_ROOT
    / "32B2C_FILTERED_OUTPUT_MANIFEST.csv"
)

output_manifest.to_csv(
    manifest_path,
    index=False,
)

shutil.copy2(
    manifest_path,
    REVIEW_ROOT / manifest_path.name,
)


###############################################################################
# 11. Plan lock
###############################################################################

plan_lock = {
    "stage": "32B2C",
    "created_at": now_string(),
    "candidate_thresholds": {
        lineage: parameters[
            "candidate_threshold"
        ]
        for lineage, parameters
        in PLAN.items()
    },
    "alternative_thresholds": {
        lineage: parameters[
            "alternative_thresholds"
        ]
        for lineage, parameters
        in PLAN.items()
    },
    "disease_labels_used": False,
    "senescence_signatures_used": False,
    "factorization_rerun": False,
    "purpose": (
        "Filtered consensus, threshold sensitivity, "
        "and primary-versus-sensitivity program matching."
    ),
}

lock_path = (
    OUTPUT_ROOT
    / "32B2C_THRESHOLD_PLAN_LOCK.json"
)

with lock_path.open(
    "w",
    encoding="utf-8",
) as handle:

    json.dump(
        plan_lock,
        handle,
        indent=2,
        ensure_ascii=False,
    )

shutil.copy2(
    lock_path,
    REVIEW_ROOT / lock_path.name,
)


###############################################################################
# 12. Completion report
###############################################################################

elapsed = time.time() - overall_started

completion_text = f"""\
32B-2C filtered consensus threshold audit completed successfully.

Completed:
{now_string()}

Elapsed seconds:
{elapsed:.2f}

Candidate thresholds:
- Astrocytes: 0.02
- Microglia/macrophage: 0.01
- OPC: 0.01
- Oligodendrocytes: 0.01

Targeted alternative thresholds:
- Astrocytes: 0.05
- OPC: 0.05

Factorization rerun:
No.

Disease labels used:
No.

Senescence signatures used:
No.

Next:
Review cross-sampling program correlations and threshold sensitivity
before freezing the final filtered consensus.
"""

completion_path = (
    OUTPUT_ROOT
    / "32B2C_COMPLETED_SUCCESSFULLY.txt"
)

completion_path.write_text(
    completion_text,
    encoding="utf-8",
)

shutil.copy2(
    completion_path,
    REVIEW_ROOT / completion_path.name,
)


###############################################################################
# 13. ZIP and SHA-256
###############################################################################

review_zip = (
    OUTPUT_ROOT
    / "32B2C_FILTERED_CONSENSUS_REVIEW.zip"
)

if review_zip.exists():
    review_zip.unlink()

with zipfile.ZipFile(
    review_zip,
    mode="w",
    compression=zipfile.ZIP_DEFLATED,
    compresslevel=9,
) as archive:

    for path in sorted(
        REVIEW_ROOT.rglob("*")
    ):
        if path.is_file():
            archive.write(
                path,
                arcname=path.relative_to(
                    REVIEW_ROOT
                ),
            )

review_hash = sha256_file(
    review_zip
)

hash_path = Path(
    str(review_zip)
    + ".sha256.txt"
)

hash_path.write_text(
    f"{review_hash}  {review_zip.name}\n",
    encoding="ascii",
)


###############################################################################
# 14. Final output
###############################################################################

print()
print("=" * 72)
print("32B-2C completed successfully")
print("=" * 72)
print(f"Elapsed: {elapsed:.2f} seconds")
print()
print("Review ZIP:")
print(review_zip)
print()
print("SHA256:")
print(review_hash)
print("=" * 72)