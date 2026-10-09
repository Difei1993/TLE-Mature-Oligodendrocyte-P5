###############################################################################
# 32B-2B
# Final fixed-K cNMF:
#   - 100 NMF replicates
#   - primary and independent-sensitivity samples
#   - unfiltered consensus diagnostic at density threshold 2.0
#
# This stage does NOT:
#   - use disease labels
#   - use senescence signatures
#   - select the final density-filter threshold
###############################################################################

from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import zipfile
from datetime import datetime
from pathlib import Path


###############################################################################
# 0. Paths and frozen parameters
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

PRIMARY_INPUT_ROOT = (
    PROJECT
    / "04_tables"
    / "32B_CNMf_WORK_20260728_215801"
    / "_h5ad_inputs"
)

SENSITIVITY_INPUT_ROOT = (
    PROJECT
    / "04_tables"
    / "32B2A_CNMf_SENSITIVITY_20260728_215801"
    / "_h5ad_inputs"
)

OUTPUT_ROOT = (
    PROJECT
    / "04_tables"
    / "32B2B_FINAL_CNMf_20260728_215801"
)

LOG_ROOT = OUTPUT_ROOT / "_logs"

REVIEW_ROOT = OUTPUT_ROOT / "_UNFILTERED_CONSENSUS_REVIEW"

FINAL_K = {
    "Astrocytes": 7,
    "Microglia_macrophage": 6,
    "OPC": 6,
    "Oligodendrocytes": 6,
}

INPUTS = {
    "primary": {
        "Astrocytes": (
            PRIMARY_INPUT_ROOT
            / "Astrocytes_primary.h5ad"
        ),
        "Microglia_macrophage": (
            PRIMARY_INPUT_ROOT
            / "Microglia_macrophage_primary.h5ad"
        ),
        "OPC": (
            PRIMARY_INPUT_ROOT
            / "OPC_primary.h5ad"
        ),
        "Oligodendrocytes": (
            PRIMARY_INPUT_ROOT
            / "Oligodendrocytes_primary.h5ad"
        ),
    },
    "sensitivity": {
        "Astrocytes": (
            SENSITIVITY_INPUT_ROOT
            / "Astrocytes_sensitivity.h5ad"
        ),
        "Microglia_macrophage": (
            SENSITIVITY_INPUT_ROOT
            / "Microglia_macrophage_sensitivity.h5ad"
        ),
        "OPC": (
            SENSITIVITY_INPUT_ROOT
            / "OPC_sensitivity.h5ad"
        ),
        "Oligodendrocytes": (
            SENSITIVITY_INPUT_ROOT
            / "Oligodendrocytes_sensitivity.h5ad"
        ),
    },
}

N_ITER = 100
SEED = 20260730
NUM_GENES = 2000
MAX_NMF_ITER = 2000

# At 2.0, no meaningful local-density filtering is applied.
DIAGNOSTIC_DENSITY_THRESHOLD = "2.0"


###############################################################################
# 1. Utilities
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


def run_command(
    *,
    label: str,
    arguments: list[str],
    log_path: Path,
) -> None:
    """Run one cNMF command and preserve its complete output."""

    log_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    started = time.time()

    print()
    print("=" * 72)
    print(label)
    print(f"Started: {now_string()}")
    print("=" * 72)

    command = [
        str(CNMF_EXE),
        *arguments,
    ]

    environment = os.environ.copy()

    # Avoid uncontrolled BLAS oversubscription and improve reproducibility.
    environment["OMP_NUM_THREADS"] = "1"
    environment["MKL_NUM_THREADS"] = "1"
    environment["OPENBLAS_NUM_THREADS"] = "1"
    environment["NUMEXPR_NUM_THREADS"] = "1"

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
        log_handle.flush()

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
            f"See log: {log_path}"
        )

    print(
        f"Completed: {label}\n"
        f"Elapsed: {elapsed / 60:.2f} minutes"
    )


def find_consensus_outputs(
    run_directory: Path,
    k_value: int,
) -> list[Path]:

    patterns = [
        f"*clustering.k_{k_value}.dt_2_0*",
        f"*gene_spectra_score.k_{k_value}.dt_2_0*",
        f"*gene_spectra_tpm.k_{k_value}.dt_2_0*",
        f"*usages.k_{k_value}.dt_2_0*",
        f"*top_genes.k_{k_value}.dt_2_0*",
    ]

    files: list[Path] = []

    for pattern in patterns:
        files.extend(
            path
            for path in run_directory.rglob(pattern)
            if path.is_file()
        )

    # Deduplicate while preserving path order.
    unique_files = list(
        dict.fromkeys(files)
    )

    return unique_files


###############################################################################
# 2. Validate environment and inputs
###############################################################################

if not CNMF_EXE.is_file():
    raise FileNotFoundError(
        f"cNMF executable not found: {CNMF_EXE}"
    )

OUTPUT_ROOT.mkdir(
    parents=True,
    exist_ok=True,
)

LOG_ROOT.mkdir(
    parents=True,
    exist_ok=True,
)

for sampling_name, lineage_inputs in INPUTS.items():
    for lineage_name, input_path in lineage_inputs.items():

        if not input_path.is_file():
            raise FileNotFoundError(
                f"Missing {sampling_name} H5AD for "
                f"{lineage_name}: {input_path}"
            )

        if input_path.stat().st_size <= 0:
            raise ValueError(
                f"Empty H5AD input: {input_path}"
            )

print("All eight H5AD inputs passed validation.")
print(f"cNMF executable: {CNMF_EXE}")


###############################################################################
# 3. Frozen K lock
###############################################################################

k_lock = {
    "stage": "32B2B",
    "created_at": now_string(),
    "final_k": FINAL_K,
    "selection_basis": {
        "primary_sampling": "32B1",
        "independent_sampling": "32B2A",
        "rule": (
            "Prefer reproducible high stability and parsimony "
            "when prediction-error improvements at larger K "
            "are small."
        ),
    },
    "nmf_replicates": N_ITER,
    "seed": SEED,
    "num_high_variance_genes": NUM_GENES,
    "diagnostic_consensus_density_threshold": 2.0,
    "disease_labels_used_for_discovery": False,
    "senescence_signatures_used_for_discovery": False,
}

with (
    OUTPUT_ROOT
    / "32B2B_FINAL_K_LOCK.json"
).open(
    "w",
    encoding="utf-8",
) as handle:
    json.dump(
        k_lock,
        handle,
        indent=2,
        ensure_ascii=False,
    )


###############################################################################
# 4. Run fixed-K primary and sensitivity analyses
###############################################################################

result_rows: list[dict[str, object]] = []

overall_started = time.time()

for lineage_name, k_value in FINAL_K.items():

    for sampling_name in (
        "primary",
        "sensitivity",
    ):

        input_path = INPUTS[
            sampling_name
        ][lineage_name]

        output_directory = (
            OUTPUT_ROOT
            / lineage_name
            / sampling_name
        )

        output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        run_name = (
            f"{lineage_name}_"
            f"{sampling_name}_"
            f"final_k{k_value}"
        )

        run_directory = (
            output_directory
            / run_name
        )

        marker_directory = (
            output_directory
            / "_markers"
        )

        marker_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        prepare_marker = (
            marker_directory
            / "PREPARE_COMPLETED.txt"
        )

        factorize_marker = (
            marker_directory
            / "FACTORIZE_COMPLETED.txt"
        )

        combine_marker = (
            marker_directory
            / "COMBINE_COMPLETED.txt"
        )

        consensus_marker = (
            marker_directory
            / "UNFILTERED_CONSENSUS_COMPLETED.txt"
        )

        log_prefix = (
            f"{lineage_name}__"
            f"{sampling_name}"
        )

        #######################################################################
        # prepare
        #######################################################################

        if not prepare_marker.exists():

            run_command(
                label=(
                    f"{lineage_name} / "
                    f"{sampling_name} : prepare"
                ),
                arguments=[
                    "prepare",
                    "--output-dir",
                    str(output_directory),
                    "--name",
                    run_name,
                    "-c",
                    str(input_path),
                    "-k",
                    str(k_value),
                    "--n-iter",
                    str(N_ITER),
                    "--total-workers",
                    "1",
                    "--seed",
                    str(SEED),
                    "--numgenes",
                    str(NUM_GENES),
                    "--max-nmf-iter",
                    str(MAX_NMF_ITER),
                ],
                log_path=(
                    LOG_ROOT
                    / f"{log_prefix}__01_prepare.log"
                ),
            )

            prepare_marker.write_text(
                f"Completed: {now_string()}\n",
                encoding="utf-8",
            )

        else:
            print(
                f"{lineage_name}/{sampling_name}: "
                "prepare already completed."
            )

        #######################################################################
        # factorize
        #######################################################################

        if not factorize_marker.exists():

            run_command(
                label=(
                    f"{lineage_name} / "
                    f"{sampling_name} : factorize"
                ),
                arguments=[
                    "factorize",
                    "--output-dir",
                    str(output_directory),
                    "--name",
                    run_name,
                    "--worker-index",
                    "0",
                    "--total-workers",
                    "1",
                    "--skip-completed-runs",
                ],
                log_path=(
                    LOG_ROOT
                    / f"{log_prefix}__02_factorize.log"
                ),
            )

            factorize_marker.write_text(
                f"Completed: {now_string()}\n",
                encoding="utf-8",
            )

        else:
            print(
                f"{lineage_name}/{sampling_name}: "
                "factorize already completed."
            )

        #######################################################################
        # combine
        #######################################################################

        if not combine_marker.exists():

            run_command(
                label=(
                    f"{lineage_name} / "
                    f"{sampling_name} : combine"
                ),
                arguments=[
                    "combine",
                    "--output-dir",
                    str(output_directory),
                    "--name",
                    run_name,
                ],
                log_path=(
                    LOG_ROOT
                    / f"{log_prefix}__03_combine.log"
                ),
            )

            combine_marker.write_text(
                f"Completed: {now_string()}\n",
                encoding="utf-8",
            )

        else:
            print(
                f"{lineage_name}/{sampling_name}: "
                "combine already completed."
            )

        #######################################################################
        # Diagnostic unfiltered consensus
        #######################################################################

        if not consensus_marker.exists():

            run_command(
                label=(
                    f"{lineage_name} / "
                    f"{sampling_name} : "
                    "unfiltered consensus"
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
                    DIAGNOSTIC_DENSITY_THRESHOLD,
                    "--show-clustering",
                ],
                log_path=(
                    LOG_ROOT
                    / f"{log_prefix}__04_consensus_dt2.log"
                ),
            )

            consensus_marker.write_text(
                f"Completed: {now_string()}\n",
                encoding="utf-8",
            )

        else:
            print(
                f"{lineage_name}/{sampling_name}: "
                "unfiltered consensus already completed."
            )

        #######################################################################
        # Validate consensus output
        #######################################################################

        consensus_outputs = (
            find_consensus_outputs(
                run_directory,
                k_value,
            )
        )

        clustering_outputs = [
            path
            for path in consensus_outputs
            if "clustering" in path.name
        ]

        if not clustering_outputs:
            raise RuntimeError(
                f"No unfiltered consensus clustergram "
                f"found for {lineage_name}/{sampling_name}. "
                f"Run directory: {run_directory}"
            )

        for output_path in consensus_outputs:
            result_rows.append(
                {
                    "lineage": lineage_name,
                    "sampling": sampling_name,
                    "k": k_value,
                    "output_file": str(output_path),
                    "size_bytes": (
                        output_path.stat().st_size
                    ),
                    "sha256": sha256_file(
                        output_path
                    ),
                }
            )


###############################################################################
# 5. Result manifest
###############################################################################

manifest_path = (
    OUTPUT_ROOT
    / "32B2B_CONSENSUS_OUTPUT_MANIFEST.csv"
)

with manifest_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "lineage",
            "sampling",
            "k",
            "output_file",
            "size_bytes",
            "sha256",
        ],
    )

    writer.writeheader()
    writer.writerows(result_rows)


###############################################################################
# 6. Build small review package
###############################################################################

if REVIEW_ROOT.exists():
    shutil.rmtree(REVIEW_ROOT)

REVIEW_ROOT.mkdir(
    parents=True,
    exist_ok=True,
)

shutil.copy2(
    OUTPUT_ROOT / "32B2B_FINAL_K_LOCK.json",
    REVIEW_ROOT / "32B2B_FINAL_K_LOCK.json",
)

shutil.copy2(
    manifest_path,
    REVIEW_ROOT / manifest_path.name,
)

for lineage_name, k_value in FINAL_K.items():

    for sampling_name in (
        "primary",
        "sensitivity",
    ):

        source_run_directory = (
            OUTPUT_ROOT
            / lineage_name
            / sampling_name
        )

        review_directory = (
            REVIEW_ROOT
            / lineage_name
            / sampling_name
        )

        review_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        diagnostic_files = [
            path
            for path in source_run_directory.rglob("*")
            if (
                path.is_file()
                and (
                    "clustering" in path.name
                    or "local_density" in path.name
                )
            )
        ]

        if not diagnostic_files:
            raise RuntimeError(
                "No consensus diagnostic files found for "
                f"{lineage_name}/{sampling_name}."
            )

        for source_path in diagnostic_files:

            destination = (
                review_directory
                / source_path.name
            )

            shutil.copy2(
                source_path,
                destination,
            )

        consensus_log = (
            LOG_ROOT
            / (
                f"{lineage_name}__"
                f"{sampling_name}__"
                "04_consensus_dt2.log"
            )
        )

        if consensus_log.exists():
            shutil.copy2(
                consensus_log,
                review_directory
                / consensus_log.name,
            )


###############################################################################
# 7. Completion summary
###############################################################################

overall_elapsed = time.time() - overall_started

summary_text = f"""\
32B-2B final fixed-K cNMF completed successfully.

Completed:
{now_string()}

Elapsed seconds:
{overall_elapsed:.2f}

Frozen K:
- Astrocytes: 7
- Microglia/macrophage: 6
- OPC: 6
- Oligodendrocytes: 6

NMF replicates:
100 per lineage per sampling.

Sampling datasets:
- Primary balanced sampling
- Independent sensitivity balanced sampling

Consensus diagnostic:
- local-density-threshold = 2.0
- show-clustering enabled
- no final outlier-filter threshold selected yet

Disease labels used:
No.

Senescence signatures used:
No.

Next stage:
Review all eight consensus clustergrams and local-density distributions,
then freeze the final density threshold before producing the filtered
consensus results.
"""

completion_path = (
    OUTPUT_ROOT
    / "32B2B_COMPLETED_SUCCESSFULLY.txt"
)

completion_path.write_text(
    summary_text,
    encoding="utf-8",
)

shutil.copy2(
    completion_path,
    REVIEW_ROOT / completion_path.name,
)


###############################################################################
# 8. ZIP and SHA-256
###############################################################################

review_zip = (
    OUTPUT_ROOT
    / "32B2B_UNFILTERED_CONSENSUS_REVIEW.zip"
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
# 9. Final console output
###############################################################################

print()
print("=" * 72)
print("32B-2B completed successfully")
print("=" * 72)
print(f"Elapsed: {overall_elapsed / 60:.2f} minutes")
print()
print("Review ZIP:")
print(review_zip)
print()
print("SHA256:")
print(review_hash)
print()
print(
    "Final density threshold has NOT yet been selected."
)
print("=" * 72)