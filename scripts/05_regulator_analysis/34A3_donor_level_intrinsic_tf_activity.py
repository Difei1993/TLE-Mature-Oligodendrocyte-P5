from __future__ import annotations

import hashlib
import itertools
import json
import os
import sys
import zipfile
from datetime import datetime
from pathlib import Path

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr


PROJECT = Path(
    os.environ.get(
        "TLE_PROJECT",
        __import__("os").environ["P5_WORK_ROOT"],
    )
)
RELEASE_ROOT = PROJECT / "99_frozen_releases"
CNMF_WORK = PROJECT / "04_tables" / "32B_CNMf_WORK_20260728_215801"
STAMP = datetime.now().strftime("%Y%m%d_%H%M%S")
OUT = RELEASE_ROOT / f"34A3_DONOR_LEVEL_INTRINSIC_TF_ACTIVITY_{STAMP}"

H5AD_NAMES = {
    "Astrocytes": "Astrocytes_primary.h5ad",
    "Microglia_macrophage": "Microglia_macrophage_primary.h5ad",
    "Oligodendrocytes": "Oligodendrocytes_primary.h5ad",
}
TARGET_SUM = 10000.0
TMIN = 5
PRIMARY_NETWORK = "PRIMARY_ORTHOGONAL_WITH_COMPLEXES"
SUPPORTIVE_NETWORK = "SUPPORTIVE_FULL_WITH_COMPLEXES"
SENSITIVITY_NETWORK = "SENSITIVITY_ORTHOGONAL_NO_COMPLEXES"
NETWORK_VARIANTS = [PRIMARY_NETWORK, SUPPORTIVE_NETWORK, SENSITIVITY_NETWORK]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def latest_release(prefix: str, required: list[str]) -> Path:
    candidates = sorted(
        [path for path in RELEASE_ROOT.glob(prefix + "*") if path.is_dir()],
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    for root in candidates:
        complete = all(
            len([path for path in root.rglob(name) if path.is_file()]) == 1
            for name in required
        )
        if complete:
            return root.resolve()
    raise FileNotFoundError(f"No complete release found: {prefix}; required={required}")


def one_file(root: Path, name: str) -> Path:
    hits = [path for path in root.rglob(name) if path.is_file()]
    if len(hits) != 1:
        raise RuntimeError(
            f"File not uniquely located: {name}; found={len(hits)}\n"
            + "\n".join(map(str, hits))
        )
    return hits[0].resolve()


def find_h5ad(name: str) -> Path:
    hits = [path for path in CNMF_WORK.rglob(name) if path.is_file()]
    if len(hits) != 1:
        raise RuntimeError(f"H5AD not uniquely located: {name}; found={len(hits)}")
    return hits[0].resolve()


def norm_gene(value) -> str:
    value = "" if value is None else str(value).strip().upper()
    if value in {"", "NA", "N/A", "NONE", "NULL", "---"}:
        return ""
    return value


def normalize_condition(value) -> str:
    value = str(value).strip().lower()
    if value in {"normal", "control", "ctrl", "con"}:
        return "Normal"
    if value in {"tle", "epilepsy", "case", "disease"}:
        return "TLE"
    raise ValueError(f"Unrecognized condition: {value!r}")


def bh_adjust(p_values: pd.Series) -> pd.Series:
    values = pd.to_numeric(p_values, errors="coerce").to_numpy(dtype=float)
    output = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return pd.Series(output, index=p_values.index)
    finite_values = values[finite]
    order = np.argsort(finite_values)
    ranked = finite_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    output[np.flatnonzero(finite)] = restored
    return pd.Series(output, index=p_values.index)


def residualize_condition(values: np.ndarray, condition: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float).copy()
    condition = np.asarray(condition).astype(str)
    for group in pd.unique(condition):
        mask = condition == group
        values[mask] -= values[mask].mean()
    return values


def pearson_from_centered(first: np.ndarray, second: np.ndarray) -> float:
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    denominator = float(np.sqrt(np.dot(first, first) * np.dot(second, second)))
    if not np.isfinite(denominator) or denominator <= 1e-15:
        return np.nan
    return float(np.dot(first, second) / denominator)


def partial_r_condition(
    first: np.ndarray,
    second: np.ndarray,
    condition: np.ndarray,
) -> float:
    return pearson_from_centered(
        residualize_condition(first, condition),
        residualize_condition(second, condition),
    )


def within_condition_permutation_indices(condition: np.ndarray) -> np.ndarray:
    condition = np.asarray(condition).astype(str)
    groups = [np.flatnonzero(condition == group) for group in sorted(pd.unique(condition))]
    if sorted(map(len, groups)) != [4, 4]:
        raise RuntimeError(f"Expected two groups of four donors, found {[len(g) for g in groups]}")
    rows = []
    base = np.arange(len(condition), dtype=int)
    for first_perm in itertools.permutations(groups[0].tolist()):
        for second_perm in itertools.permutations(groups[1].tolist()):
            row = base.copy()
            row[groups[0]] = np.asarray(first_perm, dtype=int)
            row[groups[1]] = np.asarray(second_perm, dtype=int)
            rows.append(row)
    result = np.vstack(rows)
    if result.shape != (576, 8):
        raise RuntimeError(f"Unexpected conditional permutation matrix shape: {result.shape}")
    return result


def exact_partial_permutation_p(
    activity: np.ndarray,
    program_score: np.ndarray,
    condition: np.ndarray,
    permutation_indices: np.ndarray,
) -> tuple[float, float]:
    activity_residual = residualize_condition(activity, condition)
    score_residual = residualize_condition(program_score, condition)
    observed = pearson_from_centered(activity_residual, score_residual)
    if not np.isfinite(observed):
        return np.nan, np.nan
    permuted_scores = score_residual[permutation_indices]
    activity_ss = float(np.dot(activity_residual, activity_residual))
    score_ss = np.sum(permuted_scores * permuted_scores, axis=1)
    denominators = np.sqrt(activity_ss * score_ss)
    numerators = permuted_scores @ activity_residual
    permuted_r = np.divide(
        numerators,
        denominators,
        out=np.full(numerators.shape, np.nan, dtype=float),
        where=denominators > 1e-15,
    )
    finite = np.isfinite(permuted_r)
    if not finite.any():
        return observed, np.nan
    p_value = float(
        np.mean(np.abs(permuted_r[finite]) >= abs(observed) - 1e-12)
    )
    return observed, p_value


def lodo_partial_summary(
    activity: np.ndarray,
    program_score: np.ndarray,
    condition: np.ndarray,
    full_r: float,
) -> dict:
    values = []
    for omitted in range(len(activity)):
        keep = np.arange(len(activity)) != omitted
        values.append(
            partial_r_condition(activity[keep], program_score[keep], condition[keep])
        )
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    if not finite.any() or not np.isfinite(full_r):
        return {
            "lodo_finite": int(finite.sum()),
            "lodo_same_sign_fraction": np.nan,
            "lodo_median_partial_r": np.nan,
            "lodo_minimum_partial_r": np.nan,
            "lodo_maximum_partial_r": np.nan,
        }
    same_sign = np.sign(values[finite]) == np.sign(full_r)
    return {
        "lodo_finite": int(finite.sum()),
        "lodo_same_sign_fraction": float(np.mean(same_sign)),
        "lodo_median_partial_r": float(np.median(values[finite])),
        "lodo_minimum_partial_r": float(np.min(values[finite])),
        "lodo_maximum_partial_r": float(np.max(values[finite])),
    }


def hedges_g(values: np.ndarray, condition: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    condition = np.asarray(condition).astype(str)
    normal = values[condition == "Normal"]
    tle = values[condition == "TLE"]
    if len(normal) != 4 or len(tle) != 4:
        return np.nan
    pooled_variance = (
        ((len(normal) - 1) * np.var(normal, ddof=1))
        + ((len(tle) - 1) * np.var(tle, ddof=1))
    ) / (len(normal) + len(tle) - 2)
    if not np.isfinite(pooled_variance) or pooled_variance <= 1e-15:
        return np.nan
    cohen_d = (tle.mean() - normal.mean()) / np.sqrt(pooled_variance)
    correction = 1.0 - 3.0 / (4.0 * (len(normal) + len(tle)) - 9.0)
    return float(correction * cohen_d)


def exact_group_difference_p(values: np.ndarray, condition: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    condition = np.asarray(condition).astype(str)
    tle_positions = np.flatnonzero(condition == "TLE")
    if len(tle_positions) != 4 or len(values) != 8:
        return np.nan, np.nan
    observed = float(values[tle_positions].mean() - values[condition == "Normal"].mean())
    differences = []
    all_indices = np.arange(8)
    for candidate_tle in itertools.combinations(all_indices, 4):
        candidate_tle = np.asarray(candidate_tle, dtype=int)
        candidate_normal = np.setdiff1d(all_indices, candidate_tle)
        differences.append(
            float(values[candidate_tle].mean() - values[candidate_normal].mean())
        )
    differences = np.asarray(differences, dtype=float)
    return observed, float(np.mean(np.abs(differences) >= abs(observed) - 1e-12))


def prune_network(
    network: pd.DataFrame,
    expression_genes: set[str],
    excluded_genes: set[str] | None = None,
) -> pd.DataFrame:
    required = {"source", "target", "weight"}
    missing = required - set(network.columns)
    if missing:
        raise RuntimeError(f"Network missing columns: {sorted(missing)}")
    result = network[["source", "target", "weight"]].copy()
    result["source"] = result["source"].astype(str).str.strip()
    result["target"] = result["target"].map(norm_gene)
    result["weight"] = pd.to_numeric(result["weight"], errors="coerce")
    result = result.loc[
        result["source"].ne("")
        & result["target"].isin(expression_genes)
        & np.isfinite(result["weight"])
    ].copy()
    if excluded_genes:
        result = result.loc[~result["target"].isin(excluded_genes)].copy()
    result = result.drop_duplicates(["source", "target"], keep="first")
    target_counts = result.groupby("source", observed=True)["target"].nunique()
    eligible = set(target_counts.loc[target_counts >= TMIN].index.astype(str))
    result = result.loc[result["source"].isin(eligible)].reset_index(drop=True)
    if result.empty:
        raise RuntimeError("Network became empty after expression matching and pruning")
    return result


def build_pseudobulk(
    h5ad_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    adata = ad.read_h5ad(h5ad_path)
    required_obs = {"donor_id", "condition"}
    missing = required_obs - set(adata.obs.columns)
    if missing:
        raise RuntimeError(f"H5AD missing metadata columns: {sorted(missing)}")
    donor = adata.obs["donor_id"].astype(str).to_numpy()
    condition = np.asarray([normalize_condition(v) for v in adata.obs["condition"]], dtype=str)
    donor_order = sorted(pd.unique(donor))
    if len(donor_order) != 8:
        raise RuntimeError(f"Expected eight donors, found {len(donor_order)}")
    count_rows = []
    condition_rows = []
    cell_rows = []
    for donor_id in donor_order:
        mask = donor == donor_id
        donor_conditions = pd.unique(condition[mask])
        if len(donor_conditions) != 1:
            raise RuntimeError(f"Donor {donor_id} has multiple conditions")
        subset = adata.X[mask]
        if sparse.issparse(subset):
            row = np.asarray(subset.sum(axis=0)).ravel()
        else:
            row = np.asarray(subset, dtype=float).sum(axis=0)
        count_rows.append(row)
        condition_rows.append(str(donor_conditions[0]))
        cell_rows.append(int(mask.sum()))
    raw_counts = np.vstack(count_rows).astype(float)
    genes = np.asarray([norm_gene(v) for v in adata.var_names], dtype=str)
    valid = genes != ""
    counts_frame = pd.DataFrame(raw_counts[:, valid], index=donor_order, columns=genes[valid])
    duplicate_gene_columns = int(counts_frame.columns.duplicated(keep=False).sum())
    if duplicate_gene_columns:
        counts_frame = counts_frame.T.groupby(level=0, sort=True).sum().T
    counts_frame = counts_frame.loc[:, counts_frame.sum(axis=0) > 0]
    library_size = counts_frame.sum(axis=1)
    if (library_size <= 0).any():
        raise RuntimeError("At least one donor pseudobulk has zero library size")
    lognorm = np.log1p(counts_frame.div(library_size, axis=0) * TARGET_SUM)
    metadata = pd.DataFrame(
        {
            "donor_id": donor_order,
            "condition": condition_rows,
            "n_cells": cell_rows,
            "library_size": library_size.to_numpy(dtype=float),
            "detected_genes": (counts_frame > 0).sum(axis=1).to_numpy(dtype=int),
        }
    )
    audit = {
        "h5ad_path": str(h5ad_path),
        "h5ad_sha256": sha256_file(h5ad_path),
        "n_cells": int(adata.n_obs),
        "original_genes": int(adata.n_vars),
        "collapsed_nonzero_genes": int(counts_frame.shape[1]),
        "duplicate_gene_columns_after_normalization": duplicate_gene_columns,
        "minimum_library_size": float(library_size.min()),
        "maximum_library_size": float(library_size.max()),
        "target_sum": TARGET_SUM,
        "transformation": "sum raw counts per donor; normalize_total to 1e4; log1p",
    }
    del adata
    return counts_frame, lognorm, metadata, audit


def run_ulm(expression: pd.DataFrame, network: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    estimate, padj = dc.mt.ulm(
        data=expression,
        net=network[["source", "target", "weight"]],
        tmin=TMIN,
        tval=True,
        empty=True,
        verbose=False,
    )
    estimate = pd.DataFrame(estimate).copy()
    padj = pd.DataFrame(padj).copy()
    if estimate.shape[0] != expression.shape[0]:
        raise RuntimeError(f"ULM returned unexpected observations: {estimate.shape}")
    estimate.index = expression.index.astype(str)
    padj.index = expression.index.astype(str)
    if not np.isfinite(estimate.to_numpy(dtype=float)).all():
        raise RuntimeError("ULM activity matrix contains non-finite values")
    return estimate, padj


# Locate frozen inputs.
A2 = latest_release(
    "34A2_DONOR_LEVEL_PROGRAM_MECHANISM_INPUT_AUDIT_",
    ["34A2_DONOR_PROGRAM_SCORE_TABLE.csv", "34A2_MECHANISM_INPUT_AUDIT_LOCK.json"],
)
A2B = latest_release(
    "34A2B_MECHANISM_ENVIRONMENT_RESOURCE_FREEZE_",
    ["34A2B_RESOURCE_FREEZE_LOCK.json", "34A2B_COLLECTRI_HUMAN_NO_COMPLEXES_CORE.csv"],
)
B5 = latest_release(
    "34A1B5_FINAL_CONTRASTIVEVI_CONSENSUS_",
    ["34A1B5_FINAL_PROGRAM_CONSENSUS.csv"],
)
B2 = latest_release(
    "34A1B2_FROZEN_PROGRAM_SALIENT_ALIGNMENT_",
    ["34A1B2_PROGRAM_SCORE_GENES_AND_WEIGHTS.csv"],
)

A2_SCORE_PATH = one_file(A2, "34A2_DONOR_PROGRAM_SCORE_TABLE.csv")
A2_LOCK_PATH = one_file(A2, "34A2_MECHANISM_INPUT_AUDIT_LOCK.json")
A2B_LOCK_PATH = one_file(A2B, "34A2B_RESOURCE_FREEZE_LOCK.json")
CONSENSUS_PATH = one_file(B5, "34A1B5_FINAL_PROGRAM_CONSENSUS.csv")
PROGRAM_GENES_PATH = one_file(B2, "34A1B2_PROGRAM_SCORE_GENES_AND_WEIGHTS.csv")
NO_COMPLEX_CORE_PATH = one_file(A2B, "34A2B_COLLECTRI_HUMAN_NO_COMPLEXES_CORE.csv")

with A2B_LOCK_PATH.open("r", encoding="utf-8") as handle:
    a2b_lock = json.load(handle)
if a2b_lock.get("status") != "GO_FOR_34A3_DONOR_LEVEL_INTRINSIC_TF_ACTIVITY":
    raise RuntimeError(f"34A2B did not authorize 34A3: {a2b_lock.get('status')}")
if str(dc.__version__) != str(a2b_lock.get("decoupler_version")):
    raise RuntimeError(
        f"decoupler version mismatch: expected={a2b_lock.get('decoupler_version')}, actual={dc.__version__}"
    )

consensus = pd.read_csv(CONSENSUS_PATH)
donor_scores = pd.read_csv(A2_SCORE_PATH)
program_genes = pd.read_csv(PROGRAM_GENES_PATH)
no_complex_core = pd.read_csv(NO_COMPLEX_CORE_PATH)

required_consensus = {"lineage", "program_id", "final_tier", "final_role"}
if required_consensus - set(consensus.columns):
    raise RuntimeError("Final consensus table is missing required columns")
consensus["lineage"] = consensus["lineage"].astype(str)
consensus["program_id"] = consensus["program_id"].astype(str)
if len(consensus) != 7 or consensus["program_id"].duplicated().any():
    raise RuntimeError("Expected seven unique frozen programs")

required_scores = {"lineage", "program_id", "score_definition", "donor_id", "condition", "mean_score"}
if required_scores - set(donor_scores.columns):
    raise RuntimeError("Donor score table is missing required columns")
donor_scores = donor_scores.loc[donor_scores["score_definition"] == "TOP100_WEIGHTED"].copy()

required_gene_columns = {"program_id", "score_definition", "gene"}
if required_gene_columns - set(program_genes.columns):
    raise RuntimeError("Program gene table is missing required columns")
top100_genes = program_genes.loc[
    program_genes["score_definition"] == "TOP100_WEIGHTED"
].copy()
top100_genes["gene_norm"] = top100_genes["gene"].map(norm_gene)

OUT.mkdir(parents=True, exist_ok=False)
for folder in [
    "01_pseudobulk",
    "02_networks",
    "03_tf_activities",
    "04_associations",
    "05_review",
    "06_freeze",
    "07_logs",
]:
    (OUT / folder).mkdir(parents=True, exist_ok=True)

pseudobulk_audit_rows = []
network_audit_rows = []
activity_frames = []
association_frames = []

for lineage in ["Oligodendrocytes", "Astrocytes", "Microglia_macrophage"]:
    print("=" * 78, flush=True)
    print(f"Processing lineage: {lineage}", flush=True)
    print("=" * 78, flush=True)

    h5ad_path = find_h5ad(H5AD_NAMES[lineage])
    counts, lognorm, metadata, audit = build_pseudobulk(h5ad_path)
    audit["lineage"] = lineage
    pseudobulk_audit_rows.append(audit)

    np.savez_compressed(
        OUT / "01_pseudobulk" / f"34A3_{lineage}_DONOR_PSEUDOBULK.npz",
        donor_id=metadata["donor_id"].to_numpy(dtype=str),
        condition=metadata["condition"].to_numpy(dtype=str),
        gene=np.asarray(lognorm.columns, dtype=str),
        raw_counts=counts.to_numpy(dtype=np.float64),
        lognorm=lognorm.to_numpy(dtype=np.float64),
    )
    metadata.assign(lineage=lineage).to_csv(
        OUT / "01_pseudobulk" / f"34A3_{lineage}_DONOR_METADATA.csv",
        index=False,
    )

    expression_genes = set(lognorm.columns.astype(str))
    full_network_path = one_file(A2B, f"34A2B_{lineage}_COLLECTRI_TMIN5.csv")
    full_network = prune_network(pd.read_csv(full_network_path), expression_genes)
    full_estimate, full_padj = run_ulm(lognorm, full_network)

    lineage_programs = consensus.loc[consensus["lineage"] == lineage].copy()
    for row in lineage_programs.itertuples(index=False):
        program_id = str(row.program_id)
        print(f"  Program: {program_id}", flush=True)
        excluded_genes = set(
            top100_genes.loc[
                top100_genes["program_id"].astype(str) == program_id,
                "gene_norm",
            ].astype(str)
        )
        if len(excluded_genes) != 100:
            raise RuntimeError(f"{program_id} does not have exactly 100 frozen Top100 genes")

        primary_path = one_file(A2B, f"34A2B_{program_id}_EXCLUDE_TOP100_TMIN5.csv")
        primary_network = prune_network(
            pd.read_csv(primary_path), expression_genes, excluded_genes=None
        )
        sensitivity_network = prune_network(
            no_complex_core, expression_genes, excluded_genes=excluded_genes
        )
        sensitivity_path = (
            OUT / "02_networks" / f"34A3_{program_id}_NO_COMPLEX_EXCLUDE_TOP100_TMIN5.csv"
        )
        sensitivity_network.to_csv(sensitivity_path, index=False)

        network_map = {
            PRIMARY_NETWORK: primary_network,
            SUPPORTIVE_NETWORK: full_network,
            SENSITIVITY_NETWORK: sensitivity_network,
        }
        activity_map = {}
        padj_map = {}
        for network_variant, network in network_map.items():
            if network_variant == SUPPORTIVE_NETWORK:
                estimate, padj = full_estimate, full_padj
            else:
                estimate, padj = run_ulm(lognorm, network)
            activity_map[network_variant] = estimate
            padj_map[network_variant] = padj
            network_audit_rows.append(
                {
                    "lineage": lineage,
                    "program_id": program_id,
                    "network_variant": network_variant,
                    "interactions": int(len(network)),
                    "sources": int(network["source"].nunique()),
                    "targets": int(network["target"].nunique()),
                    "minimum_targets_per_source": int(
                        network.groupby("source", observed=True)["target"].nunique().min()
                    ),
                }
            )

        score_hit = donor_scores.loc[
            (donor_scores["lineage"].astype(str) == lineage)
            & (donor_scores["program_id"].astype(str) == program_id)
        ].copy()
        if len(score_hit) != 8 or score_hit["donor_id"].astype(str).duplicated().any():
            raise RuntimeError(f"{program_id} does not have eight unique donor scores")
        score_hit["donor_id"] = score_hit["donor_id"].astype(str)
        score_hit["condition"] = score_hit["condition"].map(normalize_condition)
        score_lookup = score_hit.set_index("donor_id")
        donor_order = lognorm.index.astype(str).to_numpy()
        condition = metadata.set_index("donor_id").loc[donor_order, "condition"].to_numpy(dtype=str)
        program_score = score_lookup.loc[donor_order, "mean_score"].to_numpy(dtype=float)
        permutation_indices = within_condition_permutation_indices(condition)

        for network_variant in NETWORK_VARIANTS:
            estimate = activity_map[network_variant]
            padj = padj_map[network_variant]
            if list(estimate.index.astype(str)) != list(donor_order):
                estimate = estimate.loc[donor_order]
                padj = padj.loc[donor_order]

            activity_long = (
                estimate.reset_index(names="donor_id")
                .melt(id_vars="donor_id", var_name="tf", value_name="activity")
            )
            padj_long = (
                padj.reset_index(names="donor_id")
                .melt(id_vars="donor_id", var_name="tf", value_name="ulm_padj")
            )
            activity_long = activity_long.merge(
                padj_long,
                on=["donor_id", "tf"],
                how="left",
                validate="one_to_one",
            )
            activity_long["condition"] = activity_long["donor_id"].map(
                dict(zip(donor_order, condition))
            )
            activity_long["lineage"] = lineage
            activity_long["program_id"] = program_id
            activity_long["network_variant"] = network_variant
            activity_frames.append(activity_long)

            association_rows = []
            for tf in estimate.columns.astype(str):
                activity = estimate[tf].to_numpy(dtype=float)
                partial_r, exact_p = exact_partial_permutation_p(
                    activity,
                    program_score,
                    condition,
                    permutation_indices,
                )
                lodo = lodo_partial_summary(activity, program_score, condition, partial_r)
                unadjusted_r = float(np.corrcoef(activity, program_score)[0, 1])
                spearman_result = spearmanr(activity, program_score)
                disease_difference, disease_exact_p = exact_group_difference_p(activity, condition)
                association_rows.append(
                    {
                        "lineage": lineage,
                        "program_id": program_id,
                        "final_tier": str(row.final_tier),
                        "final_role": str(row.final_role),
                        "network_variant": network_variant,
                        "tf": tf,
                        "donors": 8,
                        "condition_adjusted_partial_r": partial_r,
                        "condition_stratified_exact_p": exact_p,
                        "unadjusted_pearson_r": unadjusted_r,
                        "unadjusted_spearman_rho": float(spearman_result.statistic),
                        "activity_tle_minus_normal": disease_difference,
                        "activity_hedges_g": hedges_g(activity, condition),
                        "activity_disease_exact_p": disease_exact_p,
                        **lodo,
                    }
                )
            association = pd.DataFrame(association_rows)
            association["condition_stratified_exact_fdr"] = bh_adjust(
                association["condition_stratified_exact_p"]
            )
            association["activity_disease_exact_fdr"] = bh_adjust(
                association["activity_disease_exact_p"]
            )
            association_frames.append(association)

pseudobulk_audit = pd.DataFrame(pseudobulk_audit_rows)
network_audit = pd.DataFrame(network_audit_rows)
activity_table = pd.concat(activity_frames, ignore_index=True)
association_table = pd.concat(association_frames, ignore_index=True)

primary = association_table.loc[
    association_table["network_variant"] == PRIMARY_NETWORK
].copy()
supportive = association_table.loc[
    association_table["network_variant"] == SUPPORTIVE_NETWORK,
    ["program_id", "tf", "condition_adjusted_partial_r", "lodo_same_sign_fraction"],
].rename(
    columns={
        "condition_adjusted_partial_r": "supportive_partial_r",
        "lodo_same_sign_fraction": "supportive_lodo_same_sign_fraction",
    }
)
sensitivity = association_table.loc[
    association_table["network_variant"] == SENSITIVITY_NETWORK,
    ["program_id", "tf", "condition_adjusted_partial_r", "lodo_same_sign_fraction"],
].rename(
    columns={
        "condition_adjusted_partial_r": "no_complex_partial_r",
        "lodo_same_sign_fraction": "no_complex_lodo_same_sign_fraction",
    }
)

concordance = (
    primary.merge(supportive, on=["program_id", "tf"], how="left", validate="one_to_one")
    .merge(sensitivity, on=["program_id", "tf"], how="left", validate="one_to_one")
)
concordance["supportive_present"] = concordance["supportive_partial_r"].notna()
concordance["no_complex_present"] = concordance["no_complex_partial_r"].notna()
concordance["supportive_sign_concordant"] = (
    np.sign(concordance["condition_adjusted_partial_r"])
    == np.sign(concordance["supportive_partial_r"])
)
concordance["no_complex_sign_concordant"] = np.where(
    concordance["no_complex_present"],
    np.sign(concordance["condition_adjusted_partial_r"])
    == np.sign(concordance["no_complex_partial_r"]),
    np.nan,
)

concordance["provisional_review_tier"] = "NO_SHORTLIST"
strong = (
    concordance["condition_adjusted_partial_r"].abs().ge(0.75)
    & concordance["condition_stratified_exact_p"].le(0.10)
    & concordance["lodo_same_sign_fraction"].ge(0.75)
    & concordance["supportive_sign_concordant"].fillna(False)
    & (
        ~concordance["no_complex_present"]
        | concordance["no_complex_sign_concordant"].fillna(False)
    )
)
moderate = (
    concordance["condition_adjusted_partial_r"].abs().ge(0.60)
    & concordance["condition_stratified_exact_p"].le(0.20)
    & concordance["lodo_same_sign_fraction"].ge(0.625)
    & concordance["supportive_sign_concordant"].fillna(False)
)
concordance.loc[moderate, "provisional_review_tier"] = "MODERATE_EXPLORATORY"
concordance.loc[strong, "provisional_review_tier"] = "STRONG_EXPLORATORY"
concordance["absolute_primary_partial_r"] = concordance[
    "condition_adjusted_partial_r"
].abs()
concordance["rank_within_program"] = (
    concordance.groupby("program_id")["absolute_primary_partial_r"]
    .rank(method="min", ascending=False)
    .astype(int)
)
concordance = concordance.sort_values(
    ["program_id", "rank_within_program", "condition_stratified_exact_p"]
)

top_review = concordance.loc[
    (concordance["rank_within_program"] <= 25)
    | (concordance["provisional_review_tier"] != "NO_SHORTLIST")
].copy()

pseudobulk_audit.to_csv(
    OUT / "01_pseudobulk" / "34A3_PSEUDOBULK_NORMALIZATION_AUDIT.csv", index=False
)
network_audit.to_csv(
    OUT / "02_networks" / "34A3_NETWORK_INPUT_AUDIT.csv", index=False
)
activity_table.to_csv(
    OUT / "03_tf_activities" / "34A3_DONOR_TF_ACTIVITY_LONG.csv.gz",
    index=False,
    compression="gzip",
)
association_table.to_csv(
    OUT / "04_associations" / "34A3_TF_PROGRAM_ASSOCIATIONS_ALL_NETWORKS.csv.gz",
    index=False,
    compression="gzip",
)
concordance.to_csv(
    OUT / "05_review" / "34A3_TF_PROGRAM_CROSS_NETWORK_CONCORDANCE.csv",
    index=False,
)
top_review.to_csv(
    OUT / "05_review" / "34A3_TOP_TF_REVIEW_TABLE.csv",
    index=False,
)

if len(pseudobulk_audit) != 3:
    raise RuntimeError("Expected three lineage pseudobulk audits")
if network_audit.groupby("program_id")["network_variant"].nunique().ne(3).any():
    raise RuntimeError("Each program must have three network variants")
if association_table.empty or concordance.empty:
    raise RuntimeError("TF activity association outputs are empty")

STATUS = "READY_FOR_34A3B_TF_ACTIVITY_REVIEW_AND_TARGET_SELECTION"
NEXT_STAGE = "34A3B_TF_ACTIVITY_REVIEW_AND_TARGET_SELECTION"

lock = {
    "stage": "34A3",
    "created_at": datetime.now().isoformat(),
    "status": STATUS,
    "source_34A2_release": str(A2),
    "source_34A2B_release": str(A2B),
    "source_34A1B5_release": str(B5),
    "programs_evaluated": consensus["program_id"].astype(str).tolist(),
    "network_variants": NETWORK_VARIANTS,
    "normalization": "donor pseudobulk raw-count sums; library-size normalize to 1e4; log1p",
    "tf_method": "decoupler ULM t-value",
    "primary_association": "Pearson correlation after within-condition centering",
    "primary_p_value": "exact 576 permutations restricted within Normal and TLE",
    "lodo_sensitivity": True,
    "primary_network_excludes_program_top100": True,
    "automatic_final_tf_selection_performed": False,
    "provisional_review_tiers_are_exploratory": True,
    "donor_is_inferential_unit": True,
    "cell_level_p_values_calculated": False,
    "causal_claim_authorized": False,
    "next_stage": NEXT_STAGE,
}
(OUT / "06_freeze" / "34A3_TF_ACTIVITY_LOCK.json").write_text(
    json.dumps(lock, ensure_ascii=False, indent=2), encoding="utf-8"
)

review_counts = (
    concordance.groupby(["program_id", "provisional_review_tier"], observed=True)
    .size()
    .reset_index(name="tf_count")
)
review_counts.to_csv(
    OUT / "05_review" / "34A3_PROVISIONAL_TIER_COUNTS.csv", index=False
)

summary_lines = [
    "# 34A3 donor-level intrinsic TF activity",
    "",
    f"## Status: {STATUS}",
    "",
    "- Donor is the inferential unit (n=8 per lineage).",
    "- ULM activity was inferred from library-size-normalized, log1p donor pseudobulk profiles.",
    "- Primary networks exclude each program's frozen Top100 genes.",
    "- Primary association is condition-adjusted and uses all 576 exact within-condition permutations.",
    "- Provisional tiers are for review only and are not final TF claims.",
    "- No causal claim is authorized.",
    "",
    "## Next stage",
    "",
    NEXT_STAGE,
]
(OUT / "34A3_EXECUTIVE_SUMMARY.md").write_text(
    "\n".join(summary_lines), encoding="utf-8"
)

provenance = {
    "stage": "34A3",
    "python": sys.version,
    "decoupler": str(dc.__version__),
    "inputs": {
        str(A2_SCORE_PATH): sha256_file(A2_SCORE_PATH),
        str(A2B_LOCK_PATH): sha256_file(A2B_LOCK_PATH),
        str(CONSENSUS_PATH): sha256_file(CONSENSUS_PATH),
        str(PROGRAM_GENES_PATH): sha256_file(PROGRAM_GENES_PATH),
        str(NO_COMPLEX_CORE_PATH): sha256_file(NO_COMPLEX_CORE_PATH),
    },
}
(OUT / "07_logs" / "34A3_PROVENANCE.json").write_text(
    json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8"
)

output_files = sorted(path for path in OUT.rglob("*") if path.is_file())
manifest = pd.DataFrame(
    [
        {
            "relative_path": path.relative_to(OUT).as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in output_files
    ]
)
manifest.to_csv(OUT / "34A3_FILE_MANIFEST.csv", index=False)

ZIP_PATH = Path(str(OUT) + ".zip")
with zipfile.ZipFile(
    ZIP_PATH,
    mode="w",
    compression=zipfile.ZIP_DEFLATED,
    compresslevel=9,
) as archive:
    for path in sorted(OUT.rglob("*")):
        if path.is_file():
            archive.write(path, arcname=path.relative_to(OUT).as_posix())
ZIP_HASH = sha256_file(ZIP_PATH)
Path(str(ZIP_PATH) + ".sha256.txt").write_text(
    f"{ZIP_HASH}  {ZIP_PATH.name}\n", encoding="ascii"
)

print("=" * 78)
print("34A3 completed")
print("Status:", STATUS)
print("Programs:", consensus["program_id"].nunique(), "/ 7")
print("Activity rows:", len(activity_table))
print("Association rows:", len(association_table))
print("Cross-network review rows:", len(concordance))
print("Next stage:", NEXT_STAGE)
print("ZIP:", ZIP_PATH)
print("SHA256:", ZIP_HASH)
print("=" * 78)
