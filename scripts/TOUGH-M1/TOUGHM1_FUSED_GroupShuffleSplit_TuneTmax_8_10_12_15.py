#!/usr/bin/env python3
"""
TOUGH-M1 | Group-based train/test evaluation with training-only t_max selection

- Uses raw FUSED curve coordinates exported over [4.8, 15.0] A.
- Uses RCSB 30% sequence-identity clusters as split groups.
- Uses 20 repeated random 80/20 train/test splits at the cluster level.
- Inside each outer training split, creates a cluster-level validation split.
- Selects t_max from {8, 10, 12, 15} A using validation AUC only.
- Evaluates the selected interval once on held-out test pairs.
- Runs both:
    1) direct split-fitted MFPC-FUSED similarity
    2) supervised Siamese FUSED model with |z_i-z_j|, z_i*z_j, cosine
- Resumable: completed split rows are skipped unless OVERWRITE=1 is set.
"""

from __future__ import annotations

import importlib.util
import math
import os
import random
import re
import sys
import time

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupShuffleSplit
from torch import nn
from torch.utils.data import DataLoader


ROOT = os.getcwd()
BASE_SCRIPT = os.path.join(
    ROOT,
    "scripts",
    "TOUGH-M1",
    "TOUGHM1_Siamese_Base.py",
)
spec = importlib.util.spec_from_file_location("toughm1_siamese_base", BASE_SCRIPT)
base = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = base
spec.loader.exec_module(base)


CACHE_DIR = os.path.join(ROOT, "data", "cache")
OUT_DIR = os.path.join(ROOT, "results", "TOUGH-M1")
INSPECT_DIR = os.path.join(ROOT, "data", "TOUGH-M1")
os.makedirs(OUT_DIR, exist_ok=True)

FEATURE_FILE = os.path.join(CACHE_DIR, "TOUGHM1_rawFUSED_features_thr4p8_15p0.csv")
POSITIVE_FILE = os.path.join(INSPECT_DIR, "TOUGH-M1_positive.list")
NEGATIVE_FILE = os.path.join(INSPECT_DIR, "TOUGH-M1_negative.list")
GROUP_FILE = os.path.join(
    INSPECT_DIR,
    "sequence_clusters",
    "TOUGHM1_code5_rcsb30_cluster_mapping.csv",
)

SPLIT_RESULTS_FILE = os.path.join(
    OUT_DIR,
    "TOUGHM1_FUSED_RCSB30_GroupShuffleSplit_TuneTmax_8_10_12_15_results.csv",
)
SUMMARY_FILE = os.path.join(
    OUT_DIR,
    "TOUGHM1_FUSED_RCSB30_GroupShuffleSplit_TuneTmax_8_10_12_15_summary.csv",
)

T_MAX_GRID = [8.0, 10.0, 12.0, 15.0]
T_MIN = 4.8
N_SPLITS = 20
OUTER_TEST_SIZE = 0.20
INNER_VALIDATION_SIZE = 0.20
SEED = 41001
VAR_TARGET = 0.95
MAX_PC = 20


class Config(base.Config):
    n_splits: int = N_SPLITS
    outer_test_size: float = OUTER_TEST_SIZE
    inner_validation_size: float = INNER_VALIDATION_SIZE
    seed: int = SEED


CFG = Config()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def completed_splits(path: str) -> set[int]:
    if not os.path.exists(path):
        return set()
    df = pd.read_csv(path)
    if "split" not in df.columns:
        return set()
    return set(df["split"].astype(int))


def append_csv_row(path: str, row: dict) -> None:
    pd.DataFrame([row]).to_csv(path, mode="a", index=False, header=not os.path.exists(path))


def load_pairs() -> pd.DataFrame:
    positive = pd.read_csv(
        POSITIVE_FILE,
        sep=r"\s+",
        header=None,
        usecols=[0, 1],
        names=["id1", "id2"],
        engine="python",
    )
    positive["label"] = 1
    negative = pd.read_csv(
        NEGATIVE_FILE,
        sep=r"\s+",
        header=None,
        usecols=[0, 1],
        names=["id1", "id2"],
        engine="python",
    )
    negative["label"] = 0
    return pd.concat([positive, negative], ignore_index=True)


def select_within_pairs(pairs: pd.DataFrame, ids: set[str]) -> pd.DataFrame:
    mask = pairs["id1"].isin(ids).to_numpy() & pairs["id2"].isin(ids).to_numpy()
    return pairs.loc[mask].reset_index(drop=True)


def has_two_classes(df: pd.DataFrame) -> bool:
    return df["label"].nunique() == 2


def make_inner_split(train_outer_ids: np.ndarray, id_to_group: dict[str, str], seed: int):
    groups = np.array([id_to_group[x] for x in train_outer_ids])
    splitter = GroupShuffleSplit(
        n_splits=1,
        test_size=INNER_VALIDATION_SIZE,
        random_state=seed,
    )
    train_idx, val_idx = next(splitter.split(train_outer_ids, groups=groups))
    return set(train_outer_ids[train_idx]), set(train_outer_ids[val_idx])


def parse_feature_threshold(col: str) -> float | None:
    match = re.search(r"_t([0-9]+)p([0-9]+)$", col)
    if match is None:
        return None
    return float(f"{match.group(1)}.{match.group(2)}")


def feature_columns_for_tmax(feature_cols: list[str], t_max: float) -> list[str]:
    keep: list[str] = []
    for col in feature_cols:
        tt = parse_feature_threshold(col)
        if tt is not None and T_MIN - 1e-8 <= tt <= t_max + 1e-8:
            keep.append(col)
    if not keep:
        raise RuntimeError(f"No feature columns found for t_max={t_max}.")
    return keep


def pair_arrays(pairs: pd.DataFrame, id_to_row: dict[str, int]):
    idx1 = pairs["id1"].map(id_to_row).to_numpy(dtype=np.int64)
    idx2 = pairs["id2"].map(id_to_row).to_numpy(dtype=np.int64)
    labels = pairs["label"].to_numpy(dtype=np.float32)
    return idx1, idx2, labels


def fit_project_train_only(
    X: np.ndarray,
    train_rows: np.ndarray,
    feature_cols: list[str],
) -> tuple[np.ndarray, int, float, float]:
    colnames = np.array(feature_cols)
    ilr_cols = np.array([name.startswith("ILR") for name in colnames])
    cdpa_cols = ~ilr_cols

    X_scaled = X.astype(np.float64, copy=True)
    ilr_energy = np.mean(np.abs(X_scaled[np.ix_(train_rows, ilr_cols)]))
    cdpa_energy = np.mean(np.abs(X_scaled[np.ix_(train_rows, cdpa_cols)]))
    scale_ilr = 1.0 if ilr_energy == 0 else cdpa_energy / ilr_energy
    X_scaled[:, ilr_cols] *= scale_ilr

    center = X_scaled[train_rows].mean(axis=0)
    X_centered = X_scaled - center
    X_train = X_centered[train_rows]

    cov = np.cov(X_train, rowvar=False)
    values, vectors = np.linalg.eigh(cov)
    order = np.argsort(values)[::-1]
    values = np.maximum(values[order], 0.0)
    vectors = vectors[:, order]

    keep = values > np.finfo(float).eps
    values = values[keep]
    vectors = vectors[:, keep]
    prop = values / values.sum()
    cum = np.cumsum(prop)
    k_use = int(np.searchsorted(cum, VAR_TARGET) + 1)
    k_use = min(k_use, MAX_PC, vectors.shape[1])

    scores = X_centered @ vectors[:, :k_use]
    explained = float(prop[:k_use].sum())
    return scores, k_use, explained, float(scale_ilr)


def auc_from_scores(scores: np.ndarray, pairs: pd.DataFrame, id_to_row: dict[str, int]) -> float:
    idx1, idx2, labels = pair_arrays(pairs, id_to_row)
    diff = scores[idx1] - scores[idx2]
    distance = np.sqrt(np.sum(diff * diff, axis=1))
    return float(roc_auc_score(labels, -distance))


def evaluate_direct_candidates(
    feature_df: pd.DataFrame,
    feature_cols_all: list[str],
    train_ids: set[str],
    val_pairs: pd.DataFrame,
    test_pairs: pd.DataFrame,
    id_to_row: dict[str, int],
) -> dict:
    train_rows = np.array([id_to_row[x] for x in train_ids], dtype=np.int64)
    best: dict | None = None

    for t_max in T_MAX_GRID:
        cols = feature_columns_for_tmax(feature_cols_all, t_max)
        X = feature_df[cols].to_numpy(dtype=np.float64)
        scores, k_use, explained, scale_ilr = fit_project_train_only(X, train_rows, cols)
        val_auc = auc_from_scores(scores, val_pairs, id_to_row)
        candidate = {
            "t_max": t_max,
            "validation_auc": val_auc,
            "scores": scores,
            "k_use": k_use,
            "explained": explained,
            "scale_ilr": scale_ilr,
            "input_features": len(cols),
        }
        if best is None or val_auc > best["validation_auc"]:
            best = candidate

    assert best is not None
    test_auc = auc_from_scores(best["scores"], test_pairs, id_to_row)
    return {
        "direct_selected_t_max": best["t_max"],
        "direct_validation_auc": best["validation_auc"],
        "direct_test_auc": test_auc,
        "direct_retained_scores": best["k_use"],
        "direct_explained_variance": best["explained"],
        "direct_scale_ilr": best["scale_ilr"],
        "direct_input_features": best["input_features"],
    }


def train_siamese_candidate(
    split_id: int,
    t_max: float,
    X_raw: np.ndarray,
    train_ids: set[str],
    train_pairs: pd.DataFrame,
    val_pairs: pd.DataFrame,
    id_to_row: dict[str, int],
    device: torch.device,
) -> tuple[base.SiameseFusedNet, float, int]:
    train_rows = np.array([id_to_row[x] for x in train_ids], dtype=np.int64)
    features = base.standardize_by_train_pockets(X_raw, train_rows)

    tr_i1, tr_i2, tr_y = pair_arrays(train_pairs, id_to_row)
    va_i1, va_i2, va_y = pair_arrays(val_pairs, id_to_row)

    train_loader = DataLoader(
        base.PairDataset(features, tr_i1, tr_i2, tr_y),
        batch_size=CFG.batch_size,
        shuffle=True,
        num_workers=CFG.num_workers,
    )
    val_loader = DataLoader(base.PairDataset(features, va_i1, va_i2, va_y), batch_size=CFG.batch_size)

    model_seed = CFG.seed + 100000 * split_id + int(t_max * 10)
    set_seed(model_seed)
    model = base.SiameseFusedNet(features.shape[1], CFG).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=CFG.learning_rate, weight_decay=CFG.weight_decay)
    loss_fn = nn.BCEWithLogitsLoss()

    best_state = None
    best_val_auc = -math.inf
    best_epoch = 0
    stale = 0

    for epoch in range(1, CFG.max_epochs + 1):
        model.train()
        total_loss = 0.0
        n_seen = 0
        for x1, x2, y in train_loader:
            x1 = x1.to(device)
            x2 = x2.to(device)
            y = y.to(device)
            opt.zero_grad(set_to_none=True)
            loss = loss_fn(model(x1, x2), y)
            loss.backward()
            opt.step()
            total_loss += float(loss.item()) * y.shape[0]
            n_seen += y.shape[0]

        val_auc = base.evaluate_auc(model, val_loader, device)
        print(
            f"Split {split_id:02d} | t_max={t_max:04.1f} | epoch {epoch:02d} | "
            f"loss={total_loss / max(1, n_seen):.5f} | val_auc={val_auc:.4f}",
            flush=True,
        )

        if val_auc > best_val_auc + 1e-4:
            best_val_auc = val_auc
            best_epoch = epoch
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
            if stale >= CFG.patience:
                break

    if best_state is not None:
        model.load_state_dict(best_state)
    return model, best_val_auc, best_epoch


def evaluate_siamese_candidates(
    split_id: int,
    feature_df: pd.DataFrame,
    feature_cols_all: list[str],
    train_ids: set[str],
    train_pairs: pd.DataFrame,
    val_pairs: pd.DataFrame,
    test_pairs: pd.DataFrame,
    id_to_row: dict[str, int],
    device: torch.device,
) -> dict:
    best: dict | None = None

    for t_max in T_MAX_GRID:
        cols = feature_columns_for_tmax(feature_cols_all, t_max)
        X_raw = feature_df[cols].to_numpy(dtype=np.float32)
        model, val_auc, best_epoch = train_siamese_candidate(
            split_id, t_max, X_raw, train_ids, train_pairs, val_pairs, id_to_row, device
        )
        candidate = {
            "t_max": t_max,
            "validation_auc": val_auc,
            "best_epoch": best_epoch,
            "model": model,
            "X_raw": X_raw,
            "input_features": len(cols),
        }
        if best is None or val_auc > best["validation_auc"]:
            best = candidate

    assert best is not None
    train_rows = np.array([id_to_row[x] for x in train_ids], dtype=np.int64)
    features = base.standardize_by_train_pockets(best["X_raw"], train_rows)
    te_i1, te_i2, te_y = pair_arrays(test_pairs, id_to_row)
    test_loader = DataLoader(base.PairDataset(features, te_i1, te_i2, te_y), batch_size=CFG.batch_size)
    test_pred = base.predict_scores(best["model"], test_loader, device)
    test_auc = float(roc_auc_score(te_y, test_pred))

    return {
        "siamese_selected_t_max": best["t_max"],
        "siamese_validation_auc": best["validation_auc"],
        "siamese_test_auc": test_auc,
        "siamese_best_epoch": best["best_epoch"],
        "siamese_input_features": best["input_features"],
    }


def run_one_split(
    split_id: int,
    train_outer_ids: np.ndarray,
    test_ids: np.ndarray,
    feature_df: pd.DataFrame,
    feature_cols_all: list[str],
    pairs: pd.DataFrame,
    id_to_row: dict[str, int],
    id_to_group: dict[str, str],
    device: torch.device,
) -> dict:
    split_start = time.time()
    print("\n" + "=" * 72, flush=True)
    print(f"Group split {split_id} of {N_SPLITS}", flush=True)
    print("=" * 72, flush=True)

    train_ids, val_ids = make_inner_split(train_outer_ids, id_to_group, SEED + 1000 * split_id)
    test_ids_set = set(test_ids)

    train_pairs = select_within_pairs(pairs, train_ids)
    val_pairs = select_within_pairs(pairs, val_ids)
    test_pairs = select_within_pairs(pairs, test_ids_set)

    if not (has_two_classes(train_pairs) and has_two_classes(val_pairs) and has_two_classes(test_pairs)):
        raise RuntimeError(f"Split {split_id} has a partition without both classes.")

    print(f"Fitting pockets: {len(train_ids)}", flush=True)
    print(f"Validation pockets: {len(val_ids)}", flush=True)
    print(f"Test pockets: {len(test_ids_set)}", flush=True)
    print(f"Training pairs: {len(train_pairs)}", flush=True)
    print(f"Validation pairs: {len(val_pairs)}", flush=True)
    print(f"Test pairs: {len(test_pairs)}", flush=True)

    direct = evaluate_direct_candidates(
        feature_df, feature_cols_all, train_ids, val_pairs, test_pairs, id_to_row
    )
    print(
        f"Split {split_id:02d} | direct selected t_max={direct['direct_selected_t_max']:.1f} | "
        f"test AUC={direct['direct_test_auc']:.4f}",
        flush=True,
    )

    siamese = evaluate_siamese_candidates(
        split_id,
        feature_df,
        feature_cols_all,
        train_ids,
        train_pairs,
        val_pairs,
        test_pairs,
        id_to_row,
        device,
    )
    print(
        f"Split {split_id:02d} | Siamese selected t_max={siamese['siamese_selected_t_max']:.1f} | "
        f"test AUC={siamese['siamese_test_auc']:.4f}",
        flush=True,
    )

    elapsed_minutes = (time.time() - split_start) / 60.0
    return {
        "split": split_id,
        "t_min": T_MIN,
        "t_max_grid": ";".join(f"{x:g}" for x in T_MAX_GRID),
        "n_fit_pockets": len(train_ids),
        "n_validation_pockets": len(val_ids),
        "n_test_pockets": len(test_ids_set),
        "n_train_pairs": len(train_pairs),
        "n_validation_pairs": len(val_pairs),
        "n_test_pairs": len(test_pairs),
        "n_positive_test_pairs": int((test_pairs["label"] == 1).sum()),
        "n_negative_test_pairs": int((test_pairs["label"] == 0).sum()),
        "auc_difference_siamese_minus_direct": siamese["siamese_test_auc"] - direct["direct_test_auc"],
        "split_method": "GroupShuffleSplit",
        "split_unit": "RCSB30_sequence_cluster",
        "outer_test_size": OUTER_TEST_SIZE,
        "inner_validation_size": INNER_VALIDATION_SIZE,
        "embedding_dim": CFG.embedding_dim,
        "dropout": CFG.dropout,
        "learning_rate": CFG.learning_rate,
        "weight_decay": CFG.weight_decay,
        "batch_size": CFG.batch_size,
        "max_epochs": CFG.max_epochs,
        "patience": CFG.patience,
        "elapsed_minutes": elapsed_minutes,
        **direct,
        **siamese,
    }


def write_summary() -> None:
    df = pd.read_csv(SPLIT_RESULTS_FILE)
    summary = pd.DataFrame(
        {
            "method": [
                "direct_MFPC_FUSED_similarity_validation_tuned_tmax",
                "supervised_Siamese_FUSED_validation_tuned_tmax",
            ],
            "split_method": ["GroupShuffleSplit", "GroupShuffleSplit"],
            "split_unit": ["RCSB30_sequence_cluster", "RCSB30_sequence_cluster"],
            "n_splits": [len(df), len(df)],
            "outer_test_size": [OUTER_TEST_SIZE, OUTER_TEST_SIZE],
            "inner_validation_size": [INNER_VALIDATION_SIZE, INNER_VALIDATION_SIZE],
            "t_min": [T_MIN, T_MIN],
            "t_max_grid": [";".join(f"{x:g}" for x in T_MAX_GRID)] * 2,
            "mean_auc": [df.direct_test_auc.mean(), df.siamese_test_auc.mean()],
            "sd_auc": [df.direct_test_auc.std(ddof=1), df.siamese_test_auc.std(ddof=1)],
            "mean_selected_t_max": [
                df.direct_selected_t_max.mean(),
                df.siamese_selected_t_max.mean(),
            ],
            "sd_selected_t_max": [
                df.direct_selected_t_max.std(ddof=1),
                df.siamese_selected_t_max.std(ddof=1),
            ],
            "mean_validation_auc": [
                df.direct_validation_auc.mean(),
                df.siamese_validation_auc.mean(),
            ],
            "sd_validation_auc": [
                df.direct_validation_auc.std(ddof=1),
                df.siamese_validation_auc.std(ddof=1),
            ],
            "mean_test_pairs": [df.n_test_pairs.mean(), df.n_test_pairs.mean()],
            "mean_positive_test_pairs": [
                df.n_positive_test_pairs.mean(),
                df.n_positive_test_pairs.mean(),
            ],
            "mean_negative_test_pairs": [
                df.n_negative_test_pairs.mean(),
                df.n_negative_test_pairs.mean(),
            ],
            "mean_auc_difference_siamese_minus_direct": [
                np.nan,
                df.auc_difference_siamese_minus_direct.mean(),
            ],
            "sd_auc_difference_siamese_minus_direct": [
                np.nan,
                df.auc_difference_siamese_minus_direct.std(ddof=1),
            ],
            "embedding_dim": [np.nan, CFG.embedding_dim],
            "dropout": [np.nan, CFG.dropout],
            "learning_rate": [np.nan, CFG.learning_rate],
            "weight_decay": [np.nan, CFG.weight_decay],
            "total_elapsed_minutes": [df.elapsed_minutes.sum(), df.elapsed_minutes.sum()],
            "mean_elapsed_minutes_per_split": [df.elapsed_minutes.mean(), df.elapsed_minutes.mean()],
        }
    )
    summary.to_csv(SUMMARY_FILE, index=False)
    print("\n===== TOUGH-M1 t_max-tuned FUSED GroupShuffleSplit summary =====", flush=True)
    print(summary.to_string(index=False), flush=True)


def main() -> None:
    run_start = time.time()
    if os.environ.get("OVERWRITE", "0") == "1":
        for path in [SPLIT_RESULTS_FILE, SUMMARY_FILE]:
            if os.path.exists(path):
                os.remove(path)

    for path in [FEATURE_FILE, POSITIVE_FILE, NEGATIVE_FILE, GROUP_FILE]:
        if not os.path.exists(path):
            raise SystemExit(f"Missing required file: {path}")

    feature_df = pd.read_csv(FEATURE_FILE)
    group_df = pd.read_csv(GROUP_FILE)
    group_df = group_df[["code5", "rcsb30_cluster_id", "mapping_status"]].drop_duplicates()
    group_df = group_df[group_df["mapping_status"] == "mapped"].copy()
    feature_df = feature_df.merge(group_df, on="code5", how="inner")

    feature_cols_all = [
        col for col in feature_df.columns
        if col not in {"code5", "rcsb30_cluster_id", "mapping_status"}
    ]
    ids = feature_df["code5"].to_numpy()
    id_to_row = {id_: i for i, id_ in enumerate(ids)}
    id_to_group = dict(zip(feature_df["code5"], feature_df["rcsb30_cluster_id"]))
    groups = np.array([id_to_group[x] for x in ids])

    pairs = load_pairs()
    pairs = pairs[pairs["id1"].isin(id_to_row) & pairs["id2"].isin(id_to_row)].reset_index(drop=True)
    pairs["label"] = pairs["label"].astype(np.int64)

    print(f"Mapped usable FUSED pockets: {len(ids)}", flush=True)
    print(f"RCSB30 groups: {len(set(groups))}", flush=True)
    print(f"Analyzed pairs after dropping unmapped/unusable pockets: {len(pairs)}", flush=True)
    print(f"Positive pairs: {int((pairs.label == 1).sum())}", flush=True)
    print(f"Negative pairs: {int((pairs.label == 0).sum())}", flush=True)
    print(f"t_min: {T_MIN}", flush=True)
    print(f"t_max grid: {T_MAX_GRID}", flush=True)
    print(f"Full exported feature dimension: {len(feature_cols_all)}", flush=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}", flush=True)

    splitter = GroupShuffleSplit(
        n_splits=N_SPLITS,
        test_size=OUTER_TEST_SIZE,
        random_state=SEED,
    )
    done = completed_splits(SPLIT_RESULTS_FILE)
    for split_id, (train_idx, test_idx) in enumerate(splitter.split(ids, groups=groups), start=1):
        if split_id in done:
            print(f"Skipping completed split {split_id}.", flush=True)
            continue
        row = run_one_split(
            split_id,
            ids[train_idx],
            ids[test_idx],
            feature_df,
            feature_cols_all,
            pairs,
            id_to_row,
            id_to_group,
            device,
        )
        append_csv_row(SPLIT_RESULTS_FILE, row)
        write_summary()

    write_summary()
    print(f"Current script-session elapsed time: {(time.time() - run_start) / 60.0:.2f} minutes", flush=True)
    print("\nSaved outputs:", flush=True)
    print(f"  {SPLIT_RESULTS_FILE}", flush=True)
    print(f"  {SUMMARY_FILE}", flush=True)


if __name__ == "__main__":
    main()
