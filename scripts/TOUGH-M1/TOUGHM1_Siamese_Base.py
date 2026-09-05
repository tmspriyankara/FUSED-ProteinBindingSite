#!/usr/bin/env python3
"""
Shared TOUGH-M1 Siamese FUSED utilities.

This file defines the neural-network architecture, pair dataset, scaling,
and AUC helper functions used by the TOUGH-M1 analysis scripts. The final
TOUGH-M1 analysis is run from:

    scripts/TOUGH-M1/TOUGHM1_FUSED_GroupShuffleSplit_TuneTmax_8_10_12_15.py
"""

from __future__ import annotations

import math
import os
import random
from dataclasses import dataclass

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score
from torch import nn
from torch.utils.data import DataLoader, Dataset


ROOT = os.getcwd()
CACHE_DIR = os.path.join(ROOT, "data", "cache")
OUT_DIR = os.path.join(ROOT, "results", "TOUGH-M1")
os.makedirs(OUT_DIR, exist_ok=True)

FEATURE_FILE = os.path.join(CACHE_DIR, "TOUGHM1_rawFUSED_features_thr4p8_10p0.csv")
PAIR_FILE = os.path.join(OUT_DIR, "TOUGHM1_FUSED_pair_scores_thr4p8_10p0.csv")
GROUP_FILE = os.path.join(
    ROOT,
    "data",
    "TOUGH-M1",
    "sequence_clusters",
    "TOUGHM1_code5_rcsb30_cluster_mapping.csv",
)

FOLD_RESULTS_FILE = os.path.join(
    OUT_DIR, "TOUGHM1_FUSED_deep_siamese_rawcurves_RCSB30CV_absProdCos_repeated5x5_fold_results.csv"
)
SUMMARY_FILE = os.path.join(
    OUT_DIR, "TOUGHM1_FUSED_deep_siamese_rawcurves_RCSB30CV_absProdCos_repeated5x5_summary.csv"
)


@dataclass
class Config:
    n_repeats: int = 5
    k_folds: int = 5
    seed: int = 31001
    batch_size: int = 4096
    max_epochs: int = 20
    patience: int = 4
    learning_rate: float = 1e-3
    weight_decay: float = 1e-4
    hidden1: int = 256
    hidden2: int = 128
    embedding_dim: int = 64
    dropout: float = 0.20
    num_workers: int = 0


CFG = Config()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def make_group_folds(
    ids: np.ndarray,
    id_to_group: dict[str, str],
    k: int,
    seed: int,
) -> list[np.ndarray]:
    group_df = pd.DataFrame({"id": ids, "group": [id_to_group[x] for x in ids]})
    sizes = group_df.groupby("group", as_index=False).size().rename(columns={"size": "n"})
    rng = np.random.default_rng(seed)
    sizes["jitter"] = rng.random(len(sizes))
    sizes = sizes.sort_values(["n", "jitter"], ascending=[False, True])

    fold_groups: list[list[str]] = [[] for _ in range(k)]
    fold_counts = np.zeros(k, dtype=np.int64)
    for row in sizes.itertuples(index=False):
        fold_id = int(np.argmin(fold_counts))
        fold_groups[fold_id].append(row.group)
        fold_counts[fold_id] += int(row.n)

    folds: list[np.ndarray] = []
    for fold_id in range(k):
        groups_here = set(fold_groups[fold_id])
        folds.append(group_df.loc[group_df["group"].isin(groups_here), "id"].to_numpy())
    return folds


class PairDataset(Dataset):
    def __init__(
        self,
        features: np.ndarray,
        idx1: np.ndarray,
        idx2: np.ndarray,
        labels: np.ndarray,
    ) -> None:
        self.features = features
        self.idx1 = idx1.astype(np.int64)
        self.idx2 = idx2.astype(np.int64)
        self.labels = labels.astype(np.float32)

    def __len__(self) -> int:
        return len(self.labels)

    def __getitem__(self, i: int):
        return (
            torch.from_numpy(self.features[self.idx1[i]]),
            torch.from_numpy(self.features[self.idx2[i]]),
            torch.tensor(self.labels[i], dtype=torch.float32),
        )


class SiameseFusedNet(nn.Module):
    def __init__(self, input_dim: int, cfg: Config):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, cfg.hidden1),
            nn.ReLU(),
            nn.BatchNorm1d(cfg.hidden1),
            nn.Dropout(cfg.dropout),
            nn.Linear(cfg.hidden1, cfg.hidden2),
            nn.ReLU(),
            nn.BatchNorm1d(cfg.hidden2),
            nn.Dropout(cfg.dropout),
            nn.Linear(cfg.hidden2, cfg.embedding_dim),
            nn.ReLU(),
        )
        pair_dim = cfg.embedding_dim * 2 + 1
        self.head = nn.Sequential(
            nn.Linear(pair_dim, 128),
            nn.ReLU(),
            nn.Dropout(cfg.dropout),
            nn.Linear(128, 32),
            nn.ReLU(),
            nn.Linear(32, 1),
        )

    def forward(self, x1: torch.Tensor, x2: torch.Tensor) -> torch.Tensor:
        z1 = self.encoder(x1)
        z2 = self.encoder(x2)
        absdiff = torch.abs(z1 - z2)
        prod = z1 * z2
        cosine = nn.functional.cosine_similarity(z1, z2).unsqueeze(1)
        pair = torch.cat([absdiff, prod, cosine], dim=1)
        return self.head(pair).squeeze(1)


def pair_arrays(pairs: pd.DataFrame, id_to_row: dict[str, int]):
    idx1 = pairs["id1"].map(id_to_row).to_numpy()
    idx2 = pairs["id2"].map(id_to_row).to_numpy()
    labels = pairs["label"].to_numpy(dtype=np.float32)
    return idx1, idx2, labels


def standardize_by_train_pockets(X: np.ndarray, train_pocket_rows: np.ndarray) -> np.ndarray:
    mu = X[train_pocket_rows].mean(axis=0)
    sd = X[train_pocket_rows].std(axis=0)
    sd[sd < 1e-8] = 1.0
    return ((X - mu) / sd).astype(np.float32)


def evaluate_auc(model: nn.Module, loader: DataLoader, device: torch.device) -> float:
    model.eval()
    preds: list[np.ndarray] = []
    labels: list[np.ndarray] = []
    with torch.no_grad():
        for x1, x2, y in loader:
            logits = model(x1.to(device), x2.to(device))
            preds.append(torch.sigmoid(logits).cpu().numpy())
            labels.append(y.numpy())
    return float(roc_auc_score(np.concatenate(labels), np.concatenate(preds)))


def predict_scores(model: nn.Module, loader: DataLoader, device: torch.device) -> np.ndarray:
    model.eval()
    preds: list[np.ndarray] = []
    with torch.no_grad():
        for x1, x2, _ in loader:
            logits = model(x1.to(device), x2.to(device))
            preds.append(torch.sigmoid(logits).cpu().numpy())
    return np.concatenate(preds)


def append_csv_row(path: str, row: dict) -> None:
    pd.DataFrame([row]).to_csv(path, mode="a", index=False, header=not os.path.exists(path))


def completed_repeat_folds(path: str) -> set[tuple[int, int]]:
    if not os.path.exists(path):
        return set()
    df = pd.read_csv(path)
    if not {"repeat", "fold"}.issubset(df.columns):
        return set()
    return set(zip(df["repeat"].astype(int), df["fold"].astype(int)))


def train_one_fold(
    repeat_id: int,
    fold_id: int,
    features_raw: np.ndarray,
    pairs: pd.DataFrame,
    folds: list[np.ndarray],
    id_to_row: dict[str, int],
    id_to_group: dict[str, str],
    device: torch.device,
) -> dict:
    print("\n" + "=" * 70, flush=True)
    print(f"Repeat {repeat_id + 1} of {CFG.n_repeats} | sequence-cluster fold {fold_id + 1} of {CFG.k_folds}", flush=True)
    print("=" * 70, flush=True)

    test_ids = set(folds[fold_id])
    train_outer_ids = np.array([x for i, fold in enumerate(folds) if i != fold_id for x in fold])
    inner_seed = CFG.seed + 100000 * (repeat_id + 1) + 100 + fold_id
    inner_folds = make_group_folds(train_outer_ids, id_to_group, CFG.k_folds - 1, inner_seed)
    val_fold_local = (fold_id + 1) % (CFG.k_folds - 1)
    val_ids = set(inner_folds[val_fold_local])
    train_ids = set(x for i, fold in enumerate(inner_folds) if i != val_fold_local for x in fold)

    id1 = pairs["id1"]
    id2 = pairs["id2"]
    is_train = id1.isin(train_ids).to_numpy() & id2.isin(train_ids).to_numpy()
    is_val = id1.isin(val_ids).to_numpy() & id2.isin(val_ids).to_numpy()
    is_test = id1.isin(test_ids).to_numpy() & id2.isin(test_ids).to_numpy()

    train_pairs = pairs.loc[is_train].reset_index(drop=True)
    val_pairs = pairs.loc[is_val].reset_index(drop=True)
    test_pairs = pairs.loc[is_test].reset_index(drop=True)

    print(f"Training pairs: {len(train_pairs)}", flush=True)
    print(f"Validation pairs: {len(val_pairs)}", flush=True)
    print(f"Test pairs: {len(test_pairs)}", flush=True)
    print(f"Test positives: {int((test_pairs.label == 1).sum())}", flush=True)
    print(f"Test negatives: {int((test_pairs.label == 0).sum())}", flush=True)

    train_pocket_rows = np.array([id_to_row[x] for x in train_ids], dtype=np.int64)
    features = standardize_by_train_pockets(features_raw, train_pocket_rows)

    tr_i1, tr_i2, tr_y = pair_arrays(train_pairs, id_to_row)
    va_i1, va_i2, va_y = pair_arrays(val_pairs, id_to_row)
    te_i1, te_i2, te_y = pair_arrays(test_pairs, id_to_row)

    train_loader = DataLoader(
        PairDataset(features, tr_i1, tr_i2, tr_y),
        batch_size=CFG.batch_size,
        shuffle=True,
        num_workers=CFG.num_workers,
        pin_memory=False,
    )
    val_loader = DataLoader(PairDataset(features, va_i1, va_i2, va_y), batch_size=CFG.batch_size)
    test_loader = DataLoader(PairDataset(features, te_i1, te_i2, te_y), batch_size=CFG.batch_size)

    model_seed = CFG.seed + 1000000 * (repeat_id + 1) + 1000 * (fold_id + 1)
    set_seed(model_seed)
    model = SiameseFusedNet(features.shape[1], CFG).to(device)
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

        val_auc = evaluate_auc(model, val_loader, device)
        print(
            f"Repeat {repeat_id + 1} | fold {fold_id + 1} | epoch {epoch:02d} | "
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
                print(f"Early stopping at epoch {epoch}.", flush=True)
                break

    if best_state is not None:
        model.load_state_dict(best_state)

    test_pred = predict_scores(model, test_loader, device)
    deep_auc = float(roc_auc_score(te_y, test_pred))
    zero_auc = float(roc_auc_score(test_pairs["label"], test_pairs["fused_similarity"]))

    print(
        f"Repeat {repeat_id + 1} | fold {fold_id + 1} | zero-shot AUC={zero_auc:.4f} | "
        f"deep AUC={deep_auc:.4f} | diff={deep_auc - zero_auc:.4f}",
        flush=True,
    )

    return {
        "repeat": repeat_id + 1,
        "fold": fold_id + 1,
        "n_train_pairs": len(train_pairs),
        "n_validation_pairs": len(val_pairs),
        "n_test_pairs": len(test_pairs),
        "n_positive_test_pairs": int((test_pairs.label == 1).sum()),
        "n_negative_test_pairs": int((test_pairs.label == 0).sum()),
        "zero_shot_auc": zero_auc,
        "deep_siamese_auc": deep_auc,
        "auc_difference": deep_auc - zero_auc,
        "best_validation_auc": best_val_auc,
        "best_epoch": best_epoch,
        "input_features": features.shape[1],
        "embedding_dim": CFG.embedding_dim,
        "dropout": CFG.dropout,
        "learning_rate": CFG.learning_rate,
        "weight_decay": CFG.weight_decay,
        "pairwise_features": "absolute_difference_product_cosine",
        "max_epochs": CFG.max_epochs,
        "batch_size": CFG.batch_size,
        "split_unit": "RCSB30_sequence_cluster",
    }


def write_summary() -> None:
    fold_df = pd.read_csv(FOLD_RESULTS_FILE)
    summary_df = pd.DataFrame(
        {
            "setting": [
                "zero_shot_FUSED_distance",
                "fine_tuned_FUSED_siamese_abs_product_cosine",
            ],
            "split_unit": ["RCSB30_sequence_cluster", "RCSB30_sequence_cluster"],
            "n_repeats": [fold_df["repeat"].nunique(), fold_df["repeat"].nunique()],
            "n_folds": [CFG.k_folds, CFG.k_folds],
            "n_outer_evals": [len(fold_df), len(fold_df)],
            "mean_auc": [fold_df.zero_shot_auc.mean(), fold_df.deep_siamese_auc.mean()],
            "sd_auc": [fold_df.zero_shot_auc.std(ddof=1), fold_df.deep_siamese_auc.std(ddof=1)],
            "mean_auc_difference": [np.nan, fold_df.auc_difference.mean()],
            "sd_auc_difference": [np.nan, fold_df.auc_difference.std(ddof=1)],
            "mean_validation_auc": [np.nan, fold_df.best_validation_auc.mean()],
            "sd_validation_auc": [np.nan, fold_df.best_validation_auc.std(ddof=1)],
            "mean_test_pairs": [fold_df.n_test_pairs.mean(), fold_df.n_test_pairs.mean()],
            "mean_positive_test_pairs": [
                fold_df.n_positive_test_pairs.mean(),
                fold_df.n_positive_test_pairs.mean(),
            ],
            "mean_negative_test_pairs": [
                fold_df.n_negative_test_pairs.mean(),
                fold_df.n_negative_test_pairs.mean(),
            ],
            "embedding_dim": [np.nan, CFG.embedding_dim],
            "dropout": [np.nan, CFG.dropout],
            "learning_rate": [np.nan, CFG.learning_rate],
            "weight_decay": [np.nan, CFG.weight_decay],
        }
    )
    summary_df.to_csv(SUMMARY_FILE, index=False)
    print("\n===== Repeated 5x5 TOUGH-M1 Siamese FUSED summary =====", flush=True)
    print(summary_df.to_string(index=False), flush=True)
    print("\nSaved outputs:", flush=True)
    print(f"  {FOLD_RESULTS_FILE}", flush=True)
    print(f"  {SUMMARY_FILE}", flush=True)


def main() -> None:
    if os.environ.get("OVERWRITE", "0") == "1":
        for path in [FOLD_RESULTS_FILE, SUMMARY_FILE]:
            if os.path.exists(path):
                os.remove(path)

    for path, msg in [
        (FEATURE_FILE, "Run: Rscript scripts/TOUGH-M1/TOUGHM1_Export_RawFUSED_FeatureMatrix_thr4p8_15p0.R"),
        (PAIR_FILE, "Run: Rscript scripts/TOUGH-M1/TOUGHM1_FUSED_PairwiseAUC_FixedTend15.R"),
        (GROUP_FILE, "Run: python3 scripts/TOUGH-M1/TOUGHM1_Build_RCSB30_ClusterMapping.py"),
    ]:
        if not os.path.exists(path):
            raise SystemExit(f"Missing file: {path}\n{msg}")

    set_seed(CFG.seed)
    feature_df = pd.read_csv(FEATURE_FILE)
    group_df = pd.read_csv(GROUP_FILE)
    group_df = group_df[["code5", "rcsb30_cluster_id", "mapping_status"]].drop_duplicates()
    feature_df = feature_df.merge(group_df, on="code5", how="inner")

    ids = feature_df["code5"].to_numpy()
    id_to_group = dict(zip(feature_df["code5"], feature_df["rcsb30_cluster_id"]))
    features_raw = feature_df.drop(
        columns=["code5", "rcsb30_cluster_id", "mapping_status"]
    ).to_numpy(dtype=np.float32)
    id_to_row = {id_: i for i, id_ in enumerate(ids)}

    pairs = pd.read_csv(PAIR_FILE)
    pairs["label"] = pairs["label"].astype(np.int64)
    pairs = pairs[pairs["id1"].isin(id_to_row) & pairs["id2"].isin(id_to_row)].reset_index(drop=True)

    print(f"Loaded {len(ids)} raw FUSED pocket feature vectors.", flush=True)
    print(f"Input feature dimension: {features_raw.shape[1]}", flush=True)
    print(f"Sequence-cluster split groups: {len(set(id_to_group.values()))}", flush=True)
    print(f"Loaded {len(pairs)} TOUGH-M1 pairs.", flush=True)
    print(f"Positive pairs: {int((pairs.label == 1).sum())}", flush=True)
    print(f"Negative pairs: {int((pairs.label == 0).sum())}", flush=True)
    print(
        f"Fixed neural settings: embedding={CFG.embedding_dim}, dropout={CFG.dropout}, "
        f"lr={CFG.learning_rate}, weight_decay={CFG.weight_decay}",
        flush=True,
    )

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}", flush=True)

    done = completed_repeat_folds(FOLD_RESULTS_FILE)
    for repeat_id in range(CFG.n_repeats):
        repeat_seed = CFG.seed + 10000 * (repeat_id + 1)
        folds = make_group_folds(ids, id_to_group, CFG.k_folds, repeat_seed)
        for fold_id in range(CFG.k_folds):
            key = (repeat_id + 1, fold_id + 1)
            if key in done:
                print(f"Skipping completed repeat {key[0]}, fold {key[1]}.", flush=True)
                continue
            row = train_one_fold(
                repeat_id,
                fold_id,
                features_raw,
                pairs,
                folds,
                id_to_row,
                id_to_group,
                device,
            )
            append_csv_row(FOLD_RESULTS_FILE, row)
            done.add(key)
            write_summary()

    write_summary()


if __name__ == "__main__":
    raise SystemExit(
        "This file is a helper module. Run "
        "scripts/TOUGH-M1/TOUGHM1_FUSED_GroupShuffleSplit_TuneTmax_8_10_12_15.py "
        "for the final TOUGH-M1 analysis."
    )
