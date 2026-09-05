#!/usr/bin/env python3
"""
Build TOUGH-M1 -> RCSB 30% sequence-cluster mapping.

The current RCSB sequence-cluster file is entity-based
(`clusters-by-entity-30.txt`), so this script maps TOUGH-M1 code5 IDs
such as 1lkxD to polymer entity IDs using the RCSB Data API. Unmapped entries are marked with singleton identifiers so they can be excluded.
"""

from __future__ import annotations

import json
import os
import time
from typing import Iterable

import pandas as pd
import requests


ROOT = os.getcwd()
INSPECT_DIR = os.path.join(ROOT, "data", "TOUGH-M1")
CLUSTER_DIR = os.path.join(INSPECT_DIR, "sequence_clusters")
POCKET_FILE = os.path.join(INSPECT_DIR, "TOUGH-M1_pocket.list")
CLUSTER_FILE = os.path.join(CLUSTER_DIR, "clusters-by-entity-30.txt")
OUT_FILE = os.path.join(CLUSTER_DIR, "TOUGHM1_code5_rcsb30_cluster_mapping.csv")

RCSB_GRAPHQL = "https://data.rcsb.org/graphql"
BATCH_SIZE = 250


def chunks(xs: list[str], n: int) -> Iterable[list[str]]:
    for i in range(0, len(xs), n):
        yield xs[i : i + n]


def load_entity_to_cluster(path: str) -> dict[str, str]:
    entity_to_cluster: dict[str, str] = {}
    with open(path, "r", encoding="utf-8") as f:
        for cluster_idx, line in enumerate(f, start=1):
            members = line.strip().split()
            cluster_id = f"RCSB30_{cluster_idx}"
            for m in members:
                entity_to_cluster[m.upper()] = cluster_id
    return entity_to_cluster


def query_instance_to_entity(instance_ids: list[str]) -> dict[str, str]:
    query = """
    query($ids: [String!]!) {
      polymer_entity_instances(instance_ids: $ids) {
        rcsb_id
        polymer_entity {
          rcsb_id
        }
      }
    }
    """
    response = requests.post(
        RCSB_GRAPHQL,
        json={"query": query, "variables": {"ids": instance_ids}},
        timeout=60,
    )
    response.raise_for_status()
    payload = response.json()
    if "errors" in payload:
        raise RuntimeError(json.dumps(payload["errors"], indent=2))

    out: dict[str, str] = {}
    for item in payload["data"]["polymer_entity_instances"] or []:
        if item and item.get("polymer_entity"):
            out[item["rcsb_id"].upper()] = item["polymer_entity"]["rcsb_id"].upper()
    return out


def main() -> None:
    if not os.path.exists(POCKET_FILE):
        raise SystemExit(f"Missing pocket list: {POCKET_FILE}")
    if not os.path.exists(CLUSTER_FILE):
        raise SystemExit(f"Missing cluster file: {CLUSTER_FILE}")

    pockets = pd.read_csv(
        POCKET_FILE,
        sep=r"\s+",
        header=None,
        names=["code5", "selected_fpocket_number", "overlap_score"],
    )
    pockets["pdb_id"] = pockets["code5"].str[:4].str.upper()
    pockets["chain_id"] = pockets["code5"].str[4:5]
    pockets["instance_id"] = pockets["pdb_id"] + "." + pockets["chain_id"]

    entity_to_cluster = load_entity_to_cluster(CLUSTER_FILE)
    print(f"Loaded {len(entity_to_cluster)} entity-to-cluster memberships.")

    unique_instances = sorted(pockets["instance_id"].unique())
    instance_to_entity: dict[str, str] = {}

    for batch_id, batch in enumerate(chunks(unique_instances, BATCH_SIZE), start=1):
      # RCSB can occasionally throttle. A short pause keeps this polite.
        got = query_instance_to_entity(batch)
        instance_to_entity.update(got)
        print(
            f"Mapped batch {batch_id}: {len(instance_to_entity)} / {len(unique_instances)} instances",
            flush=True,
        )
        time.sleep(0.1)

    pockets["rcsb_entity_id"] = pockets["instance_id"].map(instance_to_entity)
    pockets["rcsb30_cluster_id"] = pockets["rcsb_entity_id"].map(entity_to_cluster)
    pockets["mapping_status"] = "mapped"

    missing_entity = pockets["rcsb_entity_id"].isna()
    missing_cluster = pockets["rcsb30_cluster_id"].isna() & ~missing_entity
    pockets.loc[missing_entity, "mapping_status"] = "missing_entity"
    pockets.loc[missing_cluster, "mapping_status"] = "missing_cluster"

    singleton = pockets["rcsb30_cluster_id"].isna()
    pockets.loc[singleton, "rcsb30_cluster_id"] = (
        "SINGLETON_" + pockets.loc[singleton, "code5"].astype(str)
    )

    pockets.to_csv(OUT_FILE, index=False)

    print("\n===== TOUGH-M1 RCSB30 cluster mapping summary =====")
    print(pockets["mapping_status"].value_counts(dropna=False).to_string())
    print("Unique split groups:", pockets["rcsb30_cluster_id"].nunique())
    print("Output:", OUT_FILE)


if __name__ == "__main__":
    main()
