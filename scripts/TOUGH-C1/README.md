# TOUGH-C1 Three-Class Classification

We train each FUSED classifier on CONTROL, HEME,
and NUC together. HEME and NUC metrics are then calculated separately
against the remaining classes.

Run commands from the repository root. R packages: `fda`, `ranger`,
`glmnet`, and `MASS` (plus the standard `parallel` package).

## Full Cohort: 4,002 Pockets

```bash
Rscript scripts/TOUGH-C1/TOUGHC1_RF_ThreeClass_RepeatedNestedCV.R
Rscript scripts/TOUGH-C1/TOUGHC1_MLR_NMM_ThreeClass_RepeatedNestedCV.R
```

The second script combines its MLR-l1/NMM results with RF results from the
first script. Outputs go to `results/TOUGH-C1/threeclass_full4002/RF/`
and `results/TOUGH-C1/threeclass_full4002/all_classifiers/`.

The cached input is `data/cache/TOUGH_C1_LBS_ligand_ILR_CDPA_to20A.rds`.
The saved outer-fold CSV in `folds/` preserves the assignments used in the
DeepDrug3D comparison: each pocket is tested once in each of five repeats.




## Common Cohort: 4,000 Pockets

For the BindSiteS-CNN comparison, exclude `1ukmB` and `2o8cD`, whose required
precomputed pocket representations were unavailable. The retained pockets
keep their original outer-fold assignments. Run both scripts with:

```bash
export FUSED_ALLOW_SUBSET=1
export FUSED_FOLD_FILE=scripts/TOUGH-C1/folds/TOUGHC1_threeclass_5x5_outer_fold_assignments_common4000.csv
Rscript scripts/TOUGH-C1/TOUGHC1_RF_ThreeClass_RepeatedNestedCV.R
Rscript scripts/TOUGH-C1/TOUGHC1_MLR_NMM_ThreeClass_RepeatedNestedCV.R
unset FUSED_ALLOW_SUBSET FUSED_FOLD_FILE
```

Outputs are separated under `results/TOUGH-C1/threeclass_common4000/`.
Neural-baseline training code and HPCC environments are not included in
this focused FUSED update.

## Auditing and Resuming

```bash
FUSED_AUDIT_ONLY=1 Rscript scripts/TOUGH-C1/TOUGHC1_RF_ThreeClass_RepeatedNestedCV.R
FUSED_AUDIT_ONLY=1 FUSED_SKIP_RF_IMPORT=1 Rscript scripts/TOUGH-C1/TOUGHC1_MLR_NMM_ThreeClass_RepeatedNestedCV.R
Rscript scripts/TOUGH-C1/test_threeclass.R
```

Audits check fold membership, class labels, and protein/ligand identities
against the cached descriptor rows without fitting models. The subset
audit uses the same environment settings as the subset run above.

Completed folds are resumed from checkpoints. Input/source signatures
prevent reuse after scripts, descriptors, or fold assignments change.
Use a fresh `FUSED_OUTPUT_DIR` to rerun rather than deleting old results.
`FUSED_N_CORES` controls parallel inner evaluation (default 4).
`FUSED_REPETITION` and `FUSED_OUTER_FOLD` optionally select indices 1-5.
When running individual folds in separate output directories, aggregate
them before reporting a full 25-fold summary. `FUSED_SKIP_RF_IMPORT=1`
runs MLR-l1/NMM without requiring or importing RF results.
`FUSED_RF_RESULTS` overrides the RF metrics CSV path; use results from
the same cohort and requested folds. `FUSED_DESCRIPTOR_FILE`,
`FUSED_FOLD_FILE`, and `FUSED_OUTPUT_DIR` override input/output paths.

## Reference Results

`results/TOUGH-C1/threeclass_reference/full4002/` contains the original
full-cohort outer metrics, aggregate summary, and manuscript-formatted
table. These CSVs were copied from the submitted analysis, not regenerated
by this repository update. Seeds, model settings, and computational helper
functions are retained from that analysis. Package versions and numerical
libraries can affect reproducibility. Large predictions, checkpoints,
prepared neural inputs, and unrelated experiments are not published here.
