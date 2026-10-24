# TOUGH-C1 Three-Class Classification

We train each FUSED classifier on CONTROL, HEME,
and NUC together. HEME and NUC metrics are then calculated separately
against the remaining classes.

Run commands from the repository root. R packages: `fda`, `ranger`,
`glmnet`, and `MASS`.

## Full Cohort: 4,002 Pockets

```bash
Rscript scripts/TOUGH-C1/TOUGHC1_RF_ThreeClass_RepeatedNestedCV.R
Rscript scripts/TOUGH-C1/TOUGHC1_MLR_NMM_ThreeClass_RepeatedNestedCV.R
```

 Outputs go to `results/TOUGH-C1/threeclass_full4002/RF/`
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

