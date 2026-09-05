# FUSED: A Functional Representation for Joint Structural and Elemental Analysis of Protein Ligand Binding Sites

This repository contains the code, cached inputs, selected results, and figures used for the FUSED. FUSED represents each binding site as a threshold-dependent functional descriptor combining structural and  compositional features.

The repository is organized so that the main analyses can be rerun from cached descriptor objects, while raw structure data can be downloaded separately when full descriptor reconstruction is desired.

## Repository Structure

```text
FUSED-ProteinBindingSite/
├── data/
│   ├── EK/
│   │   └── Extended Kahraman Proteins Sorted/
│   ├── TOUGH-C1/
│   │   └── TOUGH_C1_protein_ligand_pairs.rds
│   ├── TOUGH-M1/
│   │   ├── TOUGH-M1_pocket.list
│   │   ├── TOUGH-M1_positive.list
│   │   ├── TOUGH-M1_negative.list
│   │   └── sequence_clusters/
│   └── cache/
├── figures/
│   ├── EK/
│   ├── TOUGH-C1/
│   └── TOUGH-M1/
├── results/
│   ├── EK/
│   ├── TOUGH-C1/
│   └── TOUGH-M1/
└── scripts/
    ├── EK/
    ├── TOUGH-C1/
    └── TOUGH-M1/
```

## Software Requirements

The R analyses were run with R and the main packages below:

```r
install.packages(c(
  "fda",
  "ranger",
  "glmnet",
  "MASS",
  "osfr"
))
```

The TOUGH-M1 Siamese neural-network analysis uses Python:

```bash
python3 -m pip install numpy pandas scikit-learn torch requests
```

## Data Sources

### Extended Kahraman Dataset

The Extended Kahraman dataset was introduced by Hoffmann et al. as an extension of the ligand binding-site dataset of Kahraman et al.

Relevant papers:

- Kahraman, Abdullah, et al. "Shape variation in protein binding pockets and their ligands." Journal of molecular biology 368.1 (2007): 283-301.
- Hoffmann et al. (2010), *A new protein binding pocket similarity measure based on comparison of clouds of atoms in 3D: application to ligand prediction*: https://pmc.ncbi.nlm.nih.gov/articles/PMC2838872/

The EK structures are PDB-derived. The cleaned repository expects the processed EK protein files under:

```text
data/EK/Extended Kahraman Proteins Sorted/
```

The cached FUSED descriptor object used by the main EK analyses is:

```text
data/cache/EK_raw_ILR_CDPA_thr4p8_20p0_by0p1.rds
```

### TOUGH-C1 Dataset

TOUGH-C1 is an established ligand binding-site classification benchmark introduced in DeepDrug3D paper.

Relevant paper:

- Pu, Limeng, et al. "DeepDrug3D: classification of ligand-binding pockets in proteins with a convolutional neural network." PLoS computational biology 15.2 (2019): e1006718.
Data set :
- https://osf.io/enz69/

The repository includes the processed protein-ligand pair object used by the final TOUGH-C1 nested cross-validation script:

```text
data/TOUGH-C1/TOUGH_C1_protein_ligand_pairs.rds
```

The cached FUSED descriptor object used by the final TOUGH-C1 nested cross-validation analysis is:

```text
data/cache/TOUGH_C1_LBS_ligand_ILR_CDPA_to20A.rds
```

If rebuilding TOUGH-C1 descriptors or rerunning the element-composition audit from raw structures, place the raw sorted files under:

```text
data/TOUGH-C1/TOUGH-C1_Proteins_Sorted/
data/TOUGH-C1/TOUGH-C1_Ligands_Sorted/
```

The script below can be used to download, organize, and build the TOUGH-C1 protein-ligand pair object:

```bash
Rscript scripts/TOUGH-C1/TOUGHC1_Download_Organize_BuildPairs.R
```

### TOUGH-M1 Dataset

TOUGH-M1 is a pairwise pocket-matching benchmark introduced by Govindaraj and Brylinski.

Relevant paper:

- Govindaraj and Brylinski (2018), *Comparative assessment of strategies to identify similar ligand-binding pockets in proteins*: https://pmc.ncbi.nlm.nih.gov/articles/PMC5845264/

The TOUGH-M1 positive, negative, and pocket-list files are expected at:

```text
data/TOUGH-M1/TOUGH-M1_pocket.list
data/TOUGH-M1/TOUGH-M1_positive.list
data/TOUGH-M1/TOUGH-M1_negative.list
```

For convenience, TOUGH-M1 data can also be obtained through the DeeplyTough repository and associated dataset resources:

- DeeplyTough GitHub: https://github.com/BenevolentAI/DeeplyTough
- DeeplyTough/TOUGH-M1 dataset resource: https://zenodo.org/records/3687316

The TOUGH-M1 sequence-cluster evaluation uses current RCSB PDB 30% sequence-identity clusters.

RCSB sequence-cluster file:

- https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-30.txt



## Cached Descriptor Files

The main analyses use cached descriptor objects stored in:

```text
data/cache/
```

Key cached files include:

```text
EK_raw_ILR_CDPA_thr4p8_20p0_by0p1.rds
EK_sulfur_present_thr4p8_20p0_by0p1.rds
TOUGH_C1_LBS_ligand_ILR_CDPA_to20A.rds
TOUGH_C1_sulfur_present_at10A.rds
TOUGHM1_raw_ILR_CDPA_thr4p8_15p0_by0p1.rds
TOUGHM1_rawFUSED_features_thr4p8_15p0.csv
```



## Running the Main Analyses

Run all commands from the repository root.

### EK Nested Cross-Validation

```bash
Rscript scripts/EK/ExtendedKahraman_RepeatedNestedCV_TuneTend_MaxInner_to20A_by1p0_parallel.R
```

### TOUGH-C1 Nested Cross-Validation

```bash
Rscript scripts/TOUGH-C1/TOUGHC1_RepeatedNestedCV_TuneTend_MaxInner_to20A_by1p0_parallel_Resumable.R
```

To force a complete rerun from scratch:

```bash
RR_START_FROM_SCRATCH=1 Rscript scripts/TOUGH-C1/TOUGHC1_RepeatedNestedCV_TuneTend_MaxInner_to20A_by1p0_parallel_Resumable.R
```

### TOUGH-M1 Sequence-Cluster Mapping

```bash
python3 scripts/TOUGH-M1/TOUGHM1_Build_RCSB30_ClusterMapping.py
```

### TOUGH-M1 Descriptor Construction

```bash
Rscript scripts/TOUGH-M1/TOUGHM1_Build_RawFUSED_Descriptors_thr4p8_15p0.R
```

### TOUGH-M1 Feature Matrix Export

```bash
Rscript scripts/TOUGH-M1/TOUGHM1_Export_RawFUSED_FeatureMatrix_thr4p8_15p0.R
```

### TOUGH-M1 Pairwise Pocket-Matching Analysis

```bash
python3 scripts/TOUGH-M1/TOUGHM1_FUSED_GroupShuffleSplit_TuneTmax_8_10_12_15.py
```

The TOUGH-M1 script is resumable. If all splits have already been completed, it will load the completed split outputs and regenerate the summary quickly. To rerun the full analysis from scratch:

```bash
OVERWRITE=1 python3 scripts/TOUGH-M1/TOUGHM1_FUSED_GroupShuffleSplit_TuneTmax_8_10_12_15.py
```

## Supporting Analyses

Additional EK scripts reproduce feature-set comparisons, smoothing sensitivity, variance-retention sensitivity, sulfur sensitivity, and exploratory figures:

```text
scripts/EK/
```

Additional TOUGH-C1 scripts reproduce sulfur and element-composition sensitivity analyses:

```text
scripts/TOUGH-C1/
```

Outputs are written to:

```text
results/EK/
results/TOUGH-C1/
results/TOUGH-M1/
figures/EK/
figures/TOUGH-C1/
figures/TOUGH-M1/
```


