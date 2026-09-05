# ================================================================
# Download TOUGH-C1 from OSF (node "enz69")
# ================================================================

rm(list = ls()); invisible(gc())
options(stringsAsFactors = FALSE)

# 1) Install & load osfr -----------------(Data set is in Open Science Framework (OSF))( https://osf.io/enz69/.).osfr package allows R to communicate with the Open Science Framework (OSF
if (!requireNamespace("osfr", quietly = TRUE)) {
  install.packages("osfr")
}
library(osfr)

# 2) Clean old downloads  ---------------
raw_root <- file.path("data", "TOUGH-C1_raw")

if (dir.exists(raw_root)) {
  cat("Removing existing TOUGH-C1_raw directory...\n")
  unlink(raw_root, recursive = TRUE, force = TRUE)
}
dir.create(raw_root, showWarnings = FALSE, recursive = TRUE)

# 3) Retrieve OSF project node ------------------------------------
# TOUGH-C1 OSF: https://osf.io/enz69/
cat("Retrieving OSF node 'enz69' (TOUGH-C1)...\n")
proj  <- osf_retrieve_node("enz69")

# 4) List ALL files in the project --------------------------------
# n_max = Inf to avoid truncation
cat("Listing all files in the project...\n")
files <- osf_ls_files(proj, n_max = Inf)

print(files[, c("name", "id")])

# 5) Download everything into TOUGH-C1_raw ------------------------
cat("Downloading all files to", raw_root, "...\n")
osf_download(
  files,
  path      = raw_root,
  conflicts = "overwrite"   # overwrite if you re-run
)

cat("\nDownload finished. Top-level contents of TOUGH-C1_raw:\n")
print(list.files(raw_root))

# 6) Unpack ALL .tar.gz archives ---------------------------------
tgz_paths <- list.files(
  raw_root,
  pattern   = "\\.tar\\.gz$",
  full.names = TRUE
)

cat("\nFound", length(tgz_paths), "tar.gz files:\n")
print(basename(tgz_paths))

for (f in tgz_paths) {
  out_dir <- sub("\\.tar\\.gz$", "", f)   # e.g. "TOUGH-C1_raw/fpocket-control"
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Unpacking", basename(f), "->", out_dir, "\n")
  untar(f, exdir = out_dir)
}

cat("\nAfter unpacking, top-level contents of TOUGH-C1_raw:\n")
print(list.files(raw_root))

# Optional peeks ----------------------------------------------
if (dir.exists(file.path(raw_root, "protein-nucleotide-pdbqt"))) {
  cat("\nExample files in protein-nucleotide-pdbqt:\n")
  print(head(list.files(file.path(raw_root, "protein-nucleotide-pdbqt"),
                        recursive = TRUE), 20))
}

if (dir.exists(file.path(raw_root, "protein-heme-pdbqt"))) {
  cat("\nExample files in protein-heme-pdbqt:\n")
  print(head(list.files(file.path(raw_root, "protein-heme-pdbqt"),
                        recursive = TRUE), 20))
}

if (dir.exists(file.path(raw_root, "fpocket-control"))) {
  cat("\nExample files in fpocket-control:\n")
  print(head(list.files(file.path(raw_root, "fpocket-control"),
                        recursive = TRUE), 20))
}

if (dir.exists(file.path(raw_root, "ligand-nucleotide-pdbqt"))) {
  cat("\nExample files in ligand-nucleotide-pdbqt:\n")
  print(head(list.files(file.path(raw_root, "ligand-nucleotide-pdbqt"),
                        recursive = TRUE), 20))
}

if (dir.exists(file.path(raw_root, "ligand-heme-pdbqt"))) {
  cat("\nExample files in ligand-heme-pdbqt:\n")
  print(head(list.files(file.path(raw_root, "ligand-heme-pdbqt"),
                        recursive = TRUE), 20))
}

if (dir.exists(file.path(raw_root, "ligand-control-pdbqt"))) {
  cat("\nExample files in ligand-control-pdbqt:\n")
  print(head(list.files(file.path(raw_root, "ligand-control-pdbqt"),
                        recursive = TRUE), 20))
}




#########################################################
##########################################################
##########################################################


# ================================================================
# Build sorted folders for TOUGH-C1:
#   Proteins: NUC / HEME / CONTROL
#   Ligands : NUC / HEME / CONTROL
# ================================================================

rm(list = ls()); invisible(gc())

prot_root <- file.path("data", "TOUGH-C1_Proteins_Sorted")
lig_root  <- file.path("data", "TOUGH-C1_Ligands_Sorted")

dir.create(prot_root, showWarnings = FALSE)
dir.create(lig_root,  showWarnings = FALSE)

# Create class subfolders for proteins and ligands
for (cls in c("NUC", "HEME", "CONTROL")) {
  dir.create(file.path(prot_root, cls), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(lig_root,  cls), showWarnings = FALSE, recursive = TRUE)
}

# Generic helper to copy all .pdbqt / .pdb files
copy_structures <- function(src_dir, dst_dir,
                            pattern = "\\.(pdbqt|pdb)$",
                            label   = "") {
  if (!dir.exists(src_dir)) {
    stop("Source directory not found: ", src_dir,
         if (nzchar(label)) paste0(" (", label, ")"))
  }
  files <- list.files(src_dir, pattern = pattern,
                      full.names = TRUE, recursive = TRUE)
  cat("Found", length(files), "files in", src_dir,
      if (nzchar(label)) paste0(" [", label, "]"), "\n")
  if (length(files) == 0) return(invisible(0))
  res <- file.copy(files, dst_dir, overwrite = TRUE)
  cat("Copied", sum(res), "files to", dst_dir, "\n\n")
  invisible(sum(res))
}

# ---------------------- Proteins ---------------------------------

# 1) NUC proteins
n_nuc_prot <- copy_structures(
  src_dir = file.path(raw_root, "protein-nucleotide-pdbqt"),
  dst_dir = file.path(prot_root, "NUC"),
  label   = "NUC proteins"
)

# 2) HEME proteins
n_heme_prot <- copy_structures(
  src_dir = file.path(raw_root, "protein-heme-pdbqt"),
  dst_dir = file.path(prot_root, "HEME"),
  label   = "HEME proteins"
)

# 3) CONTROL proteins
n_ctrl_prot <- copy_structures(
  src_dir = file.path(raw_root, "protein-control-pdbqt"),
  dst_dir = file.path(prot_root, "CONTROL"),
  label   = "CONTROL proteins"
)

# ---------------------- Ligands ----------------------------------

# ---------------------- Ligands ----------------------------------

# 4) NUC ligands  (nucleotide)
n_nuc_lig <- copy_structures(
  src_dir = file.path(raw_root, "ligand-nucleotide"),  # <- changed
  dst_dir = file.path(lig_root, "NUC"),
  label   = "NUC ligands"
)

# 5) HEME ligands
n_heme_lig <- copy_structures(
  src_dir = file.path(raw_root, "ligand-heme"),        # <- changed
  dst_dir = file.path(lig_root, "HEME"),
  label   = "HEME ligands"
)

# 6) CONTROL ligands
n_ctrl_lig <- copy_structures(
  src_dir = file.path(raw_root, "ligand-control"),     # <- changed
  dst_dir = file.path(lig_root, "CONTROL"),
  label   = "CONTROL ligands"
)

# ---------------------- Summary ----------------------------------

cat("Summary (proteins):\n")
cat("  NUC     proteins:", n_nuc_prot, "\n")
cat("  HEME    proteins:", n_heme_prot, "\n")
cat("  CONTROL proteins:", n_ctrl_prot, "\n\n")

cat("Summary (ligands):\n")
cat("  NUC     ligands :", n_nuc_lig, "\n")
cat("  HEME    ligands :", n_heme_lig, "\n")
cat("  CONTROL ligands :", n_ctrl_lig, "\n\n")

cat("Check protein folders:\n")
print(sapply(c("NUC","HEME","CONTROL"), function(cls) {
  length(list.files(file.path(prot_root, cls)))
}))

cat("\nCheck ligand folders:\n")
print(sapply(c("NUC","HEME","CONTROL"), function(cls) {
  length(list.files(file.path(lig_root, cls)))
}))


rm(list = ls()); invisible(gc())

prot_root <- file.path("data", "TOUGH-C1_Proteins_Sorted")
lig_root  <- file.path("data", "TOUGH-C1_Ligands_Sorted")

classes <- c("NUC", "HEME", "CONTROL")

# Helper to build pairs for one class ------------------------------
build_pairs_for_class <- function(cls,
                                  prot_root,
                                  lig_root,
                                  pick_one_ligand = TRUE) {
  p_dir <- file.path(prot_root, cls)
  l_dir <- file.path(lig_root,  cls)
  
  # list all protein & ligand files
  p_files <- list.files(p_dir, full.names = TRUE)
  l_files <- list.files(l_dir, full.names = TRUE)
  
  p_base <- tools::file_path_sans_ext(basename(p_files))
  l_base <- tools::file_path_sans_ext(basename(l_files))
  
  # ligands: split into core ID (protein ID) + 2-digit suffix
  lig_core   <- substr(l_base, 1, nchar(l_base) - 2)
  lig_suffix <- substr(l_base, nchar(l_base) - 1, nchar(l_base))
  
  lig_df <- data.frame(
    core   = lig_core,
    suffix = lig_suffix,
    l_base = l_base,
    l_file = l_files,
    stringsAsFactors = FALSE
  )
  lig_df$suffix_num <- suppressWarnings(as.integer(lig_df$suffix))
  
  # protein table
  prot_df <- data.frame(
    p_base = p_base,
    p_file = p_files,
    stringsAsFactors = FALSE
  )
  
  # merge by core ID (protein base)
  merged <- merge(
    prot_df,
    lig_df,
    by.x = "p_base",
    by.y = "core",
    all.x = TRUE,
    sort  = FALSE
  )
  
  # If multiple ligands per protein -> either keep all or pick one
  if (pick_one_ligand) {
    # For each protein, keep the ligand with the smallest suffix_num
    merged <- merged[order(merged$p_base, merged$suffix_num), ]
    merged <- merged[!duplicated(merged$p_base), ]
  }
  
  # Add class column
  merged$class <- cls
  
  # Reorder columns nicely
  merged <- merged[, c("class", "p_base", "p_file",
                       "l_base", "l_file", "suffix", "suffix_num")]
  merged
}

# Build pairs for all three classes -------------------------------
pairs_list <- lapply(classes, build_pairs_for_class,
                     prot_root = prot_root,
                     lig_root  = lig_root,
                     pick_one_ligand = TRUE)  # set FALSE if you want all ligands

paired_tbl <- do.call(rbind, pairs_list)

# Summary ----------------------------------------------------------
cat("Total protein–ligand pairs:", nrow(paired_tbl), "\n\n")

for (cls in classes) {
  sub <- paired_tbl[paired_tbl$class == cls, ]
  n_prot <- length(unique(sub$p_base))
  n_with_lig <- sum(!is.na(sub$l_file))
  cat("Class:", cls, "\n")
  cat("  Unique proteins:         ", n_prot, "\n")
  cat("  Proteins with a ligand:  ", n_with_lig, "\n")
  cat("  Example rows:\n")
  print(head(sub, 3))
  cat("\n")
}

# Save the mapping for downstream analysis
saveRDS(
  paired_tbl,
  file = file.path("data", "TOUGH-C1", "TOUGH_C1_protein_ligand_pairs.rds")
)

cat("Saved mapping to TOUGH_C1_protein_ligand_pairs.rds\n")
