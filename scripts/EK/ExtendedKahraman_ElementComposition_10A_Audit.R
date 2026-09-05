# =====================================================================
# Extended Kahraman | Element composition audit at 10 Angstrom
#   - Saves all element summaries, not only sulfur/phosphorus.
# =====================================================================

rm(list = ls()); invisible(gc())

# ------------------- paths -------------------
cache_file <- file.path("data", "cache", "EK_raw_ILR_CDPA_thr4p8_20p0_by0p1.rds")
out_dir <- file.path("results", "EK")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cache_file)) {
  stop("Descriptor cache not found. Run the nested-CV script first to create: ", cache_file)
}

raw <- readRDS(cache_file)
files_df <- raw$files_df
labels <- raw$labels

# ------------------- settings -------------------
dthr <- 10.0
core_elements <- c("C", "O", "N")
elements_to_force_include <- c("S", "P")

ligand_aliases <- list(
  AMP = c("AMP"),
  ATP = c("ATP", "AGS"),
  FAD = c("FAD"),
  FMN = c("FMN"),
  GLC = c("GLC", "IMD", "LAT", "NAG", "BGC"),
  HEM = c("HEM"),
  NAD = c("NAD", "NAI"),
  PO4 = c("PO4")
)

# =====================================================================
# helpers
# =====================================================================
resolve_existing_path <- function(path_value) {
  candidates <- unique(c(
    path_value,
    file.path("data", "EK", sub("^data/", "", path_value)),
    file.path("data", path_value)
  ))
  hit <- candidates[file.exists(candidates)][1]
  ifelse(is.na(hit), candidates[1], hit)
}

read_pdb_minimal <- function(pdb_path) {
  ln <- readLines(pdb_path, warn = FALSE)
  ln <- ln[grepl("^(ATOM  |HETATM)", ln)]
  if (length(ln) == 0) return(NULL)

  rec   <- substr(ln, 1, 6)
  resn  <- trimws(substr(ln, 18, 20))
  chain <- substr(ln, 22, 22)
  x     <- as.numeric(substr(ln, 31, 38))
  y     <- as.numeric(substr(ln, 39, 46))
  z     <- as.numeric(substr(ln, 47, 54))
  elem_raw <- trimws(substr(ln, 77, 78))
  aname    <- trimws(substr(ln, 13, 16))
  elem <- ifelse(elem_raw == "", toupper(substr(aname, 1, 1)), toupper(elem_raw))

  keepH <- elem != "H"
  data.frame(
    rec = rec[keepH], resn = resn[keepH], chain = chain[keepH],
    x = x[keepH], y = y[keepH], z = z[keepH], elem = elem[keepH],
    stringsAsFactors = FALSE
  )
}

get_site_atoms <- function(atom_df, ligand_vec, dthr, min_atoms = 10) {
  lig <- atom_df[atom_df$resn %in% ligand_vec, , drop = FALSE]
  if (nrow(lig) == 0) return(NULL)

  lig_chain <- lig$chain[1]
  prot <- atom_df[
    atom_df$rec == "ATOM  " &
      atom_df$chain == lig_chain &
      !(atom_df$resn %in% ligand_vec),
    , drop = FALSE
  ]
  if (nrow(prot) == 0) return(NULL)

  pick <- logical(nrow(prot))
  for (i in seq_len(nrow(lig))) {
    dx <- prot$x - lig$x[i]
    dy <- prot$y - lig$y[i]
    dz <- prot$z - lig$z[i]
    pick <- pick | (sqrt(dx*dx + dy*dy + dz*dz) <= dthr)
  }

  site <- prot[pick, , drop = FALSE]
  if (nrow(site) < min_atoms) return(NULL)
  site
}

five_number_named <- function(x, suffix = "") {
  qs <- as.numeric(quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
  names(qs) <- paste0(c("min", "q1", "median", "q3", "max"), suffix)
  qs
}

summarize_counts <- function(g) {
  all_counts <- g$count
  present_counts <- g$count[g$count > 0]
  present_five <- if (length(present_counts) > 0) {
    five_number_named(present_counts, "_among_present")
  } else {
    out <- rep(NA_real_, 5)
    names(out) <- paste0(c("min", "q1", "median", "q3", "max"), "_among_present")
    out
  }

  data.frame(
    n_sites = nrow(g),
    n_sites_present = sum(g$count > 0),
    pct_sites_present = 100 * mean(g$count > 0),
    mean_count_all_sites = mean(all_counts),
    mean_count_among_present = ifelse(length(present_counts) > 0, mean(present_counts), NA_real_),
    as.list(five_number_named(all_counts)),
    as.list(present_five),
    stringsAsFactors = FALSE
  )
}

# =====================================================================
# collect site-level element counts
# =====================================================================
cat("Auditing protein binding-site element composition at", dthr, "Angstrom\n")

raw_count_records <- list()
site_total_records <- list()

for (i in seq_len(nrow(files_df))) {
  fpath <- resolve_existing_path(files_df$path[i])
  lbl <- as.character(labels[i])
  aliases <- ligand_aliases[[lbl]]
  adf <- tryCatch(read_pdb_minimal(fpath), error = function(e) NULL)
  if (is.null(adf)) next

  site <- get_site_atoms(adf, aliases, dthr, min_atoms = 10)
  if (is.null(site)) next

  site_id <- tools::file_path_sans_ext(basename(fpath))
  tab <- table(site$elem)

  for (el in names(tab)) {
    raw_count_records[[length(raw_count_records) + 1L]] <- data.frame(
      site_index = i,
      site_id = site_id,
      label = lbl,
      element = el,
      count = as.integer(tab[[el]]),
      stringsAsFactors = FALSE
    )
  }

  site_total_records[[length(site_total_records) + 1L]] <- data.frame(
    site_index = i,
    site_id = site_id,
    label = lbl,
    total_atoms = nrow(site),
    stringsAsFactors = FALSE
  )

  if (i %% 100 == 0) cat("processed", i, "of", nrow(files_df), "\n")
}

raw_counts <- do.call(rbind, raw_count_records)
site_totals <- do.call(rbind, site_total_records)

all_elements <- sort(unique(c(raw_counts$element, elements_to_force_include)))
all_sites_elements <- expand.grid(
  site_index = site_totals$site_index,
  element = all_elements,
  stringsAsFactors = FALSE
)
all_sites_elements <- merge(all_sites_elements, site_totals, by = "site_index")
counts_df <- merge(
  all_sites_elements,
  raw_counts[, c("site_index", "element", "count")],
  by = c("site_index", "element"),
  all.x = TRUE
)
counts_df$count[is.na(counts_df$count)] <- 0L
counts_df$proportion <- counts_df$count / counts_df$total_atoms
counts_df$threshold <- dthr
counts_df <- counts_df[
  order(counts_df$site_index, counts_df$element),
  c("site_index", "site_id", "label", "threshold", "element", "count", "total_atoms", "proportion")
]

write.csv(counts_df, file.path(out_dir, "EK_element_counts_10A_by_site.csv"), row.names = FALSE)

# =====================================================================
# summaries for all elements and by class
# =====================================================================
element_summary <- do.call(rbind, lapply(split(counts_df, counts_df$element), function(g) {
  cbind(
    data.frame(element = g$element[1], stringsAsFactors = FALSE),
    summarize_counts(g)
  )
}))
element_summary <- element_summary[order(element_summary$element), ]

by_class_key <- interaction(counts_df$label, counts_df$element, drop = TRUE, sep = "|")
class_element_summary <- do.call(rbind, lapply(split(counts_df, by_class_key), function(g) {
  cbind(
    data.frame(label = g$label[1], element = g$element[1], stringsAsFactors = FALSE),
    summarize_counts(g)
  )
}))
class_element_summary <- class_element_summary[
  order(class_element_summary$label, class_element_summary$element),
]

site_non_con <- do.call(rbind, lapply(split(counts_df, counts_df$site_index), function(g) {
  non_con <- g[!(g$element %in% core_elements), , drop = FALSE]
  data.frame(
    site_index = g$site_index[1],
    site_id = g$site_id[1],
    label = g$label[1],
    threshold = dthr,
    total_atoms = g$total_atoms[1],
    n_non_CON = sum(non_con$count),
    pct_non_CON = 100 * sum(non_con$count) / g$total_atoms[1],
    elements_non_CON_present = paste(non_con$element[non_con$count > 0], collapse = ";"),
    stringsAsFactors = FALSE
  )
}))

pct_non_con_five <- five_number_named(site_non_con$pct_non_CON, "_pct_non_CON")
site_non_con_summary <- data.frame(
  threshold = dthr,
  n_sites = nrow(site_non_con),
  n_sites_with_any_non_CON = sum(site_non_con$n_non_CON > 0),
  pct_sites_with_any_non_CON = 100 * mean(site_non_con$n_non_CON > 0),
  mean_pct_non_CON = mean(site_non_con$pct_non_CON),
  as.list(pct_non_con_five),
  stringsAsFactors = FALSE
)

write.csv(element_summary, file.path(out_dir, "EK_element_summary_10A_all_elements.csv"), row.names = FALSE)
write.csv(class_element_summary, file.path(out_dir, "EK_element_summary_10A_by_class.csv"), row.names = FALSE)
write.csv(site_non_con, file.path(out_dir, "EK_non_CON_elements_10A_by_site.csv"), row.names = FALSE)
write.csv(site_non_con_summary, file.path(out_dir, "EK_non_CON_summary_10A.csv"), row.names = FALSE)

cat("\n===== Elements observed or explicitly checked at 10 Angstrom =====\n")
print(element_summary[, c(
  "element", "n_sites", "n_sites_present", "pct_sites_present",
  "mean_count_all_sites", "mean_count_among_present",
  "min", "q1", "median", "q3", "max"
)])

cat("\n===== Non-C/O/N site-level summary at 10 Angstrom =====\n")
print(site_non_con_summary)

cat("\nSaved 10A element audit outputs under:", out_dir, "\n")
