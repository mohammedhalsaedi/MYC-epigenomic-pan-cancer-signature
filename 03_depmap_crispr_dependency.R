# ==============================================================================
# Script: 03_depmap_crispr_dependency.R
# Purpose: Functional assessment of gene essentiality using DepMap CRISPR-Cas9
#          knockout data (Chronos scores) across breast, colon, and liver lineages.
# Inputs:  CRISPRGeneEffect.csv, CRISPRGeneDependency.csv, 
#          CRISPRInferredCommonEssentials.csv, and Model.csv.
# Outputs: Line-level and cancer-specific dependency statistics, waterfall plots, 
#          and selectivity summaries.
# ==============================================================================

# SECTION 0 - INSTALL / LOAD PACKAGES
# =============================================================================
# install.packages(c("tidyverse", "openxlsx", "ggplot2", "ggrepel",
#                    "patchwork", "pheatmap", "RColorBrewer"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
})

# =============================================================================
# SECTION 1 - CONFIGURATION
# =============================================================================

# -- 1a. File paths ------------------------------------------------------------
data_dir <- "~/Documents/DepMap"   # same folder used for the CRISPR analysis

# AllGenes release (log2(TPM+1), wide: ModelID x gene) -- includes non-coding genes
expression_file <- file.path(data_dir, "OmicsExpressionTPMLogp1HumanAllGenes.csv")

model_file <- file.path(data_dir, "Model.csv")   # same metadata file as CRISPR script

# -- 1b. Output directories ----------------------------------------------------
results_dir <- file.path(data_dir, "R2_CCLE_results")
figures_dir <- file.path(results_dir, "Figures")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# -- 1c. Genes of interest (all 15 - AllGenes release covers lncRNA/snoRNA too) --
GOI <- c("MYC", "PVT1", "CASC11", "BOP1", "SNHG17", "BYSL",
         "CHD7", "MLLT10", "FOXK1", "SNHG16", "SNORA71B",
         "SLC7A1", "ZBTB10", "PRKAR1B", "FCHSD1")

cat("Genes of interest (n =", length(GOI), "):\n")
cat(" ", paste(GOI, collapse = ", "), "\n\n")

# -- 1d. Cancer-type definitions (Oncotree-based, NOT fixed cell-line lists) --
# Each cancer type is defined by regex patterns matched against
# Model.csv$OncotreeLineage and Model.csv$OncotreePrimaryDisease
# (case-insensitive). ALL lines whose metadata matches EITHER pattern are
# pulled in -- there is no cap on the number of lines per type.
#
# Keep these identical to the CRISPR script so the two analyses cover the
# same set of cell lines. If a cancer type resolves to zero or an
# unexpectedly small number of lines, check:
#   sort(unique(model$OncotreeLineage))
#   sort(unique(model$OncotreePrimaryDisease))
CANCER_TYPE_DEFINITIONS <- list(
  Breast = list(lineage_pattern = "^Breast$",
                disease_pattern = "Breast"),
  Colon  = list(lineage_pattern = "^Bowel$",
                disease_pattern = "Colorectal|Colon|Cecum|Rectal|Rectum"),
  Liver  = list(lineage_pattern = "^Liver$",
                disease_pattern = "Hepatocellular|Liver|Cholangiocarcinoma|Bile Duct")
)

# -- 1e. Expression scale -------------------------------------------------------
# DepMap's public expression file is already log2(TPM+1).
INPUT_IS_LOG2 <- TRUE

# -- 1f. Colours: one per cancer type (individual lines are no longer a fixed,
#        small set, so colour by cancer-type group rather than by line name) --
cancer_type_colours <- c(Breast = "#E31A1C", Colon = "#1F78B4", Liver = "#33A02C",
                         Other = "#BBBBBB")

# =============================================================================
# SECTION 2 - LOAD DATA
# =============================================================================
cat("--- Loading data ---\n")

model <- read.csv(model_file, stringsAsFactors = FALSE)
cat("Model metadata:", nrow(model), "cell lines\n\n")

load_depmap_expression <- function(path) {
  cat("Loading DepMap AllGenes expression file (this can be large)...\n")
  if (!file.exists(path)) {
    stop(sprintf(
      "Expression file not found at:\n  %s\nDownload 'Omics Expression All Genes TPM Logp1' (OmicsExpressionTPMLogp1HumanAllGenes.csv) from https://depmap.org/portal/download/ and place it in %s",
      path, data_dir
    ))
  }
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  
  # --- Robust ID-column detection --------------------------------------------
  # Do NOT assume column 1 is the ModelID column. DepMap public exports
  # sometimes carry a leading unnamed pandas index column ahead of the real
  # "ModelID" (ACH-xxxxxx) column. Scan the first few columns for the ACH-
  # pattern and use whichever one actually contains it.
  scan_n <- min(5, ncol(raw))
  is_ach_col <- sapply(raw[seq_len(scan_n)], function(col) {
    mean(grepl("^ACH-", as.character(col)), na.rm = TRUE) > 0.5
  })
  if (any(is_ach_col)) {
    id_col_idx <- which(is_ach_col)[1]
  } else {
    id_col_idx <- 1
    cat("  WARNING: no column among the first", scan_n, "looked like ACH- ModelIDs;",
        "falling back to column 1 as the ID column. Verify rownames below look correct.\n")
  }
  cat(sprintf("  Using column %d ('%s') as the cell-line ID column\n",
              id_col_idx, colnames(raw)[id_col_idx]))
  
  ids <- as.character(raw[[id_col_idx]])
  
  # --- Select GOI columns ONLY (avoid building a huge full matrix) -----------
  gene_col_idx   <- setdiff(seq_len(ncol(raw)), id_col_idx)
  raw_gene_names <- colnames(raw)[gene_col_idx]
  clean_names    <- sub("\\s*\\(.*\\)$", "", trimws(raw_gene_names))
  
  match_idx <- which(clean_names %in% GOI)
  if (length(match_idx) == 0) {
    stop("None of the 15 GOI matched any column in the expression file. Check the file format / gene symbol column headers.")
  }
  
  # If a symbol appears more than once (e.g. duplicate Ensembl entries), keep
  # the first occurrence and warn, rather than silently duplicating columns.
  matched_clean <- clean_names[match_idx]
  dup <- duplicated(matched_clean)
  if (any(dup)) {
    cat("  NOTE: duplicate columns found for:", paste(unique(matched_clean[dup]), collapse = ", "),
        "-- keeping first occurrence of each.\n")
    match_idx     <- match_idx[!dup]
    matched_clean <- matched_clean[!dup]
  }
  
  sub_cols <- gene_col_idx[match_idx]
  sub_df   <- raw[, sub_cols, drop = FALSE]
  
  # Explicit numeric coercion, column by column, so any stray non-numeric
  # entry only NAs out that single cell (with a warning) instead of silently
  # corrupting the entire matrix.
  numeric_mat <- matrix(NA_real_, nrow = nrow(sub_df), ncol = ncol(sub_df),
                        dimnames = list(ids, matched_clean))
  for (j in seq_len(ncol(sub_df))) {
    v <- suppressWarnings(as.numeric(sub_df[[j]]))
    n_bad <- sum(is.na(v) & !is.na(sub_df[[j]]) & sub_df[[j]] != "")
    if (n_bad > 0) {
      cat(sprintf("  WARNING: %d non-numeric value(s) coerced to NA in column '%s'\n",
                  n_bad, matched_clean[j]))
    }
    numeric_mat[, j] <- v
  }
  
  cat("Expression matrix (GOI only):", nrow(numeric_mat), "lines x", ncol(numeric_mat), "genes\n\n")
  numeric_mat
}

expr_goi <- load_depmap_expression(expression_file)
goi_found <- colnames(expr_goi)

# =============================================================================
# SECTION 3 - GENE MATCH REPORT + DYNAMIC ALL-LINES-PER-CANCER-TYPE RESOLUTION
# =============================================================================
cat(strrep("=", 60), "\n")
cat(" SECTION 3: Gene match report + cancer-type line resolution\n")
cat(strrep("=", 60), "\n\n")

missing_goi <- setdiff(GOI, goi_found)
if (length(missing_goi) > 0) {
  cat("  [Expression] NOT found (", length(missing_goi), "): ", paste(missing_goi, collapse = ", "), "\n", sep = "")
}
cat(sprintf("  [Expression] matched %d / %d GOI\n", length(goi_found), length(GOI)))
cat("GOI present in expression matrix:", length(goi_found), "/", length(GOI), "\n\n")

expr_rn <- rownames(expr_goi)
cat("Example ModelIDs (first 5):", paste(head(expr_rn, 5), collapse = ", "), "\n")
if (!any(grepl("^ACH-", head(expr_rn, 20)))) {
  cat("  WARNING: ModelIDs do not look like ACH-xxxxxx. Cell-line matching below may fail --\n")
  cat("  double check the ID-column detection message above.\n")
}
cat("\n")

# -- 3a. Dynamically resolve ALL cell lines per cancer type --------------------
# expr_goi rownames and model$ModelID are both ACH-xxxxxx identifiers, so no
# fuzzy name matching is required here -- we match directly on ModelID after
# filtering Model.csv by Oncotree lineage / primary disease.
resolve_cancer_type_lines <- function(type_name, def, model_df, profiled_ids) {
  lineage_hit <- grepl(def$lineage_pattern, model_df$OncotreeLineage,
                       ignore.case = TRUE)
  disease_hit <- grepl(def$disease_pattern, model_df$OncotreePrimaryDisease,
                       ignore.case = TRUE)
  matched <- model_df[lineage_hit | disease_hit, , drop = FALSE]
  
  # Keep only lines that were actually expression-profiled (present in expr_goi)
  matched <- matched[matched$ModelID %in% profiled_ids, , drop = FALSE]
  matched <- matched[!duplicated(matched$ModelID), , drop = FALSE]
  
  cat(sprintf("  %-8s: %d Model.csv matches; %d with expression data\n",
              type_name,
              sum(lineage_hit | disease_hit),
              nrow(matched)))
  
  matched %>%
    transmute(ModelID, CellLineName, CancerType = type_name)
}

cat("--- Resolving ALL cell lines per cancer type (expression-profiled only) ---\n")
focal_meta <- map2_dfr(names(CANCER_TYPE_DEFINITIONS), CANCER_TYPE_DEFINITIONS,
                       resolve_cancer_type_lines,
                       model_df = model, profiled_ids = expr_rn)
cat("\n")

if (nrow(focal_meta) == 0)
  stop("No cell lines resolved for any cancer type -- check the regex patterns ",
       "in CANCER_TYPE_DEFINITIONS against sort(unique(model$OncotreeLineage)) ",
       "and sort(unique(model$OncotreePrimaryDisease)).")

cat("--- Lines resolved per cancer type ---\n")
print(focal_meta %>% count(CancerType, name = "n_lines"))
cat("\n")

focal_rows <- focal_meta$ModelID
focal_rows <- focal_rows[focal_rows %in% rownames(expr_goi)]

cat("--- Focal line -> cancer type map (first 10 shown) ---\n")
print(head(focal_meta, 10))
cat("...\n\n")

# -- 3b. Optional: named-line lookup helper (kept for spot-checking) ----------
possible_id_cols <- grep("id|alias", colnames(model), ignore.case = TRUE, value = TRUE)
known_id_cols <- c("ModelID", "CellLineName", "StrippedCellLineName",
                   "SangerModelID", "OncotreeCode", "CCLEName", "ModelIDAlias")
id_cols <- unique(c(known_id_cols, possible_id_cols))
id_cols <- intersect(id_cols, colnames(model))

find_line_by_name <- function(name, model_df = model, id_columns = id_cols, rn = expr_rn) {
  clean_name <- gsub("[-_ ./]", "", name)
  for (col in id_columns) {
    idx <- grep(name, model_df[[col]], ignore.case = TRUE, fixed = TRUE)
    if (length(idx) > 0) {
      row_vals <- as.character(model_df[idx[1], ])
      hit <- row_vals[row_vals %in% rn]
      if (length(hit) > 0) return(hit[1])
    }
  }
  NULL
}

# =============================================================================
# SECTION 4 - EXTRACT EXPRESSION FOR FOCAL VS OTHER LINES
# =============================================================================
cat(strrep("=", 60), "\n")
cat(" SECTION 4: Extract expression values\n")
cat(strrep("=", 60), "\n\n")

expr_focal  <- expr_goi[focal_rows, goi_found, drop = FALSE]
expr_others <- expr_goi[setdiff(rownames(expr_goi), focal_rows), goi_found, drop = FALSE]

cat("Focal matrix (all breast+liver+colon lines):", nrow(expr_focal),  "lines x", ncol(expr_focal),  "genes\n")
cat("Other matrix (all remaining profiled lines):", nrow(expr_others), "lines x", ncol(expr_others), "genes\n\n")

to_log2 <- function(m) if (INPUT_IS_LOG2) m else log2(m + 1)
to_tpm  <- function(m) if (INPUT_IS_LOG2) (2^m) - 1 else m

expr_focal_log2  <- to_log2(expr_focal)
expr_others_log2 <- to_log2(expr_others)

long_focal <- expr_focal_log2 %>%
  as.data.frame() %>%
  rownames_to_column("ModelID") %>%
  pivot_longer(-ModelID, names_to = "Gene", values_to = "log2TPM") %>%
  left_join(focal_meta, by = "ModelID") %>%
  mutate(Group = "Focal")

long_others <- expr_others_log2 %>%
  as.data.frame() %>%
  rownames_to_column("ModelID") %>%
  pivot_longer(-ModelID, names_to = "Gene", values_to = "log2TPM") %>%
  mutate(Group = "Other")

long_all <- bind_rows(long_focal, long_others %>% mutate(CellLineName = NA_character_, CancerType = NA_character_))

# =============================================================================
# SECTION 5 - SUMMARY STATISTICS + WILCOXON (each line, and each cancer type, vs all others)
# =============================================================================
cat("--- Computing summary statistics ---\n")

other_stats <- long_others %>%
  group_by(Gene) %>%
  summarise(
    n_other      = n(),
    mean_other   = mean(log2TPM, na.rm = TRUE),
    median_other = median(log2TPM, na.rm = TRUE),
    sd_other     = sd(log2TPM, na.rm = TRUE),
    .groups = "drop"
  )

# Per-line summary (every individual line within every cancer type -- no cap)
summary_full <- map_dfr(unique(focal_meta$ModelID), function(model_id) {
  meta_row <- focal_meta %>% filter(ModelID == model_id)
  line_name <- meta_row$CellLineName[1]
  
  focal_vals <- long_focal %>%
    filter(ModelID == model_id) %>%
    select(Gene, CancerType, log2TPM_focal = log2TPM)
  
  wilcox_res <- map_dfr(goi_found, function(gene) {
    fc <- focal_vals %>% filter(Gene == gene) %>% pull(log2TPM_focal)
    oc <- long_others %>% filter(Gene == gene) %>% pull(log2TPM)
    if (length(fc) == 0 || all(is.na(oc))) {
      return(tibble(Gene = gene, W_statistic = NA_real_, p_value = NA_real_))
    }
    wt <- tryCatch(wilcox.test(fc, oc), error = function(e) list(statistic = NA, p.value = NA))
    tibble(Gene = gene, W_statistic = as.numeric(wt$statistic), p_value = wt$p.value)
  })
  
  focal_vals %>%
    left_join(other_stats, by = "Gene") %>%
    left_join(wilcox_res, by = "Gene") %>%
    mutate(
      ModelID       = model_id,
      CellLineName  = line_name,
      Delta_log2TPM = log2TPM_focal - mean_other,
      Z_score       = (log2TPM_focal - mean_other) / sd_other,
      TPM_focal     = (2^log2TPM_focal) - 1
    )
}) %>%
  group_by(CellLineName) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(CancerType, CellLineName, desc(Z_score)) %>%
  select(CancerType, CellLineName, ModelID, Gene, TPM_focal, log2TPM_focal, mean_other, median_other,
         sd_other, Delta_log2TPM, Z_score, W_statistic, p_value, p_adj_BH)

cat("  Summary table rows:", nrow(summary_full), "\n\n")

# =============================================================================
# SECTION 5b - WITHIN-CANCER-TYPE CONCORDANCE (across ALL lines of that type)
# =============================================================================
# With every available line per cancer type (n can now be well above 2), this
# asks: do the lines within a cancer type agree on which genes are elevated?
# This directly extends the previous "n=1/n=2 per cancer type lacks power"
# check to the full available sample size for each type.
cat("--- Computing within-cancer-type concordance ---\n")

concordance <- summary_full %>%
  group_by(CancerType, Gene) %>%
  summarise(
    n_lines_in_type  = n(),
    mean_log2TPM     = mean(log2TPM_focal, na.rm = TRUE),
    sd_log2TPM       = sd(log2TPM_focal, na.rm = TRUE),
    mean_Z           = mean(Z_score, na.rm = TRUE),
    frac_lines_Z_gt1 = mean(Z_score > 1, na.rm = TRUE),   # fraction of lines elevated (Z>1)
    all_lines_Z_gt1  = all(Z_score > 1),                  # every line in the type elevated
    frac_same_sign   = mean(sign(Z_score) == sign(mean_Z), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(CancerType, desc(mean_Z))

cat("  Concordance table rows:", nrow(concordance), "\n\n")

# =============================================================================
# SECTION 6 - SAVE EXCEL OUTPUTS
# =============================================================================
cat("--- Saving Excel outputs ---\n")

expr_wide_focal <- expr_focal_log2 %>%
  as.data.frame() %>%
  rownames_to_column("ModelID") %>%
  left_join(focal_meta, by = "ModelID")

write.xlsx(
  list(
    "Focal_log2TPM"  = expr_wide_focal,
    "Focal_Long"     = long_focal %>% select(CancerType, CellLineName, ModelID, Gene, log2TPM)
  ),
  file = file.path(results_dir, "01_Expression_focal_vs_others.xlsx"), rowNames = FALSE
)

write.xlsx(
  list("Summary" = summary_full),
  file = file.path(results_dir, "02_Summary_statistics.xlsx"), rowNames = FALSE
)

zscore_wide <- summary_full %>%
  select(CellLineName, Gene, Z_score) %>%
  pivot_wider(names_from = CellLineName, values_from = Z_score)

write.xlsx(
  list("ZScore" = zscore_wide),
  file = file.path(results_dir, "03_ZScore_table.xlsx"), rowNames = FALSE
)

write.xlsx(
  list("CancerType_Concordance" = concordance),
  file = file.path(results_dir, "04_CancerType_concordance.xlsx"), rowNames = FALSE
)

write.xlsx(
  list("CancerType_line_manifest" = focal_meta %>% arrange(CancerType, CellLineName)),
  file = file.path(results_dir, "05_CancerType_line_manifest.xlsx"), rowNames = FALSE
)

cat("  Saved 5 Excel workbooks to", results_dir, "\n\n")

# =============================================================================
# SECTION 7 - FIGURES
# =============================================================================
gene_order <- summary_full %>%
  group_by(Gene) %>%
  summarise(mean_focal = mean(log2TPM_focal, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_focal)) %>%
  pull(Gene)

# -- 7a. Violin plots: one per gene, ALL focal lines shown as points, coloured
#        by cancer type (individual-line colour coding no longer used since
#        the number of lines per type is not fixed at 2) ---------------------
cat("--- Generating violin plots ---\n")
for (gene in goi_found) {
  plot_dat <- long_all %>%
    filter(Gene == gene) %>%
    mutate(PlotGroup = if_else(Group == "Focal", CancerType, "Other"),
           PlotGroup = factor(PlotGroup, levels = c(names(CANCER_TYPE_DEFINITIONS), "Other")))
  
  focal_points <- plot_dat %>% filter(PlotGroup != "Other")
  
  p <- ggplot(plot_dat, aes(x = PlotGroup, y = log2TPM)) +
    geom_violin(data = plot_dat %>% filter(PlotGroup == "Other"),
                scale = "width", alpha = 0.45, fill = "grey70", colour = "grey60") +
    geom_boxplot(data = plot_dat %>% filter(PlotGroup == "Other"),
                 width = 0.15, outlier.shape = NA, fill = "white", colour = "grey40") +
    geom_jitter(data = focal_points, aes(colour = PlotGroup), width = 0.15, size = 2, alpha = 0.8) +
    scale_colour_manual(values = cancer_type_colours) +
    labs(title = gene,
         subtitle = sprintf(
           "log2(TPM+1) | all Breast/Liver/Colon lines (n=%d) vs all others (n=%d)",
           nrow(focal_meta), nrow(expr_others)),
         x = NULL, y = "log2(TPM + 1)", colour = "Cancer type") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave(file.path(figures_dir, paste0("Violin_", gene, ".pdf")), p, width = 8, height = 5, useDingbats = FALSE)
}
cat("  Saved", length(goi_found), "violin PDFs\n")

# -- 7b. Heatmap: GOI x focal lines (log2TPM), grouped by cancer type --------
cat("--- Generating heatmap ---\n")
heatmap_meta <- focal_meta %>% arrange(CancerType, CellLineName)
heatmap_mat  <- t(expr_focal_log2[heatmap_meta$ModelID, gene_order, drop = FALSE])
colnames(heatmap_mat) <- heatmap_meta$CellLineName

ann_col <- data.frame(CancerType = heatmap_meta$CancerType, row.names = heatmap_meta$CellLineName)
ann_colours <- list(CancerType = cancer_type_colours[c("Breast", "Colon", "Liver")])

pdf(file.path(figures_dir, "Heatmap_Expression_GOI.pdf"),
    width = max(7.5, nrow(heatmap_meta) * 0.3 + 3),
    height = max(5, length(gene_order) * 0.4 + 2), useDingbats = FALSE)
pheatmap(
  heatmap_mat,
  color = colorRampPalette(brewer.pal(9, "OrRd"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_col = ann_col,
  annotation_colors = ann_colours,
  display_numbers = TRUE, number_format = "%.2f",
  fontsize_number = 6, fontsize_row = 10, fontsize_col = 7,
  main = sprintf("log2(TPM+1) expression - GOI x %d focal cell lines (Breast/Colon/Liver)",
                 nrow(heatmap_meta))
)
invisible(dev.off())
cat("  Saved: Heatmap_Expression_GOI.pdf\n")

# -- 7c. Z-score selectivity dot plot, coloured by cancer type ---------------
cat("--- Generating Z-score dot plot ---\n")
summary_ordered <- summary_full %>%
  mutate(Gene = factor(Gene, levels = rev(gene_order)),
         CancerType = factor(CancerType, levels = names(CANCER_TYPE_DEFINITIONS)))

p_z <- ggplot(summary_ordered, aes(x = Gene, y = Z_score, colour = CancerType,
                                   size = -log10(p_value + 1e-10))) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey60") +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", linewidth = 0.4, colour = "orange") +
  geom_point(alpha = 0.75, position = position_jitter(width = 0.25, height = 0)) +
  scale_colour_manual(values = cancer_type_colours) +
  scale_size_continuous(range = c(1.5, 6), name = "-log10(p-value)") +
  coord_flip() +
  facet_wrap(~CancerType, ncol = 1, scales = "free_y") +
  labs(title = "Focal-line expression Z-score vs all other profiled cell lines",
       subtitle = "Positive Z = higher expression in focal line than average; every line per cancer type shown; faceted by cancer type",
       x = NULL, y = "Z-score (focal log2TPM vs mean/SD of all others)",
       colour = "Cancer type", caption = "Orange dashed = |Z| > 2") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")

ggsave(file.path(figures_dir, "DotPlot_Zscore_GOI.pdf"), p_z, width = 9, height = 10, useDingbats = FALSE)
cat("  Saved: DotPlot_Zscore_GOI.pdf\n\n")

# =============================================================================
# SECTION 8 - CONSOLE SUMMARY
# =============================================================================
cat(strrep("=", 60), "\n")
cat(" FINAL SUMMARY\n")
cat(strrep("=", 60), "\n\n")

cat("--- GOI found in expression matrix ---\n")
cat(" ", paste(goi_found, collapse = ", "), "\n\n")
if (length(setdiff(GOI, goi_found)) > 0) {
  cat("--- GOI NOT found (check these did not silently drop from AllGenes file) ---\n")
  cat(" ", paste(setdiff(GOI, goi_found), collapse = ", "), "\n\n")
}

cat("--- Number of expression-profiled lines included per cancer type ---\n")
print(focal_meta %>% count(CancerType, name = "n_lines"))
cat("\n")

cat("--- Focal-line log2(TPM+1) expression (GOI), first 15 rows shown ---\n")
print(expr_wide_focal %>% select(CancerType, CellLineName, all_of(gene_order)) %>%
        arrange(CancerType, CellLineName) %>% head(15))

cat("\n--- Genes with |Z| > 2 in at least one focal line ---\n")
hi <- summary_full %>% filter(abs(Z_score) > 2)
if (nrow(hi) > 0) {
  print(hi %>% select(CancerType, CellLineName, Gene, TPM_focal, Z_score, p_value, p_adj_BH))
} else {
  cat("  None at this threshold.\n")
}

cat("\n--- Genes where ALL lines within a cancer type are elevated (Z > 1) ---\n")
conc_hits <- concordance %>% filter(all_lines_Z_gt1)
if (nrow(conc_hits) > 0) {
  print(conc_hits %>% select(CancerType, Gene, n_lines_in_type, mean_log2TPM, mean_Z))
} else {
  cat("  None at this threshold.\n")
}

cat("\n--- Output files ---\n")
for (f in list.files(results_dir, recursive = TRUE, full.names = TRUE)) {
  cat(sprintf("  %6.1f KB  %s\n", file.size(f) / 1e3, sub(results_dir, "R2_CCLE_results", f)))
}

cat("\nAnalysis complete.\n")
