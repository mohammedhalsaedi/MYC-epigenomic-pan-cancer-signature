# ==============================================================================
# Script: 02_publication_multi_cancer_analysis.R
# Purpose: Downstream cross-cancer characterisation of the 15-gene MYC signature,
#          including Kaplan-Meier/Cox survival modeling, Pearson co-expression,
#          tumour-vs-normal ROC discrimination, and permutation benchmarking.
# Inputs:  Aligned TCGA expression matrices and clinical follow-up files.
# Outputs: Survival statistics, correlation matrices, ROC curves, and validation summaries.
# ==============================================================================

# --- SETUP ---
pkgs <- c("data.table","survival","survminer","ggplot2","pheatmap",
          "RColorBrewer","pROC","gridExtra","cowplot","ggpubr",
          "dplyr","tidyr","org.Hs.eg.db","httr","jsonlite")

miss <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(miss)) {
  if(!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org")
  BiocManager::install(miss, ask=FALSE, update=FALSE)
}
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only=TRUE)))

# Pathing Configuration
base_dir <- getwd()
data_dir <- file.path(base_dir, "data")
clinical_dir <- file.path(base_dir, "clinical")
out_dir <- file.path(base_dir, "results", "analysis")
fig_dir <- file.path(out_dir, "figures")
surv_fig_dir <- file.path(fig_dir, "survival")
table_dir <- file.path(out_dir, "tables")

dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(surv_fig_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(table_dir, showWarnings=FALSE, recursive=TRUE)

# Signature Genes
common_genes <- c("MLLT10", "BYSL", "MYC", "SLC7A1", "SNHG16", 
                  "FOXK1", "CHD7", "PRKAR1B", "SNHG17", "FCHSD1", 
                  "ZBTB10", "SNORA71B", "CASC11", "PVT1", "BOP1")
cancers <- c("BRCA", "LIHC", "COAD")

cat("\nPUBLICATION ANALYSIS: 15 COMMON MYC-EPIGENETIC GENES\n")
cat("Cancers:", paste(cancers, collapse=", "), "\n")
cat("Genes:", length(common_genes), "\n")
cat("Output:", out_dir, "\n\n")

# --- HELPER FUNCTIONS ---
gdc_map <- function(uuids, batch=200) {
  uuids <- unique(uuids)
  batches <- split(uuids, ceiling(seq_along(uuids)/batch))
  res <- list()
  
  for(i in seq_along(batches)) {
    filt <- jsonlite::toJSON(
      list(op="in", content=list(field="file_id", value=as.list(batches[[i]]))),
      auto_unbox=TRUE)
    
    resp <- tryCatch(
      httr::POST("https://api.gdc.cancer.gov/files",
                 body=list(filters=filt, fields="file_id,associated_entities.entity_submitter_id",
                           format="JSON", size=as.character(batch)),
                 encode="form", httr::timeout(60)),
      error=function(e) NULL)
    
    if(!is.null(resp) && httr::status_code(resp)==200) {
      parsed <- jsonlite::fromJSON(httr::content(resp, as="text", encoding="UTF-8"), flatten=TRUE)
      if(!is.null(parsed$data$hits)) res[[i]] <- parsed$data$hits
    }
    Sys.sleep(0.5)
  }
  
  comb <- do.call(rbind, res[!sapply(res, is.null)])
  if(is.null(comb) || nrow(comb)==0) stop("API failed")
  
  df <- data.frame(
    uuid = comb$file_id,
    barcode = sapply(comb$associated_entities, function(e) {
      if(is.data.frame(e) && nrow(e)>0) e$entity_submitter_id[1]
      else if(is.list(e) && length(e)>0) e[[1]]$entity_submitter_id
      else NA_character_
    }),
    stringsAsFactors=FALSE)
  df[!is.na(df$barcode),]
}

load_expr <- function(dir) {
  cat("  Loading expression from", basename(dir), "...\n")
  folders <- list.dirs(dir, recursive=FALSE, full.names=TRUE)
  files <- sapply(folders, function(f) {
    x <- list.files(f, pattern="\\.tsv(\\.gz)?$", full.names=TRUE)
    if(length(x)) x[1] else NA_character_
  })
  
  valid <- !is.na(files)
  files <- files[valid]
  uuids <- basename(folders)[valid]
  
  peek <- data.table::fread(files[1], nrows=10)
  
  gene_col <- NULL
  if("gene_name" %in% names(peek)) {
    gene_col <- "gene_name"
  } else if("gene_id" %in% names(peek)) {
    sample_id <- peek$gene_id[1]
    if(grepl("^ENSG", sample_id)) {
      if(ncol(peek) > 1) gene_col <- names(peek)[2]
    } else {
      gene_col <- "gene_id"
    }
  } else {
    gene_col <- names(peek)[1]
  }
  
  val_col <- grep("unstranded|tpm_unstranded|fpkm", names(peek), ignore.case=TRUE, value=TRUE)[1]
  if(is.na(val_col)) {
    numeric_cols <- names(peek)[sapply(peek, is.numeric)]
    val_col <- tail(numeric_cols, 1)
  }
  
  first <- data.table::fread(files[1])
  genes_raw <- first[[gene_col]]
  
  valid_genes <- !is.na(genes_raw) & nzchar(trimws(genes_raw))
  genes_raw <- genes_raw[valid_genes]
  values_raw <- first[[val_col]][valid_genes]
  
  if(all(grepl("^ENSG", head(genes_raw, 100)))) {
    gene_map <- AnnotationDbi::select(org.Hs.eg.db, 
                                      keys=genes_raw,
                                      columns="SYMBOL",
                                      keytype="ENSEMBL")
    
    gene_map <- gene_map[!is.na(gene_map$SYMBOL) & !duplicated(gene_map$ENSEMBL), ]
    
    keep_idx <- match(gene_map$ENSEMBL, genes_raw)
    genes_raw <- genes_raw[keep_idx]
    genes <- gene_map$SYMBOL
    values_raw <- values_raw[keep_idx]
  } else {
    genes <- genes_raw
  }
  
  dup_genes <- duplicated(genes)
  if(any(dup_genes)) {
    keep <- !dup_genes
    genes <- genes[keep]
    genes_raw <- genes_raw[keep]
    values_raw <- values_raw[keep]
  }
  
  mat <- matrix(NA_real_, nrow=length(genes), ncol=length(files), 
                dimnames=list(genes, uuids))
  
  cat("  Reading")
  for(i in seq_along(files)) {
    dt <- data.table::fread(files[i])
    gene_ids <- dt[[gene_col]]
    values <- dt[[val_col]]
    
    valid_idx <- !is.na(gene_ids) & nzchar(trimws(gene_ids))
    gene_ids <- gene_ids[valid_idx]
    values <- values[valid_idx]
    
    if(exists("gene_map")) {
      idx <- match(gene_map$ENSEMBL, gene_ids)
      mat[, i] <- values[idx]
    } else {
      idx <- match(genes, gene_ids)
      mat[, i] <- values[idx]
    }
  }
  cat(" Done\n")
  
  list(mat=mat, uuids=uuids, genes=rownames(mat))
}

find_column <- function(df, possible_names) {
  for(name in possible_names) {
    if(name %in% names(df)) return(name)
    
    prefixed <- paste0(c("cases.", "follow_ups.", "diagnoses.", "demographics.", "demographic."), name)
    for(p in prefixed) {
      if(p %in% names(df)) return(p)
    }
  }
  return(NULL)
}

load_clinical <- function(cancer) {
  cat("  Loading clinical data for", cancer, "...\n")
  
  clinical_folder <- file.path(clinical_dir, paste0("clinical.cart.", cancer, "_Expression"))
  
  if(!dir.exists(clinical_folder)) stop("Clinical folder not found: ", clinical_folder)
  
  clinical_file <- file.path(clinical_folder, "clinical.tsv")
  if(file.exists(clinical_file)) {
    clinical <- data.table::fread(clinical_file, data.table=FALSE)
  } else {
    stop("clinical.tsv not found in ", clinical_folder)
  }
  
  followup_file <- file.path(clinical_folder, "follow_up.tsv")
  if(file.exists(followup_file)) {
    followup <- data.table::fread(followup_file, data.table=FALSE)
  } else {
    followup <- NULL
  }
  
  cat("    Clinical records:", nrow(clinical), "\n")
  if(!is.null(followup)) {
    cat("    Follow-up records:", nrow(followup), "\n")
  }
  
  list(clinical = clinical, followup = followup)
}

# --- PART 1: DATA LOADING ---
cat("\n--- PART 1: DATA LOADING ---\n")
all_data <- list()
for(cancer in cancers) {
  cat("\n", cancer, ":\n", sep="")
  
  expr_dir <- file.path(data_dir, paste0(cancer, "_Expression"))
  if(!dir.exists(expr_dir)) {
    cat("  ERROR: Directory not found:", expr_dir, "\n")
    next
  }
  
  E <- load_expr(expr_dir)
  cat("  Mapping UUIDs to barcodes...\n")
  expr_map <- gdc_map(E$uuids)
  
  for(i in 1:nrow(expr_map)) {
    j <- which(colnames(E$mat) == expr_map$uuid[i])
    if(length(j)) colnames(E$mat)[j] <- expr_map$barcode[i]
  }
  
  E$mat <- E$mat[, colnames(E$mat) %in% expr_map$barcode, drop=FALSE]
  cat("  Expression matrix:", nrow(E$mat), "genes ×", ncol(E$mat), "samples\n")
  
  clinical_data <- load_clinical(cancer)
  clinical <- clinical_data$clinical
  followup <- clinical_data$followup
  
  barcodes <- colnames(E$mat)
  patient_ids <- substr(barcodes, 1, 12)
  
  survival_data <- data.frame(
    barcode = barcodes,
    patient_id = patient_ids,
    stringsAsFactors = FALSE
  )
  
  sample_code <- substr(barcodes, 14, 15)
  survival_data$sample_type <- ifelse(
    sample_code == "01", "Primary Tumor",
    ifelse(sample_code == "11", "Solid Tissue Normal",
           ifelse(sample_code %in% c("02", "03", "04", "05", "06", "07"), "Metastatic/Recurrent",
                  "Other"))
  )
  
  survival_data$OS_time <- NA
  survival_data$OS_status <- 0
  
  cat("  Extracting survival data...\n")
  patient_col <- find_column(clinical, c("submitter_id", "case_submitter_id", "case_id"))
  
  if(!is.null(patient_col)) {
    clinical_match <- match(patient_ids, clinical[[patient_col]])
    
    vital_col <- find_column(clinical, c("vital_status"))
    if(!is.null(vital_col)) {
      vital_status <- clinical[[vital_col]][clinical_match]
      survival_data$OS_status <- ifelse(
        !is.na(vital_status) & tolower(as.character(vital_status)) %in% c("dead", "deceased"), 1, 0
      )
    }
    
    death_col <- find_column(clinical, c("days_to_death"))
    if(!is.null(death_col)) {
      death_days <- clinical[[death_col]][clinical_match]
      death_days <- ifelse(death_days == "'--", NA, death_days)
      survival_data$OS_time <- as.numeric(death_days)
    }
    
    fu_col <- find_column(clinical, c("days_to_last_follow_up", "days_to_last_followup"))
    if(!is.null(fu_col)) {
      fu_days <- clinical[[fu_col]][clinical_match]
      fu_days <- ifelse(fu_days == "'--", NA, fu_days)
      fu_days <- as.numeric(fu_days)
      
      survival_data$OS_time <- ifelse(
        is.na(survival_data$OS_time) | survival_data$OS_status == 0,
        fu_days,
        survival_data$OS_time
      )
    }
  }
  
  if(!is.null(followup) && nrow(followup) > 0) {
    followup_patient_col <- find_column(followup, c("submitter_id", "case_submitter_id"))
    followup_time_col <- find_column(followup, c("days_to_follow_up"))
    
    if(!is.null(followup_patient_col) && !is.null(followup_time_col)) {
      for(i in 1:length(patient_ids)) {
        pt_id <- patient_ids[i]
        pt_followup <- followup[followup[[followup_patient_col]] == pt_id, ]
        
        if(nrow(pt_followup) > 0) {
          fu_times <- as.numeric(pt_followup[[followup_time_col]])
          fu_times <- fu_times[!is.na(fu_times) & fu_times != "'--"]
          
          if(length(fu_times) > 0) {
            max_fu <- max(fu_times)
            
            if(is.na(survival_data$OS_time[i]) || 
               (survival_data$OS_status[i] == 0 && max_fu > survival_data$OS_time[i])) {
              survival_data$OS_time[i] <- max_fu
            }
          }
        }
      }
    }
  }
  
  survival_data$OS_time <- as.numeric(survival_data$OS_time)
  survival_data$OS_status <- as.numeric(survival_data$OS_status)
  survival_data$OS_time[survival_data$OS_time <= 0] <- NA
  
  n_surv <- sum(!is.na(survival_data$OS_time) & survival_data$OS_time > 0)
  n_dead <- sum(survival_data$OS_status == 1, na.rm=TRUE)
  
  cat("  Survival data:", n_surv, "patients (", n_dead, "events)\n", sep="")
  
  all_data[[cancer]] <- list(
    expr = E$mat,
    survival = survival_data,
    clinical = clinical
  )
}
cat("\nData loaded successfully!\n\n")

# --- PART 2: SURVIVAL ANALYSIS ---
cat("\n--- PART 2: SURVIVAL ANALYSIS ---\n")
perform_survival_analysis <- function(cancer, gene, data, save_plot=TRUE) {
  expr <- data$expr
  surv_data <- data$survival
  
  if(!gene %in% rownames(expr)) return(NULL)
  gene_expr <- expr[gene, ]
  
  tumor_samples <- surv_data$sample_type == "Primary Tumor"
  valid_surv <- !is.na(surv_data$OS_time) & 
    !is.na(surv_data$OS_status) & 
    surv_data$OS_time > 0
  
  use_samples <- tumor_samples & valid_surv
  if(sum(use_samples) < 20) return(NULL)
  
  surv_df <- data.frame(
    time = surv_data$OS_time[use_samples],
    status = surv_data$OS_status[use_samples],
    expression = gene_expr[use_samples]
  )
  
  surv_df <- surv_df[surv_df$time > 0, ]
  if(nrow(surv_df) < 20) return(NULL)
  
  surv_df$group <- ifelse(surv_df$expression >= median(surv_df$expression), "High", "Low")
  surv_df$group <- factor(surv_df$group, levels=c("Low", "High"))
  
  fit <- survfit(Surv(time, status) ~ group, data=surv_df)
  logrank <- survdiff(Surv(time, status) ~ group, data=surv_df)
  p_value <- 1 - pchisq(logrank$chisq, 1)
  
  cox <- coxph(Surv(time, status) ~ expression, data=surv_df)
  hr <- exp(coef(cox))
  hr_ci <- exp(confint(cox))
  cox_p <- summary(cox)$coefficients[5]
  
  if(save_plot) {
    pdf(file.path(surv_fig_dir, paste0("survival_", cancer, "_", gene, ".pdf")), 
        width=8, height=6)
    
    plot(fit, col=c("blue", "red"), lwd=3,
         xlab="Time (days)", ylab="Overall Survival Probability",
         main=sprintf("%s - %s\nHR=%.2f (95%% CI: %.2f-%.2f), P=%.4f", 
                      cancer, gene, hr, hr_ci[1], hr_ci[2], p_value),
         cex.main=1.2)
    
    legend("topright", legend=c(sprintf("Low (n=%d)", sum(surv_df$group=="Low")),
                                sprintf("High (n=%d)", sum(surv_df$group=="High"))),
           col=c("blue", "red"), lwd=3, title=paste(gene, "Expression"), bty="n")
    
    text(max(surv_df$time)*0.02, 0.05, 
         sprintf("Log-rank P = %.4f\nCox P = %.4f", p_value, cox_p), adj=c(0, 0), cex=0.9)
    dev.off()
  }
  
  list(fit = fit, surv_df = surv_df, p_value = p_value, cox_hr = hr, cox_ci = hr_ci, cox_p = cox_p, n = nrow(surv_df))
}

survival_results <- list()
survival_table <- data.frame()

for(cancer in cancers) {
  if(is.null(all_data[[cancer]])) next
  for(gene in common_genes) {
    result <- perform_survival_analysis(cancer, gene, all_data[[cancer]], save_plot=TRUE)
    if(!is.null(result)) {
      survival_results[[paste(cancer, gene, sep="_")]] <- result
      survival_table <- rbind(survival_table, data.frame(
        Cancer = cancer, Gene = gene, N = result$n, LogRank_P = result$p_value,
        HR = result$cox_hr, HR_lower = result$cox_ci[1], HR_upper = result$cox_ci[2],
        Cox_P = result$cox_p, Significant = ifelse(result$p_value < 0.05, "Yes", "No")
      ))
    }
  }
}
write.csv(survival_table, file.path(table_dir, "survival_analysis_all.csv"), row.names=FALSE)
cat("Survival analysis completed\n\n")

# --- PART 3: GENE-GENE CORRELATION ---
cat("\n--- PART 3: GENE-GENE CORRELATION ---\n")
correlation_matrices <- list()
for(cancer in cancers) {
  if(is.null(all_data[[cancer]])) next
  expr <- all_data[[cancer]]$expr
  surv_data <- all_data[[cancer]]$survival
  
  tumor_idx <- surv_data$sample_type == "Primary Tumor"
  available_genes <- intersect(common_genes, rownames(expr))
  if(length(available_genes) < 2) next
  
  expr_subset <- expr[available_genes, tumor_idx]
  cor_mat <- cor(t(expr_subset), use="pairwise.complete.obs", method="pearson")
  correlation_matrices[[cancer]] <- cor_mat
  
  write.csv(cor_mat, file.path(table_dir, paste0("correlation_matrix_", cancer, ".csv")))
  
  pdf(file.path(fig_dir, paste0("correlation_heatmap_", cancer, ".pdf")), width=10, height=9)
  pheatmap(cor_mat, color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           breaks = seq(-1, 1, length.out=101), display_numbers = TRUE, number_format = "%.2f",
           fontsize_number = 8, cluster_rows = TRUE, cluster_cols = TRUE,
           main = paste(cancer, "- Gene-Gene Correlation"), fontsize = 10)
  dev.off()
}
cat("Correlation matrices calculated\n\n")

# --- PART 4: EXPRESSION HEATMAP ---
cat("\n--- PART 4: EXPRESSION HEATMAP ---\n")
combined_expr <- NULL
cancer_labels <- NULL
for(cancer in cancers) {
  if(is.null(all_data[[cancer]])) next
  expr <- all_data[[cancer]]$expr
  surv_data <- all_data[[cancer]]$survival
  tumor_idx <- surv_data$sample_type == "Primary Tumor"
  
  available_genes <- intersect(common_genes, rownames(expr))
  expr_subset <- expr[available_genes, tumor_idx]
  expr_z <- t(scale(t(expr_subset)))
  
  n_sample <- min(50, ncol(expr_z))
  sample_idx <- sample(1:ncol(expr_z), n_sample)
  
  if(is.null(combined_expr)) combined_expr <- expr_z[, sample_idx]
  else combined_expr <- cbind(combined_expr, expr_z[, sample_idx])
  
  cancer_labels <- c(cancer_labels, rep(cancer, n_sample))
}

annotation_col <- data.frame(Cancer = cancer_labels, row.names = colnames(combined_expr))
annotation_colors <- list(Cancer = c(BRCA = "#E41A1C", LIHC = "#377EB8", COAD = "#4DAF4A"))

pdf(file.path(fig_dir, "expression_heatmap_all_cancers.pdf"), width=14, height=10)
pheatmap(combined_expr, color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-3, 3, length.out=101), annotation_col = annotation_col,
         annotation_colors = annotation_colors, cluster_rows = TRUE, cluster_cols = TRUE,
         show_colnames = FALSE, fontsize_row = 10, main = "Expression of 15 Common Genes Across Cancers (Z-score)",
         border_color = NA)
dev.off()
cat("Expression heatmap created\n\n")

# --- PART 5: ROC CURVES ---
cat("\n--- PART 5: ROC CURVES ---\n")
roc_results <- data.frame()
for(cancer in cancers) {
  if(is.null(all_data[[cancer]])) next
  expr <- all_data[[cancer]]$expr
  surv_data <- all_data[[cancer]]$survival
  
  tumor_idx <- surv_data$sample_type == "Primary Tumor"
  normal_idx <- surv_data$sample_type == "Solid Tissue Normal"
  
  if(sum(normal_idx) < 5) next
  
  for(gene in common_genes) {
    if(!gene %in% rownames(expr)) next
    gene_expr <- expr[gene, ]
    labels <- ifelse(tumor_idx, 1, ifelse(normal_idx, 0, NA))
    valid <- !is.na(labels)
    
    if(sum(valid) < 10) next
    roc_obj <- roc(labels[valid], gene_expr[valid], quiet=TRUE)
    
    roc_results <- rbind(roc_results, data.frame(
      Cancer = cancer, Gene = gene, AUC = as.numeric(auc(roc_obj)),
      AUC_lower = as.numeric(ci.auc(roc_obj))[1], AUC_upper = as.numeric(ci.auc(roc_obj))[3],
      N_tumor = sum(tumor_idx), N_normal = sum(normal_idx)
    ))
  }
}
write.csv(roc_results, file.path(table_dir, "ROC_analysis.csv"), row.names=FALSE)

for(cancer in cancers) {
  cancer_roc <- roc_results[roc_results$Cancer == cancer, ]
  if(nrow(cancer_roc) == 0) next
  
  pdf(file.path(fig_dir, paste0("ROC_curves_", cancer, ".pdf")), width=10, height=10)
  plot(0, 0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="False Positive Rate",
       ylab="True Positive Rate", main=paste(cancer, "- ROC Curves (Tumor vs Normal)"))
  abline(0, 1, lty=2, col="gray")
  
  colors <- rainbow(nrow(cancer_roc))
  expr <- all_data[[cancer]]$expr
  surv_data <- all_data[[cancer]]$survival
  labels <- ifelse(surv_data$sample_type == "Primary Tumor", 1, ifelse(surv_data$sample_type == "Solid Tissue Normal", 0, NA))
  
  for(i in 1:nrow(cancer_roc)) {
    gene <- cancer_roc$Gene[i]
    roc_obj <- roc(labels[!is.na(labels)], expr[gene, !is.na(labels)], quiet=TRUE)
    lines(1 - roc_obj$specificities, roc_obj$sensitivities, col=colors[i], lwd=2)
  }
  legend("bottomright", legend=paste(cancer_roc$Gene, sprintf("(%.2f)", cancer_roc$AUC)), col=colors, lwd=2, cex=0.7, ncol=2)
  dev.off()
}
cat("ROC analysis completed\n\n")

# --- PART 6: RANDOM GENE SET VALIDATION ---
cat("\n--- PART 6: RANDOM GENE SET VALIDATION ---\n")
n_random_sets <- 1000
random_results <- list()

for(cancer in cancers) {
  if(is.null(all_data[[cancer]])) next
  expr <- all_data[[cancer]]$expr
  surv_data <- all_data[[cancer]]$survival
  
  tumor_idx <- surv_data$sample_type == "Primary Tumor"
  normal_idx <- surv_data$sample_type == "Solid Tissue Normal"
  
  random_cors <- numeric(n_random_sets)
  random_aucs <- numeric(n_random_sets)
  
  for(i in 1:n_random_sets) {
    random_genes <- sample(rownames(expr), 15)
    expr_subset <- expr[random_genes, tumor_idx]
    cor_mat <- cor(t(expr_subset), use="pairwise.complete.obs")
    random_cors[i] <- mean(abs(cor_mat[upper.tri(cor_mat)]))
    
    if(sum(normal_idx) >= 5) {
      aucs <- numeric(15)
      labels <- ifelse(tumor_idx, 1, ifelse(normal_idx, 0, NA))
      for(j in 1:15) {
        if(sum(!is.na(labels)) >= 10) {
          roc_obj <- roc(labels[!is.na(labels)], expr[random_genes[j], !is.na(labels)], quiet=TRUE)
          aucs[j] <- as.numeric(auc(roc_obj))
        }
      }
      random_aucs[i] <- mean(aucs[aucs > 0])
    }
  }
  random_results[[cancer]] <- list(cors = random_cors, aucs = random_aucs)
}

comparison_table <- data.frame()
for(cancer in cancers) {
  if(is.null(correlation_matrices[[cancer]])) next
  cor_mat <- correlation_matrices[[cancer]]
  obs_cor <- mean(abs(cor_mat[upper.tri(cor_mat)]))
  random_cors <- random_results[[cancer]]$cors
  
  comparison_table <- rbind(comparison_table, data.frame(
    Cancer = cancer, Metric = "Mean_Absolute_Correlation", Observed = obs_cor,
    Random_Mean = mean(random_cors), Random_SD = sd(random_cors),
    P_value = sum(random_cors >= obs_cor) / n_random_sets, Percentile = sum(random_cors < obs_cor) / n_random_sets * 100
  ))
  
  cancer_aucs <- roc_results$AUC[roc_results$Cancer == cancer]
  if(length(cancer_aucs) > 0) {
    obs_auc <- mean(cancer_aucs)
    random_aucs <- random_results[[cancer]]$aucs[random_results[[cancer]]$aucs > 0]
    if(length(random_aucs) > 0) {
      comparison_table <- rbind(comparison_table, data.frame(
        Cancer = cancer, Metric = "Mean_AUC", Observed = obs_auc,
        Random_Mean = mean(random_aucs), Random_SD = sd(random_aucs),
        P_value = sum(random_aucs >= obs_auc) / length(random_aucs), Percentile = sum(random_aucs < obs_auc) / length(random_aucs) * 100
      ))
    }
  }
}
write.csv(comparison_table, file.path(table_dir, "random_gene_set_comparison.csv"), row.names=FALSE)

pdf(file.path(fig_dir, "random_gene_set_comparison.pdf"), width=12, height=8)
par(mfrow=c(2, 3))
for(cancer in cancers) {
  if(is.null(random_results[[cancer]])) next
  random_cors <- random_results[[cancer]]$cors
  if(!is.null(correlation_matrices[[cancer]]) && length(random_cors) > 0 && all(is.finite(random_cors))) {
    cor_mat <- correlation_matrices[[cancer]]
    obs_cor <- mean(abs(cor_mat[upper.tri(cor_mat)]))
    if(is.finite(obs_cor)) {
      hist(random_cors, breaks=50, col="lightblue", main=paste(cancer, "- Correlation"),
           xlab="Mean Absolute Correlation", xlim=range(c(random_cors, obs_cor)))
      abline(v=obs_cor, col="red", lwd=3, lty=2)
      legend("topright", "Observed", col="red", lwd=3, lty=2, cex=0.8)
    }
  }
  if(cancer %in% roc_results$Cancer) {
    random_aucs <- random_results[[cancer]]$aucs
    random_aucs <- random_aucs[random_aucs > 0 & is.finite(random_aucs)]
    if(length(random_aucs) > 0) {
      obs_auc <- mean(roc_results$AUC[roc_results$Cancer == cancer])
      if(is.finite(obs_auc) && all(is.finite(random_aucs))) {
        hist(random_aucs, breaks=50, col="lightgreen", main=paste(cancer, "- AUC"),
             xlab="Mean AUC", xlim=range(c(random_aucs, obs_auc)))
        abline(v=obs_auc, col="red", lwd=3, lty=2)
        legend("topright", "Observed", col="red", lwd=3, lty=2, cex=0.8)
      }
    }
  }
}
dev.off()
cat("Random validation completed\n\n")

saveRDS(list(
  common_genes = common_genes,
  all_data = all_data,
  survival_results = survival_results,
  survival_table = survival_table,
  correlation_matrices = correlation_matrices,
  roc_results = roc_results,
  random_results = random_results,
  comparison_table = comparison_table
), file.path(out_dir, "workspace.rds"))

cat("ANALYSIS COMPLETE!\n")
