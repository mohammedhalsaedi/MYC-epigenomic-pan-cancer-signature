# ==============================================================================
# Script: 01_dmrcate_methylation_tcga.R
# Purpose: Identification of upregulated and hypomethylated genes in MYC-high 
#          tumours across TCGA cohorts using limma and DMRcate.
# Inputs:  TCGA GDC unstranded RNA-seq counts, Illumina 450K methylation beta values,
#          and TCGAbiolinks Tumor.purity estimates.
# Outputs: Differentially expressed genes, hypomethylated DMRs, candidate overlap lists,
#          and functional enrichment (GO/KEGG) tables.
# ==============================================================================

cancer <- "COAD"    # Change: BRCA, COAD, LIHC

# --- PARAMETERS ---
FC_THRESHOLD <- 0.5       
METH_FDR <- 0.05          
PCUTOFF <- 0.05           
GENE_WINDOW <- 10000      
GO_PVAL <- 0.05

cat("\nTCGA MYC - DMRcate APPROACH:", cancer, "\n")
cat("FC >=", FC_THRESHOLD, "| Meth FDR =", METH_FDR, "| Gene window = ±", GENE_WINDOW/1000, "kb\n\n")

# --- Setup ---
pkgs <- c("data.table","httr","jsonlite","TCGAbiolinks","limma","edgeR",
          "DMRcate","minfi","IlluminaHumanMethylation450kanno.ilmn12.hg19",
          "ggplot2","pheatmap","RColorBrewer","org.Hs.eg.db",
          "clusterProfiler","enrichplot","DOSE",
          "TxDb.Hsapiens.UCSC.hg19.knownGene","GenomicRanges")

miss <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(miss)) {
  if(!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org")
  BiocManager::install(miss, ask=FALSE, update=FALSE)
}
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only=TRUE)))

options(EXPERIMENT_HUB_ASK=FALSE)
Sys.setenv(EXPERIMENT_HUB_ASK="FALSE")

# Pathing Configuration
base <- getwd()
expr_dir <- file.path(base, "data", paste0(cancer, "_Expression"))
meth_dir <- file.path(base, "data", paste0(cancer, "_Methylation"))
out_dir  <- file.path(base, "results", cancer, "dmrcate_approach")
plot_dir <- file.path(out_dir, "plots")

if(!dir.exists(expr_dir)) stop("Not found: ", expr_dir)
if(!dir.exists(meth_dir) || length(list.dirs(meth_dir, recursive=FALSE))==0) stop("Empty: ", meth_dir)

dir.create(out_dir,  showWarnings=FALSE, recursive=TRUE)
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

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
  cat("  Loading expression...\n")
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
    
    if(i%%100==0) cat(".")
  }
  cat(" Done (", nrow(mat), "×", ncol(mat), ")\n", sep="")
  
  list(mat=mat, uuids=uuids, genes=rownames(mat))
}

load_meth <- function(dir) {
  cat("  Loading methylation...\n")
  folders <- list.dirs(dir, recursive=FALSE, full.names=TRUE)
  files <- sapply(folders, function(f) {
    x <- list.files(f, pattern="\\.txt(\\.gz)?$", full.names=TRUE)
    if(length(x)) x[1] else NA_character_
  })
  
  valid <- !is.na(files)
  files <- files[valid]
  uuids <- basename(folders)[valid]
  
  cat("  Pass 1")
  all_probes <- character(0)
  for(i in seq_along(files)) {
    ids <- data.table::fread(files[i], header=FALSE, select=1, na.strings=c("NA",""))[[1]]
    all_probes <- union(all_probes, as.character(ids))
    if(i%%100==0) cat(".")
  }
  cat("\n")
  
  mat <- matrix(NA_real_, nrow=length(all_probes), ncol=length(files),
                dimnames=list(all_probes, uuids))
  
  cat("  Pass 2")
  for(i in seq_along(files)) {
    dt <- data.table::fread(files[i], header=FALSE, select=1:2,
                            col.names=c("probe","beta"), na.strings=c("NA",""))
    idx <- match(as.character(dt$probe), all_probes)
    ok <- !is.na(idx)
    mat[idx[ok], i] <- suppressWarnings(as.numeric(dt$beta[ok]))
    if(i%%100==0) cat(".")
  }
  cat(" Done (", nrow(mat), "×", ncol(mat), ")\n", sep="")
  
  list(mat=mat, uuids=uuids, probes=all_probes)
}

# --- STEP 1: LOAD ---
cat("\n--- STEP 1: LOAD ---\n")
E <- load_expr(expr_dir)
M <- load_meth(meth_dir)

# --- STEP 2: MAP ---
cat("\n--- STEP 2: MAP ---\n")
expr_map <- gdc_map(E$uuids)
meth_map <- gdc_map(M$uuids)

for(i in 1:nrow(expr_map)) {
  j <- which(colnames(E$mat)==expr_map$uuid[i])
  if(length(j)) colnames(E$mat)[j] <- expr_map$barcode[i]
}
for(i in 1:nrow(meth_map)) {
  j <- which(colnames(M$mat)==meth_map$uuid[i])
  if(length(j)) colnames(M$mat)[j] <- meth_map$barcode[i]
}

E$mat <- E$mat[, colnames(E$mat) %in% expr_map$barcode, drop=FALSE]
M$mat <- M$mat[, colnames(M$mat) %in% meth_map$barcode, drop=FALSE]

# --- STEP 3: TUMOR ONLY ---
cat("\n--- STEP 3: TUMOR ONLY ---\n")
tumor_samples <- colnames(E$mat)[!grepl("-1[1-3][A-Z]?$", colnames(E$mat))]
E$mat <- E$mat[, tumor_samples, drop=FALSE]
M$mat <- M$mat[, !grepl("-1[1-3][A-Z]?$", colnames(M$mat)), drop=FALSE]

# --- STEP 4: MATCH ---
cat("\n--- STEP 4: MATCH ---\n")
expr_pts <- substr(colnames(E$mat), 1, 12)
meth_pts <- substr(colnames(M$mat), 1, 12)
common_pts <- intersect(expr_pts, meth_pts)

pick <- function(pt, barcodes) {
  matches <- barcodes[substr(barcodes,1,12)==pt]
  if(length(matches)==0) return(NA_character_)
  primary <- matches[grepl("-01[A-Z]?$", matches)]
  if(length(primary)) return(primary[1])
  return(matches[1])
}

expr_sel <- sapply(common_pts, function(p) pick(p, colnames(E$mat)))
meth_sel <- sapply(common_pts, function(p) pick(p, colnames(M$mat)))

valid <- !is.na(expr_sel) & !is.na(meth_sel)
expr_sel <- expr_sel[valid]
meth_sel <- meth_sel[valid]

common_pts <- names(expr_sel)
E$mat <- E$mat[, expr_sel, drop=FALSE]
M$mat <- M$mat[, meth_sel, drop=FALSE]

colnames(E$mat) <- common_pts
colnames(M$mat) <- common_pts
pt_to_bc <- data.frame(patient=common_pts, expr_bc=expr_sel, meth_bc=meth_sel,
                       stringsAsFactors=FALSE)

cat("Aligned:", ncol(E$mat), "patients\n")

# --- STEP 5: MYC GROUPS ---
cat("\n--- STEP 5: MYC GROUPS ---\n")
myc_idx <- which(E$genes=="MYC")
if(length(myc_idx)==0) stop("MYC not found")

myc_expr <- E$mat[myc_idx,]
hi_cut <- quantile(myc_expr, 0.75, na.rm=TRUE)
lo_cut <- quantile(myc_expr, 0.25, na.rm=TRUE)

hi <- names(myc_expr)[myc_expr >= hi_cut]
lo <- names(myc_expr)[myc_expr <= lo_cut]
cat("MYC-high:", length(hi), "| MYC-low:", length(lo), "\n")

grp_samp <- c(hi, lo)
group <- factor(c(rep("hi", length(hi)), rep("lo", length(lo))), levels=c("lo","hi"))

E_grp <- E$mat[, grp_samp, drop=FALSE]
M_grp <- M$mat[, grp_samp, drop=FALSE]

e2f <- c("E2F1","E2F2","E2F3","E2F7","E2F8","CDC6","CDK1","CCNE1","CCNE2","CCNA2",
         "MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","PCNA","RRM1","RRM2","TK1","TYMS",
         "CDC25A","CDC25C","DHFR","POLA1","POLE")
e2f_idx <- which(E$genes %in% e2f)

if(length(e2f_idx)==0) {
  e2f_idx <- order(apply(E_grp,1,var,na.rm=TRUE), decreasing=TRUE)[1:50]
}
prolif <- colMeans(E_grp[e2f_idx,,drop=FALSE], na.rm=TRUE)

data("Tumor.purity", package="TCGAbiolinks", envir=environment())
pur_data <- get("Tumor.purity")
purity <- setNames(rep(NA_real_, length(grp_samp)), grp_samp)

for(i in seq_along(grp_samp)) {
  pt <- grp_samp[i]
  full_bc <- pt_to_bc$expr_bc[pt_to_bc$patient == pt]
  if(length(full_bc)==0) next
  pt_12 <- substr(full_bc, 1, 12)
  hit <- which(rownames(pur_data)==pt_12)
  if(length(hit)) purity[i] <- pur_data$CPE[hit[1]]
}

n_pur <- sum(!is.na(purity))
if(n_pur == 0) {
  use_purity <- FALSE
  complete <- !is.na(prolif)
} else if(n_pur < length(grp_samp) * 0.5) {
  purity[is.na(purity)] <- mean(purity, na.rm=TRUE)
  use_purity <- TRUE
  complete <- !is.na(prolif) & !is.na(purity)
} else {
  use_purity <- TRUE
  complete <- !is.na(prolif) & !is.na(purity)
}

E_cc <- E_grp[, complete, drop=FALSE]
M_cc <- M_grp[, complete, drop=FALSE]
group_cc <- group[complete]
prolif_cc <- prolif[complete]
if(use_purity) purity_cc <- purity[complete]

cat("Complete:", sum(complete), "\n")

# --- STEP 6: DIFFERENTIAL EXPRESSION ---
cat("\n--- STEP 6: DIFFERENTIAL EXPRESSION ---\n")
gene_means <- rowMeans(E_cc, na.rm=TRUE)
E_filt <- E_cc[gene_means > 1,]
E_log <- log2(E_filt + 1)

if(use_purity) {
  design <- model.matrix(~ group_cc + prolif_cc + purity_cc)
} else {
  design <- model.matrix(~ group_cc + prolif_cc)
}

fit <- limma::lmFit(E_log, design)
fit <- limma::eBayes(fit)
de <- limma::topTable(fit, coef=2, number=Inf, sort.by="none")
de$gene <- rownames(de)

myc_high_vs_low_up <- de$gene[de$logFC >= FC_THRESHOLD & de$adj.P.Val < 0.05]
cat("Upregulated (logFC >=", FC_THRESHOLD, "):", length(myc_high_vs_low_up), "\n")
write.csv(de, file.path(out_dir, "DE_all.csv"), row.names=FALSE)

# --- STEP 7: DMRcate ANALYSIS ---
cat("\n--- STEP 7: DMRcate ANALYSIS ---\n")
na_frac <- rowMeans(is.na(M_cc))
M_filt <- M_cc[na_frac < 0.2,]

for(i in 1:nrow(M_filt)) {
  nas <- is.na(M_filt[i,])
  if(any(nas)) M_filt[i,nas] <- mean(M_filt[i,], na.rm=TRUE)
}

anno <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
common_probes <- intersect(rownames(M_filt), rownames(anno))
M_filt <- M_filt[common_probes,]
Mvals <- log2((M_filt+1e-5)/(1-M_filt+1e-5))

cat("Running DMRcate with FDR =", METH_FDR, "\n")
cpg <- tryCatch({
  DMRcate::cpg.annotate(datatype="array", object=Mvals, what="M",
                        arraytype="450K", analysis.type="differential",
                        design=design, coef=2, fdr=METH_FDR)
}, error=function(e) {
  cat("ERROR:", e$message, "\n")
  NULL
})

use_pcutoff <- FALSE
if(!is.null(cpg)) {
  dmrs <- tryCatch({
    DMRcate::dmrcate(cpg, lambda=1000, C=2)
  }, error=function(e) {
    cat("DMRcate failed, using pcutoff =", PCUTOFF, "\n")
    use_pcutoff <- TRUE
    NULL
  })
  
  if(use_pcutoff || is.null(dmrs)) {
    cpg <- DMRcate::cpg.annotate(datatype="array", object=Mvals, what="M",
                                 arraytype="450K", analysis.type="differential",
                                 design=design, coef=2)
    dmrs <- DMRcate::dmrcate(cpg, lambda=1000, C=2, pcutoff=PCUTOFF)
  }
} else {
  cat("Using pcutoff =", PCUTOFF, "\n")
  cpg <- DMRcate::cpg.annotate(datatype="array", object=Mvals, what="M",
                               arraytype="450K", analysis.type="differential",
                               design=design, coef=2)
  dmrs <- DMRcate::dmrcate(cpg, lambda=1000, C=2, pcutoff=PCUTOFF)
}

dmr_gr <- DMRcate::extractRanges(dmrs, genome="hg19")
hypo <- dmr_gr[dmr_gr$meandiff < 0]
cat("DMRs total:", length(dmr_gr), "\n")
cat("DMRs hypomethylated:", length(hypo), "\n")

# --- STEP 8: GENE EXPANSION ---
cat("\n--- STEP 8: GENE MAPPING (±", GENE_WINDOW/1000, "kb) ---\n")
parse_genes <- function(x) {
  g <- unlist(strsplit(as.character(x), ","))
  g <- gsub("\\(.*?\\)", "", g)
  trimws(g[nzchar(trimws(g))])
}

dmr_genes_original <- unique(unlist(lapply(hypo$overlapping.genes, parse_genes)))
cat("Genes from DMR annotations:", length(dmr_genes_original), "\n")

genes_gr <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
hypo_expanded <- resize(hypo, width=width(hypo) + 2*GENE_WINDOW, fix="center")
overlaps <- findOverlaps(hypo_expanded, genes_gr)

gene_entrez <- mcols(genes_gr)$gene_id[subjectHits(overlaps)]
gene_symbols <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys=gene_entrez, 
                        columns="SYMBOL", keytype="ENTREZID")
)

dmr_genes_expanded <- unique(c(
  dmr_genes_original,
  gene_symbols$SYMBOL[!is.na(gene_symbols$SYMBOL)]
))

cat("Genes within ±", GENE_WINDOW/1000, "kb:", length(dmr_genes_expanded), "\n")

if(length(hypo)>0) {
  write.csv(as.data.frame(hypo), file.path(out_dir, "DMR_hypomethylated.csv"), row.names=FALSE)
}

# --- STEP 9: INTERSECTION ---
cat("\n--- STEP 9: INTERSECTION ---\n")
candidates <- intersect(myc_high_vs_low_up, dmr_genes_expanded)
cat("MYC up:", length(myc_high_vs_low_up), "\n")
cat("DMR genes:", length(dmr_genes_expanded), "\n")
cat("FINAL:", length(candidates), "\n")
writeLines(candidates, file.path(out_dir, "final_results.txt"))

# --- STEP 10: GO/PATHWAY ---
cat("\n--- STEP 10: GO & PATHWAY ---\n")
if(length(candidates) < 3) {
  cat("WARNING: <3 candidates - skipping enrichment\n")
  go_bp <- NULL; go_mf <- NULL; kegg <- NULL
} else {
  gene_entrez <- suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys=candidates, 
                          columns="ENTREZID", keytype="SYMBOL")
  )
  gene_entrez <- gene_entrez$ENTREZID[!is.na(gene_entrez$ENTREZID)]
  
  cat("Genes with Entrez IDs:", length(gene_entrez), "\n")
  
  if(length(gene_entrez) < 3) {
    go_bp <- NULL; go_mf <- NULL; kegg <- NULL
  } else {
    go_bp <- tryCatch(enrichGO(gene = gene_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                               pAdjustMethod = "BH", pvalueCutoff = GO_PVAL, qvalueCutoff = GO_PVAL, 
                               readable = TRUE), error=function(e) NULL)
    
    go_mf <- tryCatch(enrichGO(gene = gene_entrez, OrgDb = org.Hs.eg.db, ont = "MF",
                               pAdjustMethod = "BH", pvalueCutoff = GO_PVAL, qvalueCutoff = GO_PVAL,
                               readable = TRUE), error=function(e) NULL)
    
    kegg <- tryCatch(enrichKEGG(gene = gene_entrez, organism = 'hsa',
                                pvalueCutoff = GO_PVAL, pAdjustMethod = "BH"), error=function(e) NULL)
    
    if(!is.null(kegg) && nrow(kegg@result) > 0) {
      kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    }
    
    if(!is.null(go_bp) && nrow(go_bp@result) > 0) {
      write.csv(go_bp@result, file.path(out_dir, "GO_enrichment_BP.csv"), row.names=FALSE)
      cat("  GO BP terms:", nrow(go_bp@result), "\n")
    }
    if(!is.null(go_mf) && nrow(go_mf@result) > 0) {
      write.csv(go_mf@result, file.path(out_dir, "GO_enrichment_MF.csv"), row.names=FALSE)
    }
    if(!is.null(kegg) && nrow(kegg@result) > 0) {
      write.csv(kegg@result, file.path(out_dir, "KEGG_pathways.csv"), row.names=FALSE)
      cat("  KEGG pathways:", nrow(kegg@result), "\n")
    }
  }
}

# --- STEP 11: PLOTS ---
de$col <- ifelse(de$logFC >= FC_THRESHOLD & de$adj.P.Val<0.05, "red",
                 ifelse(de$logFC <= -FC_THRESHOLD & de$adj.P.Val<0.05, "blue", "grey"))
p <- ggplot2::ggplot(de, ggplot2::aes(logFC, -log10(adj.P.Val))) +
  ggplot2::geom_point(ggplot2::aes(color=col), size=0.8, alpha=0.6) +
  ggplot2::scale_color_identity() +
  ggplot2::geom_hline(yintercept=-log10(0.05), lty=2) +
  ggplot2::geom_vline(xintercept=c(-FC_THRESHOLD, FC_THRESHOLD), lty=2) +
  ggplot2::labs(title=paste(cancer, "MYC-high vs MYC-low (DMRcate)"),
                subtitle=paste("Candidates:", length(candidates)),
                x="log2 FC", y="-log10(adj P)") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(plot_dir, "volcano.pdf"), p, width=8, height=6)

if(!is.null(go_bp) && nrow(go_bp@result) > 0) {
  pdf(file.path(plot_dir, "GO_BP_dotplot.pdf"), width=10, height=8)
  print(dotplot(go_bp, showCategory=20))
  dev.off()
}

if(!is.null(kegg) && nrow(kegg@result) > 0) {
  pdf(file.path(plot_dir, "KEGG_dotplot.pdf"), width=10, height=8)
  print(dotplot(kegg, showCategory=20))
  dev.off()
}

pdf(file.path(plot_dir, "MYC_boxplot.pdf"), width=5, height=6)
boxplot(log2(myc_expr[grp_samp]+1) ~ group, col=c("blue","red"), 
        main=paste(cancer, "MYC"), ylab="log2(MYC+1)", las=1, outline=FALSE)
stripchart(log2(myc_expr[grp_samp]+1) ~ group, vertical=TRUE, add=TRUE, 
           method="jitter", pch=20, cex=0.6, col=rgb(0,0,0,0.3))
dev.off()

# --- SUMMARY ---
sink(file.path(out_dir, "summary.txt"))
cat("TCGA MYC - DMRcate APPROACH:", cancer, "\n\n")
cat("PARAMETERS:\n")
cat("  Expr logFC >=", FC_THRESHOLD, "\n")
cat("  Meth FDR =", METH_FDR, "\n")
cat("  Gene window = ±", GENE_WINDOW/1000, "kb\n\n")
cat("RESULTS:\n")
cat("  MYC-high vs low up:", length(myc_high_vs_low_up), "\n")
cat("  DMRs hypomethylated:", length(hypo), "\n")
cat("  Genes near DMRs:", length(dmr_genes_expanded), "\n")
cat("  FINAL CANDIDATES:", length(candidates), "\n\n")
if(!is.null(go_bp) && nrow(go_bp@result) > 0) {
  cat("GO BIOLOGICAL PROCESS (top 10):\n")
  top <- head(go_bp@result[order(go_bp@result$p.adjust), c("Description", "p.adjust", "Count")], 10)
  for(i in 1:nrow(top)) {
    cat("  ", top$Description[i], " (n=", top$Count[i], ")\n", sep="")
  }
  cat("\n")
}
cat("CANDIDATES:\n")
if(length(candidates)) cat(paste(candidates, collapse="\n"), "\n")
sink()

saveRDS(list(cancer=cancer, candidates=candidates, de=de, dmrs=hypo,
             go_bp=go_bp, go_mf=go_mf, kegg=kegg),
        file.path(out_dir, "workspace.rds"))

cat("\nCOMPLETE:", cancer, "- Candidates:", length(candidates), "\n\n")
