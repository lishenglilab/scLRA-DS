#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyverse)
  library(Seurat)
  library(ggVennDiagram)
  library(future)
  library(rtracklayer)
  library(parallel)
  library(optparse)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--seurat_tr_rds"), type="character", help="Path to Seurat transcript RDS file"),
  make_option(c("--seurat_gene_rds"), type="character", help="Path to Seurat gene RDS file"),
  make_option(c("--gtf_file"), type="character", help="Path to GTF file"),
  make_option(c("--dominant_transcripts_dir"), type="character", help="Path to dominant transcripts results directory"),
  make_option(c("--output_dir"), type="character", help="Output directory"),
  make_option(c("--ident_1"), type="character", help="First group identifier for comparison (e.g., AD)"),
  make_option(c("--ident_2"), type="character", help="Second group identifier for comparison (e.g., N or Normal)"),
  make_option(c("--comparison_name"), type="character", help="Comparison name (e.g., AD_vs_Normal)"),
  make_option(c("--logfc_threshold"), type="numeric", help="LogFC threshold (e.g., 1.5 for fold change)"),
  make_option(c("--pval_threshold"), type="numeric", help="P-value threshold for significance (e.g., 0.05)"),
  make_option(c("--min_pct"), type="numeric", help="Minimum percentage of cells expressing the feature (e.g., 0.01)"),
  make_option(c("--n_cores"), type="integer", help="Number of cores for parallel processing (e.g., 10)"),
  make_option(c("--stage_mapping_file"), type="character", help="Path to stage mapping file (TSV format: sample stage)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
required_args <- c("seurat_tr_rds", "seurat_gene_rds", "gtf_file", 
                   "dominant_transcripts_dir", "output_dir", 
                   "ident_1", "ident_2", "comparison_name",
                   "logfc_threshold", "pval_threshold", "min_pct", "n_cores",
                   "stage_mapping_file")

missing_args <- required_args[vapply(
  required_args,
  function(name) is.null(opt[[name]]) || length(opt[[name]]) == 0 ||
    (is.character(opt[[name]]) && !nzchar(opt[[name]])),
  logical(1)
)]
if (length(missing_args) > 0) {
  stop(paste("Missing required arguments:", paste(missing_args, collapse=", ")))
}

# Create output directory
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

# Set parallel options
options(future.globals.maxSize = 80 * 1024^3)
plan('multisession', workers = opt$n_cores)

# ------------------------------------------------------------------
# Part 1: Load and process GTF file
# ------------------------------------------------------------------
cat("Loading GTF file...\n")
gtf <- as.data.frame(import(opt$gtf_file))
gtf <- gtf[gtf$type == 'transcript', ]
id_map <- unique(gtf[, c('gene_id', 'gene_name', 'transcript_id')])

# Fill missing gene names
id_map <- id_map %>%
  group_by(gene_id) %>%
  mutate(gene_name = case_when(
    is.na(gene_name) & all(is.na(gene_name)) ~ gene_id,  
    is.na(gene_name) ~ dplyr::first(gene_name[!is.na(gene_name)]),  
    TRUE ~ gene_name  
  )) %>%
  ungroup()

# ------------------------------------------------------------------
# Part 2: Load stage mapping information
# ------------------------------------------------------------------
cat("Loading stage mapping information...\n")
stage_mapping <- read.table(opt$stage_mapping_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
cat("Stage mapping loaded:\n")
print(stage_mapping)

# ------------------------------------------------------------------
# Part 3: Differential Dominant Transcript Usage (DDTU)
# ------------------------------------------------------------------
cat("Analyzing differential dominant transcript usage...\n")

# Try different path patterns to find the dominant transcripts file
possible_paths <- c(
  file.path(opt$dominant_transcripts_dir, opt$comparison_name, paste0("statistics_", opt$comparison_name, ".txt")),
  file.path(opt$dominant_transcripts_dir, paste0("statistics_", opt$comparison_name, ".txt"))
)

dominant_file <- NULL
for (path in possible_paths) {
  if (file.exists(path)) {
    dominant_file <- path
    cat(paste("Found dominant transcripts file at:", path, "\n"))
    break
  }
}

if (is.null(dominant_file)) {
  cat("Could not find dominant transcripts file. Tried paths:\n")
  for (path in possible_paths) {
    cat(paste("  -", path, "\n"))
  }
  stop(paste("Dominant transcripts file not found for comparison:", opt$comparison_name))
}

dominant_transcripts_diff <- read.table(dominant_file, sep = '\t', header = TRUE)

# Save DDTU results
ddtu_output_file <- file.path(opt$output_dir, 
                              paste0("Differential_Dominant_Transcript_Usage_", 
                                     opt$comparison_name, ".txt"))
write.table(dominant_transcripts_diff, ddtu_output_file, 
            sep = '\t', row.names = FALSE, quote = FALSE)

# Get genes with significant DDTU
dominant_tr_change_ge <- unique(dominant_transcripts_diff$gene_name[dominant_transcripts_diff$type != 'notSig'])
cat(paste("Found", length(dominant_tr_change_ge), "genes with significant DDTU\n"))

# ------------------------------------------------------------------
# Part 4: Differential Transcript Expression (DTE)
# ------------------------------------------------------------------
cat("Analyzing differential transcript expression...\n")
sce.tr <- readRDS(opt$seurat_tr_rds)

# Add stage information to Seurat transcript object
cat("Adding stage information to Seurat transcript object...\n")
# Extract sample name from cell barcode (assuming format: SAMPLE_CELLBARCODE)
cell_barcodes <- colnames(sce.tr)
sample_names <- sub("_[^_]+$", "", cell_barcodes)

# Map sample names to stage using stage_mapping
if (all(sample_names %in% stage_mapping$sample)) {
  sce.tr@meta.data$stage <- stage_mapping$stage[match(sample_names, stage_mapping$sample)]
  cat("Stage information added to transcript Seurat object\n")
} else {
  missing_samples <- setdiff(unique(sample_names), stage_mapping$sample)
  stop(paste("The following samples are not found in stage mapping:", 
             paste(missing_samples, collapse=", ")))
}

DefaultAssay(sce.tr) <- 'RNA'
Idents(sce.tr) <- sce.tr$stage

# Check stage column values
cat(paste("Available stages in transcript data:", paste(unique(sce.tr$stage), collapse=", "), "\n"))

# Check if specified idents exist in the data
if (!(opt$ident_1 %in% unique(sce.tr$stage))) {
  stop(paste("ident_1", opt$ident_1, "not found in Seurat object. Available stages:", 
             paste(unique(sce.tr$stage), collapse=", ")))
}
if (!(opt$ident_2 %in% unique(sce.tr$stage))) {
  stop(paste("ident_2", opt$ident_2, "not found in Seurat object. Available stages:", 
             paste(unique(sce.tr$stage), collapse=", ")))
}

# Perform differential expression
DTE <- FindMarkers(sce.tr, 
                   ident.1 = opt$ident_1, 
                   ident.2 = opt$ident_2,
                   logfc.threshold = 0,
                   min.pct = opt$min_pct)

# Add significance column
DTE$change <- ifelse(DTE$p_val_adj < opt$pval_threshold & DTE$avg_log2FC > log2(opt$logfc_threshold), 'Up',
                     ifelse(DTE$p_val_adj < opt$pval_threshold & DTE$avg_log2FC < -log2(opt$logfc_threshold), 'Down', 'Stable'))

cat(paste("DTE results for", opt$comparison_name, ":\n"))
print(table(DTE$change))

# Add transcript and gene information
DTE$transcript_id <- rownames(DTE)
DTE <- left_join(DTE, unique(id_map[, c('transcript_id', 'gene_name')]), by = "transcript_id")

# Save DTE results
dte_output_file <- file.path(opt$output_dir, 
                             paste0("Differential_Transcript_Expression_", 
                                    opt$comparison_name, ".txt"))
write.table(DTE, dte_output_file, sep = '\t', row.names = FALSE, quote = FALSE)

# Get genes with significant DTE
DTE_gene <- unique(DTE$gene_name[DTE$change != 'Stable'])
cat(paste("Found", length(DTE_gene), "genes with significant DTE\n"))

# ------------------------------------------------------------------
# Part 5: Differential Gene Expression (DGE)
# ------------------------------------------------------------------
cat("Analyzing differential gene expression...\n")
sce.gene <- readRDS(opt$seurat_gene_rds)

# Add stage information to Seurat gene object
cat("Adding stage information to Seurat gene object...\n")
# Extract sample name from cell barcode (assuming format: SAMPLE_CELLBARCODE)
cell_barcodes <- colnames(sce.gene)
sample_names <- sub("_[^_]+$", "", cell_barcodes)

# Map sample names to stage using stage_mapping
if (all(sample_names %in% stage_mapping$sample)) {
  sce.gene@meta.data$stage <- stage_mapping$stage[match(sample_names, stage_mapping$sample)]
  cat("Stage information added to gene Seurat object\n")
} else {
  missing_samples <- setdiff(unique(sample_names), stage_mapping$sample)
  stop(paste("The following samples are not found in stage mapping:", 
             paste(missing_samples, collapse=", ")))
}

DefaultAssay(sce.gene) <- 'RNA'
Idents(sce.gene) <- sce.gene$stage

# Check stage column values
cat(paste("Available stages in gene data:", paste(unique(sce.gene$stage), collapse=", "), "\n"))

# Perform differential expression
DGE <- FindMarkers(sce.gene, 
                   ident.1 = opt$ident_1, 
                   ident.2 = opt$ident_2,
                   logfc.threshold = 0)

# Add significance column
DGE$change <- ifelse(DGE$p_val_adj < opt$pval_threshold & DGE$avg_log2FC > log2(opt$logfc_threshold), 'Up',
                     ifelse(DGE$p_val_adj < opt$pval_threshold & DGE$avg_log2FC < -log2(opt$logfc_threshold), 'Down', 'Stable'))

cat(paste("DGE results for", opt$comparison_name, ":\n"))
print(table(DGE$change))

DGE$gene <- rownames(DGE)

# Save DGE results
dge_output_file <- file.path(opt$output_dir, 
                             paste0("Differential_Gene_Expression_", 
                                    opt$comparison_name, ".txt"))
write.table(DGE, dge_output_file, sep = '\t', row.names = FALSE, quote = FALSE)

# Get genes with significant DGE
DGE_gene <- DGE$gene[DGE$change != 'Stable']
cat(paste("Found", length(DGE_gene), "genes with significant DGE\n"))

# ------------------------------------------------------------------
# Part 6: Differential Gene Usage (DGU)
# ------------------------------------------------------------------
cat("Analyzing differential gene usage...\n")
DefaultAssay(sce.gene) <- 'RNA'

# Extract expression matrix
exp <- as.data.frame(sce.gene@assays$RNA@layers$data)
rownames(exp) <- rownames(sce.gene)
colnames(exp) <- colnames(sce.gene)

# Get metadata
meta <- sce.gene@meta.data
meta$cell_id <- rownames(meta)

# Get cell IDs by stage
ident_2_cells <- meta$cell_id[meta$stage == opt$ident_2]
ident_1_cells <- meta$cell_id[meta$stage == opt$ident_1]

# Function to compute percentage of cells expressing each gene
compute_pct <- function(i) {
  ge <- rownames(exp)[i]
  exp_sub_ident1 <- exp[i, ident_1_cells]
  n_ge_ident1 <- sum(exp_sub_ident1 > 0)
  exp_sub_ident2 <- exp[i, ident_2_cells]
  n_ge_ident2 <- sum(exp_sub_ident2 > 0)
  
  res_sub <- data.frame(
    gene_name = ge,
    ident2_pct = n_ge_ident2 / length(ident_2_cells),
    ident1_pct = n_ge_ident1 / length(ident_1_cells)
  )
  
  return(res_sub)
}

# Parallel computation
cat(paste("Computing gene usage with", opt$n_cores, "cores...\n"))

# Use mclapply for parallel processing
res_list <- mclapply(1:nrow(exp), compute_pct, mc.cores = opt$n_cores)
gene_pct <- do.call(rbind, res_list)

# Calculate log fold change
gene_pct$lgFC <- log2((gene_pct$ident1_pct + 0.000001) / (gene_pct$ident2_pct + 0.000001))
gene_pct <- na.omit(gene_pct)

# Classify gene usage changes
gene_pct$type <- ifelse(gene_pct$lgFC > 1, opt$ident_1,
                        ifelse(gene_pct$lgFC < -1, opt$ident_2, 'notSig'))

cat(paste("DGU results for", opt$comparison_name, ":\n"))
print(table(gene_pct$type))

# Save DGU results
dgu_output_file <- file.path(opt$output_dir, 
                             paste0("Differential_Gene_Usage_", 
                                    opt$comparison_name, ".txt"))
write.table(gene_pct, dgu_output_file, sep = '\t', row.names = FALSE, quote = FALSE)

# Get genes with significant DGU
ge_pct_change_ge <- gene_pct$gene_name[gene_pct$type != 'notSig']
cat(paste("Found", length(ge_pct_change_ge), "genes with significant DGU\n"))

# ------------------------------------------------------------------
# Part 7: Generate Venn Diagram
# ------------------------------------------------------------------
cat("Generating Venn diagram...\n")

gene_list <- list(
  "DGU (Differential Gene Usage)" = ge_pct_change_ge,
  "DGE (Differential Gene Expression)" = DGE_gene,
  "DTE (Differential Transcript Expression)" = DTE_gene,
  "DDTU (Differential Dominant Transcript Usage)" = dominant_tr_change_ge
)

# Create Venn diagram
venn_plot <- ggVennDiagram(gene_list) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ggtitle(paste("Gene Overlap -", opt$comparison_name)) +
  theme(plot.title = element_text(hjust = 0.5))

# Save Venn diagram
venn_output_file <- file.path(opt$output_dir, 
                              paste0("Venn_Diagram_", 
                                     opt$comparison_name, ".pdf"))
ggsave(venn_output_file, venn_plot, width = 8, height = 6)
cat(paste("Venn diagram saved to:", venn_output_file, "\n"))

# Also save as PNG
venn_png_file <- file.path(opt$output_dir, 
                           paste0("Venn_Diagram_", 
                                  opt$comparison_name, ".png"))
ggsave(venn_png_file, venn_plot, width = 8, height = 6)
cat(paste("Venn diagram (PNG) saved to:", venn_png_file, "\n"))

# ------------------------------------------------------------------
# Part 8: Create Analysis Summary
# ------------------------------------------------------------------
cat("Creating analysis summary...\n")

summary_stats <- data.frame(
  Analysis = c("Differential Dominant Transcript Usage (DDTU)",
               "Differential Transcript Expression (DTE)",
               "Differential Gene Expression (DGE)",
               "Differential Gene Usage (DGU)"),
  Significant_Genes = c(length(dominant_tr_change_ge),
                        length(DTE_gene),
                        length(DGE_gene),
                        length(ge_pct_change_ge)),
  Upregulated = c(sum(dominant_transcripts_diff$type == opt$ident_1, na.rm = TRUE),
                  sum(DTE$change == 'Up', na.rm = TRUE),
                  sum(DGE$change == 'Up', na.rm = TRUE),
                  sum(gene_pct$type == opt$ident_1, na.rm = TRUE)),
  Downregulated = c(sum(dominant_transcripts_diff$type == opt$ident_2, na.rm = TRUE),
                    sum(DTE$change == 'Down', na.rm = TRUE),
                    sum(DGE$change == 'Down', na.rm = TRUE),
                    sum(gene_pct$type == opt$ident_2, na.rm = TRUE))
)

summary_file <- file.path(opt$output_dir, 
                          paste0("Analysis_Summary_", 
                                 opt$comparison_name, ".txt"))
write.table(summary_stats, summary_file, sep = '\t', row.names = FALSE, quote = FALSE)

cat(paste("Analysis summary saved to:", summary_file, "\n"))
cat("Differential analysis completed successfully!\n")

# Save gene lists for further analysis
gene_lists <- list(
  DDTU_genes = dominant_tr_change_ge,
  DTE_genes = DTE_gene,
  DGE_genes = DGE_gene,
  DGU_genes = ge_pct_change_ge
)

saveRDS(gene_lists, file.path(opt$output_dir, paste0("gene_lists_", opt$comparison_name, ".rds")))
