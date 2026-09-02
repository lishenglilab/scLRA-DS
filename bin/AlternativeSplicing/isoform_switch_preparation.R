#!/usr/bin/env Rscript
# Isoform Switch Analysis Preparation

suppressPackageStartupMessages({
  library(future)
  library(IsoformSwitchAnalyzeR)
  library(tidyverse)
  library(muscat)
  library(dplyr)
  library(Seurat)
  library(DESeq2)
  library(stringr)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(optparse)
  library(rtracklayer)
  library(GenomicRanges)
})

# Define command line options
option_list <- list(
  make_option(c("-d", "--dominant_results"), type="character", default=NULL,
              help="Path to dominant isoform results directory", metavar="DIR"),
  make_option(c("-s", "--seurat_rds"), type="character", default=NULL,
              help="Path to Seurat RDS file (transcript level)", metavar="FILE"),
  make_option(c("-g", "--seurat_gene_rds"), type="character", default=NULL,
              help="Path to Seurat RDS file (gene level)", metavar="FILE"),
  make_option(c("-t", "--gtf_file"), type="character", default=NULL,
              help="Path to transcript models GTF file", metavar="FILE"),
  make_option(c("-f", "--fasta_file"), type="character", default=NULL,
              help="Path to transcript models FASTA file", metavar="FILE"),
  make_option(c("-o", "--output_dir"), type="character", default="./IsoformSwitch",
              help="Output directory [default: %default]", metavar="DIR"),
  make_option(c("-m", "--stage_mapping"), type="character", default=NULL,
              help="Path to stage mapping file (TSV: sample, stage)", metavar="FILE"),
  make_option(c("-c", "--comparison"), type="character", default="AD_vs_Normal",
              help="Comparison name (e.g., AD_vs_Normal) [default: %default]", metavar="STR"),
  make_option(c("--reference_group"), type="character", default=NULL,
              help="Reference group; avoids inferring it from the comparison name", metavar="STR"),
  make_option(c("--test_group"), type="character", default=NULL,
              help="Test group; avoids inferring it from the comparison name", metavar="STR"),
  make_option(c("--cpc2_result"), type="character", default="",
              help="Path to CPC2 result file (optional)", metavar="FILE"),
  make_option(c("--pfam_result"), type="character", default="",
              help="Path to PFAM result file (optional)", metavar="FILE"),
  make_option(c("--signalp_result"), type="character", default="",
              help="Path to SignalP result file (optional)", metavar="FILE"),
  make_option(c("--iupred2a_result"), type="character", default="",
              help="Path to IUPred2A result file (optional)", metavar="FILE"),
  make_option(c("--n_cores"), type="integer", default=15,
              help="Number of cores for parallel processing [default: %default]", metavar="NUM"),
  make_option(c("--dIF_cutoff"), type="numeric", default=0.1,
              help="dIF cutoff for significant isoform switches [default: %default]", metavar="NUM"),
  make_option(c("--alpha"), type="numeric", default=0.05,
              help="Alpha level for significance [default: %default]", metavar="NUM")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$dominant_results) || is.null(opt$seurat_rds) || 
    is.null(opt$gtf_file) || is.null(opt$fasta_file)) {
  print_help(opt_parser)
  stop("Missing required arguments!")
}

# Assign variables
dominant_results_dir <- opt$dominant_results
seurat_tr_path <- opt$seurat_rds
seurat_gene_path <- opt$seurat_gene_rds
gtf_file <- opt$gtf_file
fasta_file <- opt$fasta_file
output_dir <- opt$output_dir
stage_mapping_file <- opt$stage_mapping
comparison_name <- opt$comparison
reference_group <- opt$reference_group
test_group <- opt$test_group
cpc2_result <- opt$cpc2_result
pfam_result <- opt$pfam_result
signalp_result <- opt$signalp_result
iupred2a_result <- opt$iupred2a_result
n_cores <- opt$n_cores
dIF_cutoff <- opt$dIF_cutoff
alpha <- opt$alpha

# Create output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
consequence_dir <- file.path(output_dir, "consequences")
dir.create(consequence_dir, showWarnings = FALSE, recursive = TRUE)
extdata_dir <- file.path(consequence_dir, "extdata")
dir.create(extdata_dir, showWarnings = FALSE, recursive = TRUE)

# Function to log messages
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(paste0("[", timestamp, "] ", msg))
  write(paste0("[", timestamp, "] ", msg), 
        file = file.path(output_dir, "analysis_prep.log"), append = TRUE)
}

log_message("Starting Isoform Switch Analysis Preparation...")
log_message(paste("Dominant results directory:", dominant_results_dir))
log_message(paste("Seurat transcript RDS:", seurat_tr_path))
log_message(paste("Seurat gene RDS:", ifelse(is.null(seurat_gene_path), "Not provided", seurat_gene_path)))
log_message(paste("GTF file:", gtf_file))
log_message(paste("FASTA file:", fasta_file))
log_message(paste("Output directory:", output_dir))
log_message(paste("Stage mapping file:", ifelse(is.null(stage_mapping_file), "Not provided", stage_mapping_file)))
log_message(paste("Comparison:", comparison_name))
log_message(paste("CPC2 result file:", ifelse(cpc2_result == "", "Not provided", cpc2_result)))
log_message(paste("PFAM result file:", ifelse(pfam_result == "", "Not provided", pfam_result)))
log_message(paste("SignalP result file:", ifelse(signalp_result == "", "Not provided", signalp_result)))
log_message(paste("IUPred2A result file:", ifelse(iupred2a_result == "", "Not provided", iupred2a_result)))
log_message(paste("Number of cores:", n_cores))
log_message(paste("dIF cutoff:", dIF_cutoff))
log_message(paste("Alpha:", alpha))

# Set parallel processing
options(future.globals.maxSize = 80 * 1024^3)
plan(multisession, workers = n_cores)

# ----------------------------------------------------------------------
# Part 1: Load data and prepare for analysis
# ----------------------------------------------------------------------

# Function for DESeq2 analysis
run_deseq2 <- function(expression_matrix, sample_groups, ref) {
  rownames(expression_matrix) <- expression_matrix[,1]
  expression_matrix[,1] <- NULL
  
  if (length(sample_groups) != ncol(expression_matrix)) {
    stop("Sample group length must match expression matrix columns")
  }
  
  # Prepare sample information
  sample_info <- data.frame(
    row.names = colnames(expression_matrix),
    condition = relevel(factor(sample_groups), ref = ref)
  )
  
  # Build DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = expression_matrix,
                                colData = sample_info,
                                design = ~ condition)
  
  # Run differential expression analysis
  dds <- DESeq(dds)
  
  # Get results
  results <- as.data.frame(results(dds))
  results$padj[is.na(results$padj)] <- 1
  return(results)
}

log_message("Loading dominant isoform results...")
# Load dominant isoform results
dominant_file <- file.path(dominant_results_dir, comparison_name, paste0("statistics_", comparison_name, ".txt"))
if (!file.exists(dominant_file)) {
  stop(paste("Dominant results file not found:", dominant_file))
}

dominant_df <- read.table(dominant_file, header = TRUE, sep = "\t")
dominant_df_sig <- dominant_df[dominant_df$type != "notSig", ]

if (nrow(dominant_df_sig) == 0) {
  log_message("No significant dominant isoforms found. Exiting.")
  quit(status = 0)
}

log_message(paste("Found", nrow(dominant_df_sig), "significant dominant isoforms"))

# Load Seurat objects
log_message("Loading Seurat objects...")
sce.tr <- readRDS(seurat_tr_path)

# Subset to significant transcripts
sce.tr <- sce.tr[unique(dominant_df_sig$transcript_id), ]

# Load gene-level Seurat object if provided
if (!is.null(seurat_gene_path) && file.exists(seurat_gene_path)) {
  sce.gene <- readRDS(seurat_gene_path)
  sce.gene <- sce.gene[na.omit(unique(dominant_df_sig$gene_name)), ]
} else {
  log_message("Gene-level Seurat object not provided or not found. Using transcript-level only.")
  sce.gene <- NULL
}

# Load stage mapping
if (!is.null(stage_mapping_file) && file.exists(stage_mapping_file)) {
  log_message("Loading stage mapping...")
  stage_mapping <- read.table(stage_mapping_file, header = TRUE, sep = "\t")
  
  # Add stage information to Seurat objects
  if (!is.null(sce.tr)) {
    # Extract sample name from orig.ident or cell_id
    if ("orig.ident" %in% colnames(sce.tr@meta.data)) {
      sce.tr$sample <- sce.tr$orig.ident
    } else {
      sce.tr$sample <- sub("_[^_]+$", "", colnames(sce.tr))
    }
    
    # Merge stage information
    sce.tr$stage <- stage_mapping$stage[match(sce.tr$sample, stage_mapping$sample)]
    
    # Check for unmapped samples
    unmapped_tr <- sum(is.na(sce.tr$stage))
    if (unmapped_tr > 0) {
      log_message(paste("Warning:", unmapped_tr, "transcripts have unmapped stage information"))
    }
  }
  
  if (!is.null(sce.gene)) {
    if ("orig.ident" %in% colnames(sce.gene@meta.data)) {
      sce.gene$sample <- sce.gene$orig.ident
    } else {
      sce.gene$sample <- sub("_[^_]+$", "", colnames(sce.gene))
    }
    
    sce.gene$stage <- stage_mapping$stage[match(sce.gene$sample, stage_mapping$sample)]
    
    unmapped_gene <- sum(is.na(sce.gene$stage))
    if (unmapped_gene > 0) {
      log_message(paste("Warning:", unmapped_gene, "genes have unmapped stage information"))
    }
  }
} else {
  log_message("No stage mapping file provided. Attempting to infer stages from sample names...")
  
  # Infer stages from sample names
  if (!is.null(sce.tr)) {
    sce.tr$sample <- ifelse("orig.ident" %in% colnames(sce.tr@meta.data), 
                            sce.tr$orig.ident, 
                            sub("_[^_]+$", "", colnames(sce.tr)))
    sce.tr$stage <- gsub("[0-9].*", "", sce.tr$sample)
  }
  
  if (!is.null(sce.gene)) {
    sce.gene$sample <- ifelse("orig.ident" %in% colnames(sce.gene@meta.data), 
                              sce.gene$orig.ident, 
                              sub("_[^_]+$", "", colnames(sce.gene)))
    sce.gene$stage <- gsub("[0-9].*", "", sce.gene$sample)
  }
}

# ----------------------------------------------------------------------
# Part 2: Prepare data for IsoformSwitchAnalyzeR
# ----------------------------------------------------------------------
log_message("Preparing data for IsoformSwitchAnalyzeR...")

# Filter GTF to include only transcripts in our analysis
gtf <- rtracklayer::import(gtf_file)
gtf_df <- as.data.frame(gtf)
filtered_gtf_df <- gtf_df[gtf_df$transcript_id %in% rownames(sce.tr), ]

filtered_gtf <- makeGRangesFromDataFrame(filtered_gtf_df, keep.extra.columns = TRUE)
filtered_gtf_path <- file.path(output_dir, "filtered_transcripts.gtf")
rtracklayer::export(filtered_gtf, filtered_gtf_path)
log_message(paste("Filtered GTF saved to:", filtered_gtf_path))

# Convert Seurat to SingleCellExperiment and aggregate data
DefaultAssay(sce.tr) <- "RNA"
sce.tr <- as.SingleCellExperiment(sce.tr)

# Use sample_id and stage for aggregation
if (!"sample_id" %in% colnames(colData(sce.tr))) {
  colData(sce.tr)$sample_id <- sce.tr$sample
}
if (!"stage" %in% colnames(colData(sce.tr))) {
  colData(sce.tr)$stage <- sce.tr$stage
}

# Remove any NA values
sce.tr <- sce.tr[, !is.na(colData(sce.tr)$stage) & !is.na(colData(sce.tr)$sample_id)]

log_message("Aggregating data by sample...")
pb <- aggregateData(sce.tr,
                    assay = "counts", 
                    fun = "sum",
                    by = c("sample_id"))

isoformCountMatrix <- as.data.frame(pb@assays@data@listData[[1]])
isoformCountMatrix <- rownames_to_column(isoformCountMatrix, var = "isoform_id")

# Load GTF for annotation
log_message("Loading GTF annotation...")
isoformExonAnnotation <- filtered_gtf_path
localSwitchList <- importGTF(pathToGTF = isoformExonAnnotation, 
                             addAnnotatedORFs = FALSE, 
                             removeTECgenes = FALSE, 
                             quiet = TRUE)

isoAnnot <- unique(localSwitchList$isoformFeatures[, c("gene_id", "isoform_id", "gene_name")])
isoAnnot <- isoAnnot %>%
  group_by(gene_id) %>%
  mutate(gene_name = ifelse(all(is.na(gene_name)), 
                            NA, 
                            gene_name[!is.na(gene_name)][1])) %>%
  ungroup()

# Filter isoform count matrix
isoformCountMatrix <- isoformCountMatrix[isoformCountMatrix$isoform_id %in% isoAnnot$isoform_id, ]

# Create gene count matrix
geneCountMatrix <- isoformToGeneExp(isoformCountMatrix, localSwitchList)
geneCountMatrix <- dplyr::select(geneCountMatrix, "gene_id", everything())

# ----------------------------------------------------------------------
# Part 3: Differential Expression Analysis
# ----------------------------------------------------------------------
log_message("Performing differential expression analysis...")

# Extract sample groups from metadata
sample_ids <- colnames(isoformCountMatrix)[-1]
sample_groups <- sapply(sample_ids, function(x) {
  unique(colData(sce.tr)$stage[colData(sce.tr)$sample_id == x])
})

# Prefer explicit group values from the workflow. Keep comparison-name parsing
# as a compatibility fallback for existing direct callers.
if (!is.null(reference_group) && nzchar(reference_group) &&
    !is.null(test_group) && nzchar(test_group)) {
  ref_group <- reference_group
  test_group <- test_group
} else {
  comparison_groups <- strsplit(comparison_name, "_vs_", fixed = TRUE)[[1]]
  if (length(comparison_groups) != 2) {
    stop(paste(
      "Cannot infer groups from comparison name:", comparison_name,
      "Provide --reference_group and --test_group"
    ))
  }
  ref_group <- comparison_groups[1]
  test_group <- comparison_groups[2]
}

log_message(paste("Reference group:", ref_group))
log_message(paste("Test group:", test_group))

# Run DESeq2 for genes and isoforms
DEG <- run_deseq2(geneCountMatrix, sample_groups, ref_group)
DET <- run_deseq2(isoformCountMatrix, sample_groups, ref_group)

# Save DE results
de_output_dir <- file.path(output_dir, "differential_expression")
dir.create(de_output_dir, showWarnings = FALSE)

write.table(DEG, file.path(de_output_dir, "gene_DE_results.txt"), 
            sep = "\t", row.names = TRUE, quote = FALSE)
write.table(DET, file.path(de_output_dir, "transcript_DE_results.txt"), 
            sep = "\t", row.names = TRUE, quote = FALSE)

log_message("Differential expression analysis completed")

# ----------------------------------------------------------------------
# Part 4: Isoform Switch Analysis
# ----------------------------------------------------------------------
log_message("Performing isoform switch analysis...")

# Prepare design matrix
designMatrix <- data.frame(
  sampleID = colnames(isoformCountMatrix)[-1],
  condition = sample_groups
)

comparisonsToMake <- data.frame(
  condition_1 = ref_group,
  condition_2 = test_group
)

# Import data into IsoformSwitchAnalyzeR
aSwitchList <- importRdata(
  isoformCountMatrix = isoformCountMatrix,
  designMatrix = designMatrix,
  isoformExonAnnoation = isoformExonAnnotation,
  comparisonsToMake = comparisonsToMake,
  showProgress = TRUE
)

# ----------------------------------------------------------------------
# Part 5: Analyze ORFs and sequences
# ----------------------------------------------------------------------
log_message("Analyzing ORFs and sequences...")

aSwitchList <- analyzeORF(aSwitchList, BSgenome.Hsapiens.UCSC.hg38)

# Extract sequences
sequence_output_dir <- consequence_dir
aSwitchList <- extractSequence(aSwitchList, 
                               genomeObject = BSgenome.Hsapiens.UCSC.hg38,
                               pathToOutput = sequence_output_dir,
                               outputPrefix = paste0("ALL_", comparison_name),
                               onlySwitchingGenes = FALSE)

# ----------------------------------------------------------------------
# Part 6: Integrate external tool results (optional)
# ----------------------------------------------------------------------
log_message("Integrating external tool results...")

# CPC2 analysis (if result file provided)
if (cpc2_result != "" && file.exists(cpc2_result)) {
  log_message(paste("Loading CPC2 results from:", cpc2_result))
  aSwitchList <- analyzeCPC2(aSwitchList,
                             pathToCPC2resultFile = cpc2_result,
                             removeNoncodinORFs = FALSE)
} else {
  log_message("CPC2 result file not provided or not found. Skipping CPC2 analysis.")
}

# PFAM analysis (if result file provided)
if (pfam_result != "" && file.exists(pfam_result)) {
  log_message(paste("Loading PFAM results from:", pfam_result))
  aSwitchList <- analyzePFAM(aSwitchList,
                             pathToPFAMresultFile = pfam_result)
} else {
  log_message("PFAM result file not provided or not found. Skipping PFAM analysis.")
}

# SignalP analysis (if result file provided)
if (signalp_result != "" && file.exists(signalp_result)) {
  log_message(paste("Loading SignalP results from:", signalp_result))
  aSwitchList <- analyzeSignalP(aSwitchList,
                                pathToSignalPresultFile = signalp_result)
} else {
  log_message("SignalP result file not provided or not found. Skipping SignalP analysis.")
}

# IUPred2A analysis (if result file provided)
if (iupred2a_result != "" && file.exists(iupred2a_result)) {
  log_message(paste("Loading IUPred2A results from:", iupred2a_result))
  aSwitchList <- analyzeIUPred2A(aSwitchList,
                                 pathToIUPred2AresultFile = iupred2a_result)
} else {
  log_message("IUPred2A result file not provided or not found. Skipping IUPred2A analysis.")
}

# ----------------------------------------------------------------------
# Save preparation results
# ----------------------------------------------------------------------
log_message("Saving preparation results...")

# Save switch list for consequences analysis
switchlist_file <- file.path(output_dir, paste0("switchList_prep_", comparison_name, ".rds"))
saveRDS(aSwitchList, switchlist_file)

# Save preparation summary
external_tools_used <- c(
  ifelse(cpc2_result != "" && file.exists(cpc2_result), "CPC2", ""),
  ifelse(pfam_result != "" && file.exists(pfam_result), "PFAM", ""),
  ifelse(signalp_result != "" && file.exists(signalp_result), "SignalP", ""),
  ifelse(iupred2a_result != "" && file.exists(iupred2a_result), "IUPred2A", "")
)
external_tools_used <- external_tools_used[nzchar(external_tools_used)]

summary_stats <- list(
  Total_Isoforms_Analyzed = nrow(aSwitchList$isoformFeatures),
  Total_Genes_Analyzed = length(unique(aSwitchList$isoformFeatures$gene_id)),
  Comparison = comparison_name,
  Reference_Group = ref_group,
  Test_Group = test_group,
  External_Tools_Used = ifelse(
    length(external_tools_used) > 0,
    paste(external_tools_used, collapse = ", "),
    "None"
  ),
  Preparation_Date = as.character(Sys.time())
)

summary_file <- file.path(output_dir, "preparation_summary.txt")
sink(summary_file)
cat("Isoform Switch Analysis Preparation Summary\n")
cat("===========================================\n\n")
for (name in names(summary_stats)) {
  cat(paste0(name, ": ", summary_stats[[name]], "\n"))
}
sink()

# Save session info
writeLines(capture.output(sessionInfo()), 
           file.path(output_dir, "preparation_session_info.txt"))

writeLines(capture.output(opt), 
           file.path(output_dir, "preparation_command_args.txt"))

log_message("Isoform switch analysis preparation completed successfully!")
log_message(paste("Switch list saved to:", switchlist_file))
log_message("Ready for consequences analysis.")
