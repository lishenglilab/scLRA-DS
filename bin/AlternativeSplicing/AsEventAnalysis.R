#!/usr/bin/env Rscript
# Alternative Splicing Event Analysis

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggsci)
  library(purrr)
  library(forcats)
  library(dplyr)
  library(tidyr)
  library(UpSetR)
  library(optparse)
})

# Define command line options
option_list <- list(
  make_option(c("-p", "--psi_file"), type="character", default=NULL,
              help="Path to merged PSI file", metavar="FILE"),
  make_option(c("-i", "--ioe_file"), type="character", default=NULL,
              help="Path to IOE events file", metavar="FILE"),
  make_option(c("-s", "--seurat_rds"), type="character", default=NULL,
              help="Path to Seurat RDS file", metavar="FILE"),
  make_option(c("-n", "--novel_vs_known"), type="character", default=NULL,
              help="Path to novel_vs_known file", metavar="FILE"),
  make_option(c("-o", "--output_dir"), type="character", default="./ASE_Analysis",
              help="Output directory [default: %default]", metavar="DIR"),
  make_option(c("--min_cell_prop"), type="numeric", default=0.01,
              help="Minimum proportion of cells expressing ASE [default: %default]",
              metavar="NUM")
)

# Parse arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$psi_file) || is.null(opt$ioe_file) || 
    is.null(opt$seurat_rds) || is.null(opt$novel_vs_known)) {
  print_help(opt_parser)
  stop("Missing required arguments!")
}

# Assign variables
psi_file <- opt$psi_file
ioe_file <- opt$ioe_file
seurat_rds <- opt$seurat_rds
novel_vs_known_file <- opt$novel_vs_known
output_dir <- opt$output_dir
min_cell_prop <- opt$min_cell_prop

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Set output file paths
output_plot_type_dist <- file.path(output_dir, "ASE_type_distribution.pdf")
output_plot_transcript_cat <- file.path(output_dir, "ASE_transcript_structural_category.pdf")
output_plot_prevalence <- file.path(output_dir, "ASE_prevalence_distribution.pdf")
output_plot_sample_analysis <- file.path(output_dir, "ASE_sample_analysis.pdf")
output_stats <- file.path(output_dir, "ASE_statistics.txt")

# Log start
message("Starting ASE analysis...")
message(paste("Input PSI file:", psi_file))
message(paste("Input IOE file:", ioe_file))
message(paste("Input Seurat RDS:", seurat_rds))
message(paste("Input novel_vs_known file:", novel_vs_known_file))
message(paste("Output directory:", output_dir))
message(paste("Min cell proportion:", min_cell_prop))

# Function to log messages
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(paste0("[", timestamp, "] ", msg))
  write(paste0("[", timestamp, "] ", msg), 
        file = file.path(output_dir, "analysis.log"), append = TRUE)
}

log_message("Loading data...")

# Load PSI data
psi <- fread(psi_file, sep = '\t', data.table = FALSE)
if (ncol(psi) < 2) {
  stop("PSI input must contain an ASE column and at least one cell column")
}
log_message(paste("Loaded PSI data with", nrow(psi), "events and", ncol(psi)-1, "cells"))

# Filter ASE events (expressed in at least min_cell_prop of cells)
ase_keep <- rowSums(psi[, -1, drop = FALSE] > 0, na.rm = TRUE) >=
  min_cell_prop * (ncol(psi) - 1)
psi_filtered <- psi[ase_keep, ]
log_message(paste("Filtered to", nrow(psi_filtered), 
                  "ASE events expressed in >=", min_cell_prop*100, "% of cells"))

if (nrow(psi_filtered) == 0) {
  writeLines(
    c(
      "=== ASE ANALYSIS STATISTICS ===",
      "",
      paste("Total_ASE_events:", nrow(psi)),
      "Filtered_ASE_events: 0",
      paste("Number_of_cells:", ncol(psi) - 1),
      paste("Min_cell_proportion:", min_cell_prop),
      "No ASE events passed the configured prevalence threshold."
    ),
    output_stats
  )
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
  writeLines(capture.output(opt), file.path(output_dir, "command_args.txt"))
  log_message("No ASE events passed filtering; wrote an empty-result summary.")
  quit(save = "no", status = 0)
}

# Load Seurat object and get cell metadata
log_message("Loading Seurat object...")
sce.tr <- readRDS(seurat_rds)
col_data <- data.frame(cell_id = colnames(psi_filtered)[-1])
col_data$sample <- sub('_[^_]+$', '', col_data$cell_id)
col_data <- col_data[col_data$cell_id %in% colnames(sce.tr), ]
if (nrow(col_data) == 0) {
  stop("No PSI cell IDs matched the Seurat object")
}
psi_filtered <- psi_filtered[, c('ASE', col_data$cell_id)]
log_message(paste("Matched", nrow(col_data), "cells between PSI and Seurat data"))

# Load event annotation
ioe <- fread(ioe_file, data.table = FALSE)

# Extract event types
ASE <- data.frame(ASE = psi_filtered[,1])
ASE$type <- str_split(str_split(ASE$ASE, pattern = ";", simplify = TRUE)[,2], 
                      pattern = ":", simplify = TRUE)[,1]
ASE <- left_join(ASE, ioe, by = c("ASE" = "event_id"))

# ----------------------------------------------------------------------
# Plot 1: ASE type distribution (pie chart)
# ----------------------------------------------------------------------
log_message("Creating ASE type distribution plot...")
df <- as.data.frame(table(ASE$type))
colnames(df) <- c("Type", "Num")
df <- df %>%
  mutate(Type_with_num = paste0(Type, " (", Num, ")"))

p_type_dist <- ggplot(df, aes(x = "", y = Num, fill = Type_with_num)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar(theta = "y") +
  scale_fill_nejm(name = "Event Type (Count)") +
  theme_void() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14)) +
  ggtitle("ASE Type Distribution")

ggsave(plot = p_type_dist, filename = output_plot_type_dist, width = 8, height = 6)
log_message(paste("Saved type distribution plot to:", output_plot_type_dist))

# ----------------------------------------------------------------------
# Plot 2: Transcript structural category by ASE type
# ----------------------------------------------------------------------
log_message("Creating transcript structural category plot...")
split_rows <- str_split(as.character(ASE$alternative_transcripts), pattern = ",")
df_b <- data.frame(
  alternative_transcripts = unlist(split_rows),
  Type = rep(ASE$type, sapply(split_rows, length))
)

# Load structural category data
if (file.exists(novel_vs_known_file)) {
  first_field <- scan(
    novel_vs_known_file,
    what = character(),
    sep = "\t",
    nmax = 1,
    quiet = TRUE
  )
  has_header <- length(first_field) == 1 && sub("^#", "", first_field) == "isoform"
  structural_category <- fread(
    novel_vs_known_file,
    sep = "\t",
    header = has_header
  )
  if (ncol(structural_category) >= 6) {
    structural_category <- structural_category[, c(1, 6)]
    colnames(structural_category) <- c('isoform', 'structural_category')
    df_b <- left_join(df_b, structural_category, 
                      by = c("alternative_transcripts" = "isoform"))
    
    df_b$structural_category[is.na(df_b$structural_category)] <- 'full_splice_match'
    category_select <- c("full_splice_match", "incomplete_splice_match", 
                         "novel_in_catalog", "novel_not_in_catalog")
    df_b$structural_category[!(df_b$structural_category %in% category_select)] <- "other"
    
    pt_data <- df_b %>% 
      group_by(structural_category, Type) %>% 
      summarise(Num = n(), .groups = 'drop')
    
    p_transcript_cat <- ggplot(pt_data, aes(x = Type, y = Num, fill = structural_category)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_nejm(name = "Transcript Category") +
      labs(x = "ASE Type", y = "Count", 
           title = "Transcript Structural Category by ASE Type") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(plot = p_transcript_cat, filename = output_plot_transcript_cat, 
           width = 10, height = 6)
    log_message(paste("Saved transcript category plot to:", output_plot_transcript_cat))
  } else {
    log_message("Warning: novel_vs_known file has insufficient columns, skipping transcript category plot")
  }
} else {
  log_message(paste("Warning: File not found:", novel_vs_known_file, 
                    "- skipping transcript category plot"))
}

# ----------------------------------------------------------------------
# Plot 3: ASE prevalence across samples
# ----------------------------------------------------------------------
log_message("Creating ASE prevalence plot...")
samples <- unique(col_data$sample)
res_list <- lapply(samples, function(s) {
  cell_ids <- col_data$cell_id[col_data$sample == s]
  psi_sub <- psi_filtered[, c("ASE", cell_ids), drop = FALSE]
  data.frame(
    sample = s,
    ASE = psi_sub$ASE,
    proportion = rowMeans(psi_sub[-1] > 0)
  )
})

res <- do.call(rbind, res_list)

# Create combined plot
p_prevalence_density <- ggplot(res, aes(x = proportion)) +
  geom_density(fill = "steelblue", alpha = 0.7) +
  facet_wrap(~sample) +
  labs(x = "Proportion of expressing cells", y = "Density", 
       title = "ASE Prevalence Distribution by Sample") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

p_prevalence_box <- ggplot(res, aes(x = sample, y = proportion)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  labs(x = "Sample", y = "Proportion of expressing cells",
       title = "ASE Prevalence by Sample") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Save combined plot
pdf(output_plot_prevalence, width = 12, height = 8)
print(p_prevalence_density)
print(p_prevalence_box)
dev.off()
log_message(paste("Saved prevalence plot to:", output_plot_prevalence))

# ----------------------------------------------------------------------
# Plot 4: Sample-level ASE analysis (stacked bar + upset plot)
# ----------------------------------------------------------------------
log_message("Creating sample-level ASE analysis...")
sample_ids <- col_data$sample
unique_samples <- unique(sample_ids)

# Calculate ASE statistics per sample
ase_stats <- lapply(unique_samples, function(sample) {
  sample_cols <- which(sample_ids == sample)
  
  # Get ASEs expressed in > min_cell_prop of cells for this sample
  expressed_ases <- rowSums(
    psi_filtered[, sample_cols + 1L, drop = FALSE] > 0,
    na.rm = TRUE
  ) > min_cell_prop * length(sample_cols)
  ase_names <- psi_filtered[[1]][expressed_ases]
  
  list(
    total_ases = sum(expressed_ases),
    type_counts = table(
      str_extract(ase_names, "(?<=;)[A-Z0-9]+(?=:)")
    ),
    ase_names = ase_names
  )
}) %>% setNames(unique_samples)

# Convert to data frame
ase_stats_df <- map_df(ase_stats, ~{
  if (length(.x$type_counts) > 0) {
    tibble(
      Type = names(.x$type_counts),
      Count = as.numeric(.x$type_counts),
      Proportion = Count / .x$total_ases
    )
  } else {
    tibble(Type = character(), Count = numeric(), Proportion = numeric())
  }
}, .id = "Sample")

# Order samples by total ASE count
ase_stats_df <- ase_stats_df %>%
  mutate(Sample = fct_reorder(Sample, Count, .fun = sum, .desc = TRUE))

# Plot 4a: Stacked bar chart
p_sample_stacked <- ggplot(ase_stats_df, aes(x = Sample, y = Count, fill = Type)) +
  geom_col(position = "stack") +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "ASE Type Distribution by Sample",
       x = "Sample", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Prepare data for upset plot
upset_input <- lapply(ase_stats, function(x) x$ase_names)
upset_data <- fromList(upset_input)

# Function for accurate ASE intersection analysis
get_upsetr_accurate_ases <- function(new_data, set_names) {
  all_combinations <- map(length(set_names):1, ~combn(set_names, .x, simplify = FALSE)) %>% 
    unlist(recursive = FALSE) %>% 
    set_names(map_chr(., ~paste(.x, collapse = "&")))
  
  ase_assignments <- setNames(rep("", nrow(new_data)), rownames(new_data))
  
  for (comb in all_combinations) {
    comb_name <- paste(comb, collapse = "&")
    ase_in_comb <- rownames(new_data)[
      rowSums(new_data[, comb, drop = FALSE]) == length(comb) & 
        ase_assignments == ""
    ]
    ase_assignments[ase_in_comb] <- comb_name
  }
  
  list(
    counts = table(ase_assignments) %>% sort(decreasing = TRUE),
    ase_lists = split(names(ase_assignments), ase_assignments) %>% 
      map(~list(
        ases = .x,
        types = table(str_extract(.x, "(?<=;)[A-Z0-9]+(?=:)"))
      ))
  )
}

# Prepare matrix for intersection analysis
set_names <- names(upset_input)
ase_names <- unique(unlist(upset_input))
new_data <- matrix(0, nrow = length(ase_names), ncol = length(set_names),
                   dimnames = list(ase_names, set_names))

for (sample in set_names) {
  new_data[upset_input[[sample]], sample] <- 1
}

# Get intersection analysis
intersect_ases <- get_upsetr_accurate_ases(new_data, set_names)

# Prepare data for intersection type composition plot
plot_data <- map_df(names(intersect_ases$ase_lists), function(comb_name) {
  comb_data <- intersect_ases$ase_lists[[comb_name]]
  
  if (length(comb_data$types) > 0) {
    types_df <- tibble(
      Type = names(comb_data$types),
      Count = as.numeric(comb_data$types)
    )
  } else {
    types_df <- tibble(Type = character(), Count = integer())
  }
  
  tibble(
    Combination = comb_name,
    types_df
  )
}) %>%
  filter(Count > 0)

# Define colors for ASE types
type_colors <- c(
  "A3" = "#8DD3C7", "A5" = "#FFFFB3", "AF" = "#BEBADA",
  "AL" = "#FB8072", "MX" = "#80B1D3", "RI" = "#FDB462", "SE" = "#B3DE69"
)

# Plot 4c: Intersection type composition
p_intersection_composition <- ggplot(plot_data, 
                                     aes(x = reorder(Combination, -Count, sum), 
                                         y = Count, fill = Type)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = type_colors) +
  labs(
    title = "ASE Type Composition by Sample Combination",
    x = "Sample Combination",
    y = "Number of ASEs",
    fill = "Splicing Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  scale_x_discrete(limits = names(sort(intersect_ases$counts, decreasing = TRUE)))

# Save all sample analysis plots
pdf(output_plot_sample_analysis, width = 14, height = 10)
print(p_sample_stacked)

# Upset plot
upset(upset_data, 
      nsets = length(unique_samples),
      nintersects = 20,
      order.by = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "darkred")

print(p_intersection_composition)
dev.off()
log_message(paste("Saved sample analysis plot to:", output_plot_sample_analysis))

# ----------------------------------------------------------------------
# Save statistics
# ----------------------------------------------------------------------
log_message("Saving statistics...")

stats_summary <- data.frame(
  Total_ASE_events = nrow(psi),
  Filtered_ASE_events = nrow(psi_filtered),
  Number_of_cells = ncol(psi) - 1,
  Number_of_samples = length(unique_samples),
  Min_cell_proportion = min_cell_prop,
  Analysis_date = as.character(Sys.time())
)

# Add type counts
type_counts <- as.data.frame(table(ASE$type))
colnames(type_counts) <- c("ASE_Type", "Count")
type_counts$Percentage <- round(type_counts$Count / sum(type_counts$Count) * 100, 2)

# Write statistics
sink(output_stats)
cat("=== ASE ANALYSIS STATISTICS ===\n\n")
cat("Summary:\n")
print(stats_summary)
cat("\nASE Type Distribution:\n")
print(type_counts)
cat("\nSample-level Statistics:\n")
for (sample in names(ase_stats)) {
  cat(paste0("\n", sample, ":\n"))
  cat(paste0("  Total ASEs: ", ase_stats[[sample]]$total_ases, "\n"))
  if (length(ase_stats[[sample]]$type_counts) > 0) {
    for (type in names(ase_stats[[sample]]$type_counts)) {
      cat(paste0("  ", type, ": ", ase_stats[[sample]]$type_counts[type], "\n"))
    }
  }
}
sink()

log_message(paste("Saved statistics to:", output_stats))

# Save session info
writeLines(capture.output(sessionInfo()), 
           file.path(output_dir, "session_info.txt"))

# Save command line arguments
writeLines(capture.output(opt), 
           file.path(output_dir, "command_args.txt"))

log_message("ASE analysis completed successfully!")
