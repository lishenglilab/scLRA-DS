
# Single-Cell Long-Read Transcriptomics Analysis Pipeline
![image](https://github.com/lishenglilab/Single-Cell-Long-Read-Transcriptomics-Analysis-Pipeline/blob/main/Schematic.png)
## Overview

A comprehensive Nextflow pipeline for analyzing single-cell long-read RNA sequencing data. This pipeline integrates multiple analysis modules to process raw sequencing data from FASTQ files to comprehensive biological interpretation, including alternative splicing analysis, ORF-based clustering, and LLM-assisted biological interpretation.

## Workflow Structure

The pipeline follows a sequential workflow organized into five main modules:

### 1. Input and Preprocessing (`InputAndPreprocess.nf`)
- **Metadata Parsing**: Reads metadata and whitelist files to set up initial input data
- **FastQ Processing**: Prepares sample channels and runs preprocessing steps including:
  - Blaze filtering to refine FastQ files
  - Quality control with NanoPlot
  - AI-assisted QC reports using DeepSeek
- **Barcode Handling**: Extracts and validates cell barcodes from 10x Genomics data
- **UMI Quality Control**: Filters and validates UMIs for accurate single-cell analysis

### 2. Alignment and Transcript Reconstruction (`AlignmentAndTranscriptReconstruction.nf`)
- **Read Alignment**: Aligns sequencing reads to reference genome
- **Transcript Reconstruction**: Reconstructs full-length transcripts using IsoQuant
- **ORF Prediction**: Identifies open reading frames from reconstructed transcripts
- **Quality Metrics**: Generates alignment and reconstruction statistics with AI-based quality assessment

### 3. Preliminary Single-Cell Analysis (`PreliminarySingleCellAnalysis.nf`)
- **Dual-QC Framework**:
  - QC-1: Initial quality assessment and cell filtering
  - QC-2: Comprehensive evaluation of transcript reconstruction
- **Data Processing**:
  - Data normalization and batch correction
  - Dimensionality reduction (PCA, UMAP)
  - Cell clustering and AI-assisted cell type annotation 
  - Marker gene identification


### 4. Alternative Splicing Analysis (`AlternativeSplicing.nf`)
- **Splicing Event Detection**: Identifies alternative splicing patterns
- **Intron Retention Analysis**: Detects and quantifies intron retention events
- **Dominant Transcript Identification**: Identifies major transcript isoforms
- **Isoform Switch Analysis**: Analyzes differential transcript usage and isoform switching between conditions
- **AI-Assisted Isoform Interpretation**: Uses AI model to interpret isoform switch analysis results, providing biological context and mechanistic insights into transcript switching eventsUses DeepSeek models to predict functional consequences of isoform switches, prioritize biologically relevant isoforms, and provide biological context for switching events

### 5. Integrative Analysis (`IntegrativeAnalysis.nf`)
- **Differential Analysis**:
  - Differential gene expression
  - Differential transcript usage
  - ORF-based clustering and expression analysis
- **Multi-Omics Integration**: Combines results from multiple analysis modules
- **Biological Interpretation**: AI-assisted interpretation using DeepSeek models
- **ORF Clustering**: Advanced clustering based on open reading frame sequences

## Quick Start

### Prerequisites

- Nextflow (>=22.10.0)
- Java (>=11)

### Installation

1. Clone the repository:
```bash
git clone https://github.com/lishenglilab/Single-Cell-Long-Read-Transcriptomics-Analysis-Pipeline.git
cd Single-Cell-Long-Read-Transcriptomics-Analysis-Pipeline
```

2. Test the installation:
```bash
nextflow run main.nf --help
```

### Basic Usage
Prepare your input files:
- Raw long-read FASTQ files (PacBio or Oxford Nanopore)
- 10x Genomics barcode metadata
- Reference genome assembly (FASTA format)
- Gene annotation (GTF format)
- Sample metadata CSV file

```nextflow
// Configuration file for Single-Cell Long-Read Transcriptomics Analysis Pipeline
nextflow.enable.dsl = 2

params {
  // Study description for AI-assisted interpretation
  study_description = "Your study description here..."
  
  // Output directory
  outdir = "results"
  
  // DeepSeek API configuration
  api_key = "your_api_key_here"
  species = 'Human'
  tissue = 'Skin'
  
  // Input files
  metadata = "metadata.tsv"
  reference_fa = "reference.fa"
  genedb_gtf = "annotation.gtf"
  whitelist = "barcodes.txt"
  
  // QC parameters
  expect_cells = 300
  min_cells = 10
  min_features = 200
  
  // Analysis parameters
  comparisons = [
    ["Condition_A", "Condition_B", "Condition_A_vs_Condition_B"]
  ]
  
  // AI-assisted features
  enable_ai_cell_annotation = true
  enable_ai_quality_assessment = true
  enable_ai_isoform_interpretation = true
}

// Conda environments for reproducible analysis
conda.enabled = true
conda.cacheDir = "${projectDir}/.conda"

process {
  conda = "${projectDir}/envs/sclong.yml"
}

Run the complete pipeline:
```bash
nextflow run main.nf -c params.config 
```

## Input Requirements

### Required Files
- **Metadata CSV**: Contains at minimum sample and bam columns, plus additional sample information
- **Reference Genome**: FASTA format, indexed if required by alignment tools
- **Gene Annotation**: GTF format with transcript models
- **FastQ Files**: Raw sequencing reads, can be compressed (.fastq.gz)

### File Structure Example
```
input/
├── fastq/
│   ├── sample1.fastq.gz
│   └── sample2.fastq.gz
├── metadata.csv
├── reference.fa
├── reference.fa.fai
└── annotation.gtf
```

### Output Structure
```
results/
├── InputAndPreprocess/
│   ├── filtered_fastq/
│   ├── qc_reports/
│   ├── nanoplot_results/
│   └── deepseek_qc/
│
├── AlignmentAndTranscriptReconstruction/
│   ├── aligned_bam/
│   ├── transcript_models.gtf
│   ├── orf_predictions.fa
│   └── alignment_stats/
│
├── PreliminarySingleCellAnalysis/
│   ├── seurat_gene.rds
│   ├── seurat_tr.rds
│   ├── clustering_results/
│   ├── marker_genes/
│   └── qc_plots/
│
├── AlternativeSplicing/
│   ├── splicing_events/
│   ├── intron_retention/
│   ├── dominant_transcripts/
│   └── differential_splicing/
│
└── IntegrativeAnalysis/
    ├── differential_expression/
    ├── orf_clustering/
    ├── biological_interpretation/
    └── integrated_reports/
```



### Key Parameters
```nextflow
// Basic input parameters
params.metadata = "metadata.csv"
params.reference_fa = "reference.fa"
params.genedb_gtf = "annotation.gtf"
params.whitelist = "barcodes.tsv"

// QC parameters
params.qc_min_cells = 3
params.qc_min_features = 200
params.qc_max_mt_percent = 20

// Single-cell analysis parameters
params.resolution = 0.8
params.pca_dims = 30
params.umap_dims = 2

// Differential analysis
params.comparisons = [["AD", "Normal", "AD_vs_Normal"]]
params.logfc_threshold = 1.5
params.pval_threshold = 0.05

// ORF clustering
params.orf_resolution = 0.4
params.orf_dims = "1:30"

// DeepSeek integration
params.api_key = ""  // Required for AI-assisted interpretation
params.study_description = "Single-cell long-read study"
params.species = "Human"  // For DeepSeek context
params.tissue = "Skin"    // For DeepSeek context
```

## Running the Pipeline

### Complete Analysis
```bash
# Run all modules sequentially
nextflow run main.nf 
```


