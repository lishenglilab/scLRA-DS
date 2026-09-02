nextflow.enable.dsl=2

/*
 * Prepare gene expression matrix from IsoQuant output
 */
process PREPARE_GENE_MATRIX {
  label 'sclong'
  tag "prepare_gene_matrix"
  
  publishDir "${params.outdir ?: 'Result'}/PreliminarySingleCellAnalysis/gene_matrix", mode: 'copy'
  
  input:
    path isoquant_dir
  
  output:
    path "gene_matrix", emit: gene_matrix_dir
  
  script:
    """
    # Create output directory
    mkdir -p gene_matrix
    
    # Check if required files exist
    if [[ ! -f "${isoquant_dir}/OUT/OUT.discovered_gene_grouped_counts.barcodes.tsv" ]]; then
      echo "ERROR: Gene barcodes file not found" >&2
      exit 1
    fi
    
    if [[ ! -f "${isoquant_dir}/OUT/OUT.discovered_gene_grouped_counts.features.tsv" ]]; then
      echo "ERROR: Gene features file not found" >&2
      exit 1
    fi
    
    if [[ ! -f "${isoquant_dir}/OUT/OUT.discovered_gene_grouped_counts.matrix.mtx" ]]; then
      echo "ERROR: Gene matrix file not found" >&2
      exit 1
    fi
    
    # Copy and rename files
    cp "${isoquant_dir}/OUT/OUT.discovered_gene_grouped_counts.barcodes.tsv" gene_matrix/barcodes.tsv
    cp "${isoquant_dir}/OUT/OUT.discovered_gene_grouped_counts.features.tsv" gene_matrix/genes.tsv
    cp "${isoquant_dir}/OUT/OUT.discovered_gene_grouped_counts.matrix.mtx" gene_matrix/matrix.mtx
    
    echo "Gene matrix prepared in: gene_matrix"
    echo "Files:"
    ls -la gene_matrix/
    """
}

/*
 * Prepare transcript expression matrix from IsoQuant output
 */
process PREPARE_TRANSCRIPT_MATRIX {
  label 'sclong'
  tag "prepare_transcript_matrix"
  
  publishDir "${params.outdir ?: 'Result'}/PreliminarySingleCellAnalysis/transcript_matrix", mode: 'copy'
  
  input:
    path isoquant_dir
  
  output:
    path "transcript_matrix", emit: transcript_matrix_dir
  
  script:
    """
    # Create output directory
    mkdir -p transcript_matrix
    
    # Check if required files exist
    if [[ ! -f "${isoquant_dir}/OUT/OUT.discovered_transcript_grouped_counts.barcodes.tsv" ]]; then
      echo "ERROR: Transcript barcodes file not found" >&2
      exit 1
    fi
    
    if [[ ! -f "${isoquant_dir}/OUT/OUT.discovered_transcript_grouped_counts.features.tsv" ]]; then
      echo "ERROR: Transcript features file not found" >&2
      exit 1
    fi
    
    if [[ ! -f "${isoquant_dir}/OUT/OUT.discovered_transcript_grouped_counts.matrix.mtx" ]]; then
      echo "ERROR: Transcript matrix file not found" >&2
      exit 1
    fi
    
    # Copy and rename files
    cp "${isoquant_dir}/OUT/OUT.discovered_transcript_grouped_counts.barcodes.tsv" transcript_matrix/barcodes.tsv
    cp "${isoquant_dir}/OUT/OUT.discovered_transcript_grouped_counts.features.tsv" transcript_matrix/genes.tsv
    cp "${isoquant_dir}/OUT/OUT.discovered_transcript_grouped_counts.matrix.mtx" transcript_matrix/matrix.mtx
    
    echo "Transcript matrix prepared in: transcript_matrix"
    echo "Files:"
    ls -la transcript_matrix/
    """
}

/*
 * Run single-cell analysis with R
 */
process RUN_SC_ANALYSIS {
  label 'sclong'
  tag "sc_analysis"
  
  publishDir "${params.outdir ?: 'Result'}/PreliminarySingleCellAnalysis", mode: 'copy'
  
  input:
    path gene_matrix_dir
    path transcript_matrix_dir
    path pipeline_bin
  
  output:
    path "sce.gene.rds", emit: seurat_gene
    path "sce.tr.rds", emit: seurat_tr
    path "*.pdf", emit: plots
    path "*.txt", emit: marker_tables
    path "qc_metrics.tsv", emit: qc_metrics
    path "*.csv", optional: true, emit: annotation_files
  
  script:
    """
    # Run R analysis script
    "${params.rscript_bin ?: 'Rscript'}" "${pipeline_bin}/preliminary_single_cell_analysis.R" \
      --gene_matrix "${gene_matrix_dir}" \
      --transcript_matrix "${transcript_matrix_dir}" \
      --output_dir "." \
      --species "${params.species ?: 'Human'}" \
      --tissue "${params.tissue ?: 'Skin'}" \
      --project_dir "." \
      --min_cells ${params.min_cells ?: 10} \
      --min_features ${params.min_features ?: 200} \
      --llm_enabled "${params.llm_enabled}" \
      --llm_base_url "${params.llm_base_url}" \
      --llm_model "${params.llm_model}"
    
    echo "Single-cell analysis completed"
    echo "Generated files:"
    ls -la *.rds *.pdf *.txt *.tsv *.csv 2>/dev/null || true
    """
}

/*
 * Generate an optional single-cell QC report using an LLM API
 */
process LLM_SC_QC_REPORT {
  label 'preprocess'
  tag "llm_sc_qc_report"
  
  // Keep the historical directory name so existing result consumers do not break.
  publishDir "${params.outdir ?: 'Result'}/PreliminarySingleCellAnalysis/deepseek_qc", mode: 'copy'
  
  input:
    path qc_metrics
    path pipeline_bin
  
  output:
    path "sc_qc_report.txt", emit: sc_qc_report

  when:
    params.llm_enabled
  
  script:
    """
    # Check if qc_metrics file exists
    if [[ ! -f "${qc_metrics}" ]]; then
      echo "ERROR: QC metrics file not found: ${qc_metrics}" >&2
      exit 1
    fi
    
    # Run Python script to generate QC report
    python "${pipeline_bin}/llm/cell_qc.py" \
      --csv_file_path "${qc_metrics}" \
      --output_text_file "sc_qc_report.txt" \
      --species "${params.species ?: 'Human'}" \
      --tissue "${params.tissue ?: 'Skin'}" \
      --base_url "${params.llm_base_url}" \
      --model "${params.llm_model}"
    
    # Check if report was generated
    if [[ ! -s "sc_qc_report.txt" ]]; then
      echo "WARNING: Single-cell QC report is empty" >&2
    else
      echo "Single-cell QC report generated successfully"
      echo "Report preview:"
      head -50 sc_qc_report.txt
    fi
    """
}

/*
 * Main workflow for single-cell analysis
 */
workflow PreliminarySingleCellAnalysis {
  take:
    isoquant_dir
  
  main:
    pipeline_bin = file("${projectDir}/bin", checkIfExists: true)

    // Prepare gene and transcript matrices
    gene_matrix = PREPARE_GENE_MATRIX(isoquant_dir)
    transcript_matrix = PREPARE_TRANSCRIPT_MATRIX(isoquant_dir)
    
    // Run single-cell analysis
    sc_results = RUN_SC_ANALYSIS(
      gene_matrix.gene_matrix_dir,
      transcript_matrix.transcript_matrix_dir,
      pipeline_bin
    )
    
    // Generate QC report (optional, requires API key)
    qc_report = LLM_SC_QC_REPORT(sc_results.qc_metrics, pipeline_bin)
  
  emit:
    seurat_gene = sc_results.seurat_gene     // sce.gene.rds
    seurat_tr   = sc_results.seurat_tr       // sce.tr.rds
    plots = sc_results.plots
    marker_tables = sc_results.marker_tables
    qc_metrics = sc_results.qc_metrics
    annotation_files = sc_results.annotation_files  // optional
    llm_qc_report = qc_report.sc_qc_report
    // Backward-compatible output name.
    sc_qc_report = qc_report.sc_qc_report
}
