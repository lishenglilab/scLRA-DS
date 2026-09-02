nextflow.enable.dsl=2

/*
 * Alternative Splicing Analysis Workflow Module
 */

process SPLIT_EXPRESSION {
  label 'sclong'
  tag "split"
  
  publishDir "${params.outdir ?: 'Result'}/AlternativeSplicing", mode: 'copy'
  
  input:
    path seurat_rds
    path pipeline_bin
  
  output:
    path "split_files/*.txt", emit: split_files 
    path "split_info.txt", emit: split_info
  
  script:
    """
    # Create output directory
    mkdir -p split_files
    
    RSCRIPT_PATH="${pipeline_bin}/AlternativeSplicing/split_expression.R"
    "${params.rscript_bin ?: 'Rscript'}" "\${RSCRIPT_PATH}" \
      --seurat_file "${seurat_rds}" \
      --split_size ${params.split_size ?: 100} \
      --output_dir "split_files"
    
    # Create split info file
    echo "Split completed at: \$(date)" > split_info.txt
    echo "Input: ${seurat_rds}" >> split_info.txt
    echo "Split files created:" >> split_info.txt
    ls -1 split_files/*.txt | wc -l >> split_info.txt
    """
}

process GENERATE_AS_EVENTS {  
  label 'alternativesplicing'
  tag "events"
  
  publishDir "${params.outdir ?: 'Result'}/AlternativeSplicing", mode: 'copy'
  
  input:
    path gtf_file
  
  output:
    path "merged_events.ioe", emit: merged_events
  
  script:
  """
  mkdir -p events_output

  echo "=== Generating SUPPA2 splicing events ==="
  echo "Input GTF file: ${gtf_file}"

  suppa.py generateEvents \
    -i "${gtf_file}" \
    -o events_output/isoform_annotated.events \
    -e SE SS MX RI FL \
    -f ioe \
    2>&1 | tee events_output/generate_events.log

  cd events_output

  ioe_files=\$(ls *.ioe 2>/dev/null)
  if [ -z "\$ioe_files" ]; then
    echo "ERROR: No ioe files generated." >&2
    cat generate_events.log >&2
    exit 1
  fi

  # SUPPA writes the same header to every IOE file. Keep it only once.
  awk 'FNR == 1 && NR != 1 { next } { print }' *.ioe > merged_events.ioe

  if [[ ! -s "merged_events.ioe" ]]; then
    echo "ERROR: Failed to create merged events file" >&2
    exit 1
  fi

  cp *.ioe ../
  cp generate_events.log ../

  cd ..
  echo "Event generation completed"
  """
}

process CALCULATE_PSI {
  label 'alternativesplicing'
  tag { exp_file.baseName }
  
  publishDir "${params.outdir ?: 'Result'}/AlternativeSplicing/psi", mode: 'copy'
  
  input:
    path merged_events
    path exp_file
  
  output:
    path "${exp_file.baseName}.psi", emit: psi_file
    path "${exp_file.baseName}.log", emit: psi_log
  
  script:
    exp_name = exp_file.baseName
    """
    echo "Calculating PSI for: ${exp_name}"
    
    suppa.py psiPerEvent \
      -i "${merged_events}" \
      -e "${exp_file}" \
      -o "${exp_name}_project_events" \
      > "${exp_name}.log" 2>&1
    
    if [[ -f "${exp_name}_project_events.psi" ]]; then
      mv "${exp_name}_project_events.psi" "${exp_name}.psi"
    fi
    """
}

process MERGE_AS_RESULTS {  
  label 'sclong'  // 确保使用包含 R 和 data.table 的环境
  tag "merge_psi"
  
  publishDir "${params.outdir ?: 'Result'}/AlternativeSplicing", mode: 'copy'
  
  input:
    path psi_files
  
  output:
    path "merged_project_events.psi", emit: merged_psi
  
  script:
    """
    "${params.rscript_bin ?: 'Rscript'}" --vanilla - << 'RSCRIPT'
    

    library(data.table)
    files <- list.files(pattern = "\\\\.psi\$")
    files <- sort(files) 
    
    if (length(files) == 0) {
        stop("No PSI files found to merge.")
    }
    
    message(paste("Found", length(files), "files to merge."))
    

    message(paste("Processing base file:", files[1]))
    merged_dt <- fread(files[1])
    

    colnames(merged_dt)[1] <- "ASE"
    

    if (ncol(merged_dt) > 1) {
        for (j in seq.int(2, ncol(merged_dt))) {
            set(merged_dt, which(is.na(merged_dt[[j]])), j, 0)
        }
    }
    
   
    if (length(files) > 1) {
        for (i in 2:length(files)) {
            f <- files[i]
            message(paste("Merging:", f))
            
            dt <- fread(f)
            colnames(dt)[1] <- "ASE"
            

            if (nrow(dt) != nrow(merged_dt)) {
                stop(paste("Error: Row count mismatch in file:", f))
            }

            if (!identical(dt[[1]], merged_dt[["ASE"]])) {
                stop(paste("Error: ASE row order mismatch in file:", f))
            }
            

            if (ncol(dt) > 1) {
                for (j in seq.int(2, ncol(dt))) {
                    set(dt, which(is.na(dt[[j]])), j, 0)
                }
            }
            

            dt[, ASE := NULL]
            

            merged_dt <- cbind(merged_dt, dt)
        }
    }
    
    message("Writing merged output...")
    fwrite(merged_dt, "merged_project_events.psi", sep = "\\t", quote = FALSE, row.names = FALSE)
    message("Done.")
    RSCRIPT
    """
}

process ASE_ANALYSIS {
  label 'sclong'
  tag "ase_analysis"
  
  publishDir "${params.outdir ?: 'Result'}/AlternativeSplicing/ASE_Analysis", mode: 'copy'
  
  input:
    path merged_psi
    path merged_events
    path seurat_rds
    path isoquant_dir
    path pipeline_bin
  
  output:
    path "ase_analysis_output/*", emit: analysis_results
    path "ase_analysis_summary.txt", emit: summary
  
  script:
    output_dir = "ase_analysis_output"
    """
    # Create output directory
    mkdir -p ${output_dir}
    
    # Run ASE analysis with optparse
    NOVEL_VS_KNOWN="${isoquant_dir}/OUT/OUT.novel_vs_known.SQANTI-like.tsv"
    if [[ ! -f "\${NOVEL_VS_KNOWN}" ]]; then
      echo "ERROR: Novel-vs-known annotation not found: \${NOVEL_VS_KNOWN}" >&2
      exit 1
    fi

    RSCRIPT_PATH="${pipeline_bin}/AlternativeSplicing/AsEventAnalysis.R"
    
    "${params.rscript_bin ?: 'Rscript'}" "\${RSCRIPT_PATH}" \
      --psi_file "${merged_psi}" \
      --ioe_file "${merged_events}" \
      --seurat_rds "${seurat_rds}" \
      --novel_vs_known "\${NOVEL_VS_KNOWN}" \
      --output_dir "${output_dir}" \
      --min_cell_prop ${params.ase_min_cell_prop ?: 0.01}
    
    # Create summary file
    echo "ASE Analysis Summary" > ase_analysis_summary.txt
    echo "====================" >> ase_analysis_summary.txt
    echo "Timestamp: \$(date)" >> ase_analysis_summary.txt
    echo "Input files:" >> ase_analysis_summary.txt
    echo "- PSI file: ${merged_psi}" >> ase_analysis_summary.txt
    echo "- Events file: ${merged_events}" >> ase_analysis_summary.txt
    echo "- Seurat object: ${seurat_rds}" >> ase_analysis_summary.txt
    echo "- Novel vs Known: \${NOVEL_VS_KNOWN}" >> ase_analysis_summary.txt
    echo "" >> ase_analysis_summary.txt
    echo "Output generated in: ${output_dir}" >> ase_analysis_summary.txt
    echo "Number of plots generated: \$(find ${output_dir} -name '*.pdf' | wc -l)" >> ase_analysis_summary.txt
    """
}

process DOMINANT_ISOFORM_DETECTION {
  label 'sclong'
  tag "dominant_isoform"
  
  publishDir "${params.outdir ?: 'Result'}/AlternativeSplicing/DominantIsoform", mode: 'copy'
  
  input:
    path seurat_rds
    path transcript_gtf
    path pipeline_bin
  
  output:
    path "dominant_isoform_output", emit: analysis_results_dir
    path "dominant_summary.txt", emit: summary
  
  script:
    output_dir = "dominant_isoform_output"
    
    // Get absolute path for the GTF file
    def gtf_path = transcript_gtf.toString()
    
    // Create stage mapping file
    stage_mapping_content = ""
    params.sample_stage_mapping.each { mapping ->
      stage_mapping_content += "${mapping[0]}\t${mapping[1]}\n"
    }
    stage_mapping_file = "sample_stage_mapping.tsv"
    
    // Create comparisons file
    comparisons_content = ""
    (params.comparisons).each { comparison ->
      comparisons_content += "${comparison[0]}\t${comparison[1]}\t${comparison[2]}\n"
    }
    comparisons_file = "comparisons_config.tsv"
    
    """
    # Create output directory
    mkdir -p ${output_dir}
    
    # Create stage mapping file
    echo -e "sample\\tstage" > ${stage_mapping_file}
    cat << 'EOF' >> ${stage_mapping_file}
${stage_mapping_content}
EOF
    
    # Create comparisons file
    echo -e "group1\\tgroup2\\tname" > ${comparisons_file}
    cat << 'EOF' >> ${comparisons_file}
${comparisons_content}
EOF
    
    # Run dominant isoform analysis
    RSCRIPT_PATH="${pipeline_bin}/AlternativeSplicing/dominant_isoform_detection.R"
    
    # Pass the absolute path to R script
    "${params.rscript_bin ?: 'Rscript'}" "\${RSCRIPT_PATH}" \\
      --seurat_rds "${seurat_rds}" \\
      --gtf_file "${gtf_path}" \\
      --stage_mapping "${stage_mapping_file}" \\
      --comparisons "${comparisons_file}" \\
      --output_dir "${output_dir}" \\
      --min_cells ${params.dominant_min_cells ?: 10} \\
      --n_cores ${params.dominant_n_cores ?: 15} \\
      --logfc_threshold ${params.dominant_logfc_threshold ?: 1.0}
    
    # Create summary file
    echo "Dominant Isoform Analysis Summary" > dominant_summary.txt
    echo "==================================" >> dominant_summary.txt
    echo "Timestamp: \$(date)" >> dominant_summary.txt
    echo "Input files:" >> dominant_summary.txt
    echo "- Seurat RDS: ${seurat_rds}" >> dominant_summary.txt
    echo "- Transcript GTF: ${gtf_path}" >> dominant_summary.txt
    echo "" >> dominant_summary.txt
    echo "Parameters:" >> dominant_summary.txt
    echo "- Min cells: ${params.dominant_min_cells ?: 10}" >> dominant_summary.txt
    echo "- Number of cores: ${params.dominant_n_cores ?: 15}" >> dominant_summary.txt
    echo "- LogFC threshold: ${params.dominant_logfc_threshold ?: 1.0}" >> dominant_summary.txt
    echo "" >> dominant_summary.txt
    echo "Output generated in: ${output_dir}" >> dominant_summary.txt
    """
}

process ISOFORM_SWITCH_PREPARATION {
  label 'sclong'
  tag "isoform_switch_prep"
  
  publishDir "${params.outdir ?: 'Result'}/AlternativeSplicing/IsoformSwitch/preparation", mode: 'copy'
  
  input:
    path dominant_results_dir
    path seurat_tr
    path seurat_gene
    path transcript_gtf
    path transcript_fasta
    path pipeline_bin
  
  output:
    path "isoform_switch_prep_output/*", emit: preparation_results
    path "switchList_prep.rds", optional: true, emit: switchlist_rds
    path "preparation_summary.txt", emit: prep_summary
  
  script:
    output_dir = "isoform_switch_prep_output"
    
    // Get absolute paths
    def dominant_dir = dominant_results_dir.toString()
    def gtf_path = transcript_gtf.toString()
    def fasta_path = transcript_fasta.toString()
    def seurat_tr_path = seurat_tr.toString()
    def seurat_gene_path = seurat_gene.toString()
    
    // Create stage mapping file from params
    stage_mapping_content = ""
    params.sample_stage_mapping.each { mapping ->
      stage_mapping_content += "${mapping[0]}\t${mapping[1]}\n"
    }
    stage_mapping_file = "sample_stage_mapping.tsv"
    
    // Get comparison from config
    comparisons = params.comparisons 
    comparison = comparisons[0]  // Use first comparison
    comparison_group1 = comparison[0]
    comparison_group2 = comparison[1]
    comparison_name = comparison[2]
    
    // Get external tool result files (empty string if not provided)
    cpc2_result = params.cpc2_result_file ?: ""
    pfam_result = params.pfam_result_file ?: ""
    signalp_result = params.signalp_result_file ?: ""
    iupred2a_result = params.iupred2a_result_file ?: ""
    
    """
    # Create output directory
    mkdir -p ${output_dir}
    
    # Create stage mapping file
    echo -e "sample\\tstage" > ${stage_mapping_file}
    cat << 'EOF' >> ${stage_mapping_file}
${stage_mapping_content}
EOF
    
    # Run isoform switch preparation
    RSCRIPT_PATH="${pipeline_bin}/AlternativeSplicing/isoform_switch_preparation.R"
    
    "${params.rscript_bin ?: 'Rscript'}" "\${RSCRIPT_PATH}" \\
      --dominant_results "${dominant_dir}" \\
      --seurat_rds "${seurat_tr_path}" \\
      --seurat_gene_rds "${seurat_gene_path}" \\
      --gtf_file "${gtf_path}" \\
      --fasta_file "${fasta_path}" \\
      --output_dir "${output_dir}" \\
      --stage_mapping "${stage_mapping_file}" \\
      --comparison "${comparison_name}" \\
      --reference_group "${comparison_group1}" \\
      --test_group "${comparison_group2}" \\
      --cpc2_result "${cpc2_result}" \\
      --pfam_result "${pfam_result}" \\
      --signalp_result "${signalp_result}" \\
      --iupred2a_result "${iupred2a_result}" \\
      --n_cores ${params.isoform_n_cores ?: 15} \\
      --dIF_cutoff ${params.isoform_dIF_cutoff ?: 0.1} \\
      --alpha ${params.isoform_alpha ?: 0.05}
    
    # Find and rename the switchlist file
    find ${output_dir} -name "switchList_prep_*.rds" -exec cp {} switchList_prep.rds \\;
    
    # Create summary file
    echo "Isoform Switch Preparation Summary" > preparation_summary.txt
    echo "===================================" >> preparation_summary.txt
    echo "Timestamp: \$(date)" >> preparation_summary.txt
    echo "Input files:" >> preparation_summary.txt
    echo "- Dominant results directory: ${dominant_dir}" >> preparation_summary.txt
    echo "- Seurat transcript: ${seurat_tr_path}" >> preparation_summary.txt
    echo "- Seurat gene: ${seurat_gene_path}" >> preparation_summary.txt
    echo "- Transcript GTF: ${gtf_path}" >> preparation_summary.txt
    echo "- Transcript FASTA: ${fasta_path}" >> preparation_summary.txt
    echo "" >> preparation_summary.txt
    echo "Parameters:" >> preparation_summary.txt
    echo "- Comparison: ${comparison_name}" >> preparation_summary.txt
    echo "- Groups: ${comparison_group1} vs ${comparison_group2}" >> preparation_summary.txt
    echo "- Number of cores: ${params.isoform_n_cores ?: 15}" >> preparation_summary.txt
    echo "- dIF cutoff: ${params.isoform_dIF_cutoff ?: 0.1}" >> preparation_summary.txt
    echo "- Alpha: ${params.isoform_alpha ?: 0.05}" >> preparation_summary.txt
    echo "" >> preparation_summary.txt
    echo "Output generated in: ${output_dir}" >> preparation_summary.txt
    """
}

process ISOFORM_SWITCH_CONSEQUENCES {
  label 'sclong'
  tag "isoform_switch_conseq"
  
  publishDir "${params.outdir ?: 'Result'}/AlternativeSplicing/IsoformSwitch/consequences", mode: 'copy'
  
  input:
    path switchlist_rds
    path pipeline_bin
  
  output:
    path "isoform_switch_conseq_output/", emit: consequences_results
    path "consequences_summary.txt", emit: conseq_summary
  
  script:
    output_dir = "isoform_switch_conseq_output"
    
    // Get absolute path
    def switchlist_path = switchlist_rds.toString()
    
    """
    # Create output directory
    mkdir -p ${output_dir}
    
    # Run consequences analysis
    RSCRIPT_PATH="${pipeline_bin}/AlternativeSplicing/isoform_switch_consequences.R"
    
    "${params.rscript_bin ?: 'Rscript'}" "\${RSCRIPT_PATH}" \\
      --switchlist_file "${switchlist_path}" \\
      --output_dir "${output_dir}" \\
      --dIF_cutoff ${params.isoform_dIF_cutoff ?: 0.1} \\
      --alpha ${params.isoform_alpha ?: 0.05}
    
    # Create summary file
    echo "Isoform Switch Consequences Analysis Summary" > consequences_summary.txt
    echo "=============================================" >> consequences_summary.txt
    echo "Timestamp: \$(date)" >> consequences_summary.txt
    echo "Input switch list: ${switchlist_path}" >> consequences_summary.txt
    echo "" >> consequences_summary.txt
    echo "Parameters:" >> consequences_summary.txt
    echo "- dIF cutoff: ${params.isoform_dIF_cutoff ?: 0.1}" >> consequences_summary.txt
    echo "- Alpha: ${params.isoform_alpha ?: 0.05}" >> consequences_summary.txt
    echo "" >> consequences_summary.txt
    echo "Output generated in: ${output_dir}" >> consequences_summary.txt
    """
}


process LLM_ISOFORM_SWITCH_INTERPRETATION {
  label 'sclong'
  tag "llm_isoform_switch_interpretation"
  
  publishDir "${params.outdir ?: 'Result'}/AlternativeSplicing/IsoformSwitch/interpretation", mode: 'copy'
  
  input:
    path consequences_results_dir
    path pipeline_bin
  
  output:
    path "biological_interpretation.txt", emit: interpretation_report
    // Historical filename retained for backward compatibility.
    path "deepseek_api_response.json", emit: api_response
    path "gene_summary.txt", emit: gene_summary

  when:
    params.llm_enabled
  
  script:
    // Find the switch_consequences.txt file in the input directory
    def consequences_file = "${consequences_results_dir}/switch_consequences.txt"
    
    // Set study description
    def study_description = params.study_description 
    
    """
    # Check if consequences file exists
    if [[ ! -f "${consequences_file}" ]]; then
      echo "ERROR: Consequences file not found: ${consequences_file}" >&2
      exit 1
    fi
    
    # Run the provider-agnostic interpretation script.
    PYTHON_SCRIPT="${pipeline_bin}/llm/isoform_switch_interpretation.py"
    
    python3 \${PYTHON_SCRIPT} \
      --consequences_file "${consequences_file}" \
      --output_report "biological_interpretation.txt" \
      --output_api "deepseek_api_response.json" \
      --study_description "${study_description}" \
      --base_url "${params.llm_base_url}" \
      --model "${params.llm_model}" \
      --save_prompt
    
    # Check if report was generated
    if [[ ! -s "biological_interpretation.txt" ]]; then
      echo "WARNING: Biological interpretation report is empty or skipped" >&2
    else
      echo "Biological interpretation generated successfully"
      echo "Report preview:"
      head -10 biological_interpretation.txt
    fi
    """
}

workflow AlternativeSplicing {
  take:
    seurat_tr
    seurat_gene
    gtf_file
    transcript_fasta
    isoquant_out
  
  main:
    pipeline_bin = file("${projectDir}/bin", checkIfExists: true)

    // ------------------------------------------------------------------
    // Part 1: Alternative Splicing Event Detection
    // ------------------------------------------------------------------
    split_results = SPLIT_EXPRESSION(seurat_tr, pipeline_bin)
    events_results = GENERATE_AS_EVENTS(gtf_file)
    split_files_channel = split_results.split_files.flatten()
    psi_results = CALCULATE_PSI(events_results.merged_events, split_files_channel)
    all_psi_files = psi_results.psi_file.collect()
    merge_results = MERGE_AS_RESULTS(all_psi_files)
    
    // ------------------------------------------------------------------
    // Part 2: ASE Analysis
    // ------------------------------------------------------------------
    ase_analysis = ASE_ANALYSIS(
      merge_results.merged_psi,
      events_results.merged_events,
      seurat_tr,
      isoquant_out,
      pipeline_bin
    )
    
    // ------------------------------------------------------------------
    // Part 3: Dominant Isoform Detection
    // ------------------------------------------------------------------
    
    dominant_isoform = DOMINANT_ISOFORM_DETECTION(
      seurat_tr,
      gtf_file,
      pipeline_bin
    )
    
    // ------------------------------------------------------------------
    // Part 4: Isoform Switch Analysis - Split into two parts
    // ------------------------------------------------------------------
    def dominant_results_dir = dominant_isoform.analysis_results_dir
    
    // 4a: Preparation phase
    isoform_switch_prep = ISOFORM_SWITCH_PREPARATION(
      dominant_results_dir,
      seurat_tr,
      seurat_gene,
      gtf_file,
      transcript_fasta,
      pipeline_bin
    )
    
    // 4b: Consequences phase
    isoform_switch_conseq = ISOFORM_SWITCH_CONSEQUENCES(
      isoform_switch_prep.switchlist_rds,
      pipeline_bin
    )

    // 4c: Optional provider-agnostic LLM interpretation.
    biological_interpretation = LLM_ISOFORM_SWITCH_INTERPRETATION(
      isoform_switch_conseq.consequences_results,
      pipeline_bin
    )
  
  emit:
    // ASE-related outputs
    merged_psi = merge_results.merged_psi
    merged_events = events_results.merged_events
    ase_analysis_results = ase_analysis.analysis_results
    ase_summary = ase_analysis.summary
    
    // Dominant isoform-related outputs
    dominant_results = dominant_isoform.analysis_results_dir
    dominant_summary = dominant_isoform.summary
    
    // Isoform switch preparation outputs
    switch_prep_results = isoform_switch_prep.preparation_results
    switchlist_file = isoform_switch_prep.switchlist_rds
    switch_prep_summary = isoform_switch_prep.prep_summary
    
    // Isoform switch consequences outputs
    switch_conseq_results = isoform_switch_conseq.consequences_results
    switch_conseq_summary = isoform_switch_conseq.conseq_summary

    // Biological interpretation outputs (may be empty if API key not provided)
    interpretation_report = biological_interpretation.interpretation_report
    api_response = biological_interpretation.api_response
    gene_summary = biological_interpretation.gene_summary

}
