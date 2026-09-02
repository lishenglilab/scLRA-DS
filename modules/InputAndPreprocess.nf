nextflow.enable.dsl=2

process PARSE_METADATA {
  label 'preprocess'
  tag "parse_metadata"

  input:
    path metadata
    path pipeline_bin

  output:
    path "samples.tsv", emit: samples

  script:
    """
    python "${pipeline_bin}/parse_metadata.py" --metadata "${metadata}" --out samples.tsv
    """
}

process BLAZE_RUN {
  label 'preprocess'
  tag { sample_id }

  publishDir "${params.outdir ?: 'Result'}/InputAndPreprocess/blaze/${sample_id}", mode: 'copy'

  input:
    tuple val(sample_id), path(fastq), path(whitelist), path(metadata)

  output:
    tuple val(sample_id),
          path("${sample_id}.blaze*"),
          path("${sample_id}.meta.tsv"),
          emit: blaze_out
    tuple val(sample_id), path("${sample_id}.blaze.matched_reads.fastq.gz"), emit: matched_reads

  script:
    """
    if [[ -z "${params.expect_cells}" ]]; then
      echo "ERROR: params.expect_cells is not set" >&2
      exit 1
    fi
    if [[ -z "${params.threads}" ]]; then
      echo "ERROR: params.threads is not set" >&2
      exit 1
    fi

    cp "${metadata}" "${sample_id}.meta.tsv"

    blaze \
      --expect-cells ${params.expect_cells} \
      --output-prefix "${sample_id}.blaze." \
      --threads ${params.threads} \
      --high-sensitivity-mode \
      --minimal_stdout \
      --full-bc-whitelist "${whitelist}" \
      "${fastq}"
    """
}

process FILTER_READS {
  // NanoFilt is provided by envs/preprocess.yml, not the main sclong env.
  label 'preprocess'
  tag { sample_id }

  publishDir "${params.outdir ?: 'Result'}/InputAndPreprocess/filtered_reads/${sample_id}", mode: 'copy'

  input:
    tuple val(sample_id), path(fastq_gz)

  output:
    tuple val(sample_id), path("${sample_id}.hq.fastq.gz"), emit: hq_fastq

  script:
    """
    gunzip -c "${fastq_gz}" \
      | NanoFilt -q ${params.nanofilt_q ?: 10} \
      | gzip > "${sample_id}.hq.fastq.gz"
    """
}

process NANOPLOT_QC {
  label 'preprocess'
  tag "NANOPLOT_QC"
  
  publishDir "${params.outdir ?: 'Result'}/InputAndPreprocess/nanoplot_qc", mode: 'copy'

  input:
    path fastq_files

  output:
    path("all_samples_NanoStats.txt"), emit: nanoplot_report
    path("all_samples_NanoPlot_output/*"), emit: nanoplot_output
  
  script:
    def fastq_args = fastq_files.collect { fastq -> "\"${fastq}\"" }.join(' ')
    """
    mkdir -p all_samples_NanoPlot_output
    
    NanoPlot \
      -t ${params.threads ?: 2} \
      --fastq ${fastq_args} \
      -o all_samples_NanoPlot_output \
      --title "All Samples Combined NanoPlot Report"
    
    cp all_samples_NanoPlot_output/NanoStats.txt all_samples_NanoStats.txt
    """
}

process LLM_READ_QC {
  label 'preprocess'
  tag "llm_read_qc"
  
  // Keep the historical directory name so existing result consumers do not break.
  publishDir "${params.outdir ?: 'Result'}/InputAndPreprocess/deepseek_qc", mode: 'copy'
  
  input:
    path nanoplot_stats_file
    path pipeline_bin
  
  output:
    path "qc_report.txt", emit: llm_qc_report

  when:
    params.llm_enabled
  
  script:
    """
    # Run analysis
    python "${pipeline_bin}/llm/read_qc.py" \
      --file_path "${nanoplot_stats_file}" \
      --output_text_file "qc_report.txt" \
      --species "${params.species ?: 'Unknown'}" \
      --tissue "${params.tissue ?: 'Unknown'}" \
      --base_url "${params.llm_base_url}" \
      --model "${params.llm_model}"
    """
}

workflow InputAndPreprocess {
  take:
    metadata_path
    whitelist_path

  main:
    meta_path = file(metadata_path)
    wl_path   = file(whitelist_path)
    pipeline_bin = file("${projectDir}/bin", checkIfExists: true)

    parsed = PARSE_METADATA(meta_path, pipeline_bin)

    // 准备原始数据通道给BLAZE_RUN
    ch_original_samples = parsed.samples
      .splitCsv(header:true, sep:"\t")
      .map { row ->
        tuple(row.sample_id, file(row.fastq), wl_path, meta_path)
      }

    
    blaze_results = BLAZE_RUN(ch_original_samples)


    filtered = FILTER_READS(blaze_results.matched_reads)
    
    ch_all_fastq_files = filtered.hq_fastq
      .map { _sample_id, fastq -> fastq }
      .collect()
    nanoplot_qc = NANOPLOT_QC(ch_all_fastq_files)
    
    // Optional provider-agnostic LLM report.
    llm_qc = LLM_READ_QC(nanoplot_qc.nanoplot_report, pipeline_bin)

  emit:
    blaze_out = blaze_results.blaze_out
    matched_reads = blaze_results.matched_reads
    filtered_fastq = filtered.hq_fastq
    qc_reports = nanoplot_qc.nanoplot_report
    qc_outputs = nanoplot_qc.nanoplot_output
    llm_qc_report = llm_qc.llm_qc_report
    // Backward-compatible alias for existing consumers.
    deepseek_qc_report = llm_qc.llm_qc_report
}
