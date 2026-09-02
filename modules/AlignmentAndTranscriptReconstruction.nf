nextflow.enable.dsl=2

/*
 * This process creates read groups file for IsoQuant
 */
process MAKE_READ_GROUPS {
  label 'sclong'
  tag "read_groups"

  publishDir "${params.outdir ?: 'Result'}/AlignmentAndTranscriptReconstruction", mode: 'copy'

  input:
    val sample_ids
    path fastqs
    path pipeline_bin

  output:
    path "all_read_groups.tsv", emit: read_groups

  script:
    def fq_args  = fastqs.collect { fq -> "\"${fq}\"" }.join(' ')
    def sid_args = sample_ids.collect { sid -> "\"${sid}\"" }.join(' ')
    """
    python "${pipeline_bin}/make_read_groups.py" \
      --fastq ${fq_args} \
      --sample-id ${sid_args} \
      --out all_read_groups.tsv
    """
}

/*
 * IsoQuant for single-cell long-read data
 */
process ISOQUANT_SC {
  label 'sclong'
  tag "isoquant_sc"

  publishDir "${params.outdir ?: 'Result'}/AlignmentAndTranscriptReconstruction/isoquant", mode: 'copy'

  input:
    path reference_fa
    path genedb_gtf
    path read_groups_tsv
    path fastqs

  output:
    path "IsoQuant_sc", emit: isoquant_dir

  script:
    def fastq_args = fastqs.collect { fq -> "\"${fq}\"" }.join(' ')
    """
    mkdir -p IsoQuant_sc

    isoquant.py \
      --reference "${reference_fa}" \
      --genedb "${genedb_gtf}" \
      --complete_genedb \
      --fastq ${fastq_args} \
      --data_type nanopore \
      -o IsoQuant_sc \
      --threads ${params.threads ?: 10} \
      --sqanti_output \
      --read_group "file:${read_groups_tsv}"
    """
}

/*
 * Extract transcript sequences from IsoQuant output
 */
process EXTRACT_TRANSCRIPTS {
  label 'sclong'
  tag "extract_transcripts"

  publishDir "${params.outdir ?: 'Result'}/AlignmentAndTranscriptReconstruction/orf", mode: 'copy'

  input:
    path reference_fa
    path isoquant_dir

  output:
    path "OUT.transcript_models.gtf",      emit: models_gtf
    path "OUT.transcript_models.fa",       emit: models_fa
    path "OUT.transcript_models.1line.fa", emit: models_fa_1line

  script:
    """
    GTF="${isoquant_dir}/OUT/OUT.transcript_models.gtf"

    if [[ ! -s "\$GTF" ]]; then
      echo "ERROR: missing or empty: \$GTF" >&2
      exit 1
    fi

    cp "\$GTF" OUT.transcript_models.gtf

    gffread OUT.transcript_models.gtf \
      -g "${reference_fa}" \
      -w OUT.transcript_models.fa

    if [[ ! -s OUT.transcript_models.fa ]]; then
      echo "ERROR: gffread output is missing or empty: OUT.transcript_models.fa" >&2
      exit 1
    fi

    python - << 'PY'
    in_fa = "OUT.transcript_models.fa"
    out_fa = "OUT.transcript_models.1line.fa"

    with open(in_fa, "r") as fin, open(out_fa, "w") as fout:
      header = None
      seq_parts = []
      for line in fin:
        line = line.rstrip("\\n")
        if not line:
          continue
        if line.startswith(">"):
          if header is not None:
            fout.write(header + "\\n")
            fout.write("".join(seq_parts) + "\\n")
          header = line
          seq_parts = []
        else:
          seq_parts.append(line.strip())
      if header is not None:
        fout.write(header + "\\n")
        fout.write("".join(seq_parts) + "\\n")
    PY

    test -s OUT.transcript_models.1line.fa
    """
}

/*
 * ORF calling from transcript sequences
 */
process ORF_CALL {
  label 'sclong'
  tag "orf_call"

  publishDir "${params.outdir ?: 'Result'}/AlignmentAndTranscriptReconstruction/orf", mode: 'copy'

  input:
    path transcript_fa_1line
    path pipeline_bin

  output:
    path "ORF.fa", emit: orf_fa

  script:
    """
    python "${pipeline_bin}/ORF_generate.py" \
      -i "${transcript_fa_1line}" \
      -o . \
      -p ORF

    test -s ORF.fa
    """
}

/*
 * Optional transcriptome assembly quality assessment using an LLM API
 */
process LLM_ASSEMBLY_QC {
  label 'preprocess'
  tag "llm_assembly_qc"
  
  publishDir "${params.outdir ?: 'Result'}/AlignmentAndTranscriptReconstruction/assembly_qc", mode: 'copy'
  
  input:
    path isoquant_dir
    path pipeline_bin
  
  output:
    path "assembly_qc_report.txt", emit: assembly_qc_report
    path "assembly_qc_report_metrics.json", emit: assembly_qc_metrics

  when:
    params.llm_enabled
  
  script:
    """
    # Debug: show directory contents
    echo "IsoQuant directory: ${isoquant_dir}"
    echo "Directory contents:"
    ls -la "${isoquant_dir}/OUT/" 2>/dev/null || echo "OUT directory not found"
    
    # Run the Python script
    python "${pipeline_bin}/llm/assembly_qc.py" \
      --isoquant_dir "${isoquant_dir}" \
      --output_text_file "assembly_qc_report.txt" \
      --species "${params.species ?: 'Unknown'}" \
      --tissue "${params.tissue ?: 'Unknown'}" \
      --base_url "${params.llm_base_url}" \
      --model "${params.llm_model}" \
      --save_metrics
    
    # Check if report was generated
    if [[ ! -s "assembly_qc_report.txt" ]]; then
      echo "WARNING: Assembly QC report is empty" >&2
    else
      echo "Assembly QC report generated successfully"
    fi
    """
}

/*
 * Main workflow for transcript reconstruction
 * Takes filtered FASTQ files from InputAndPreprocess module
 */
workflow AlignmentAndTranscriptReconstruction {
  take:
    filtered_fastq  // Expects tuple(sample_id, fastq) from InputAndPreprocess
    reference_fa
    genedb_gtf

  main:
    pipeline_bin = file("${projectDir}/bin", checkIfExists: true)

    // Collect sample IDs and fastq files separately for MAKE_READ_GROUPS
    sample_ids = filtered_fastq.map { sid, _fq -> sid }.collect()
    hq_fastqs  = filtered_fastq.map { _sid, fq -> fq }.collect()

    // Create read groups file
    rg = MAKE_READ_GROUPS(sample_ids, hq_fastqs, pipeline_bin)
    
    // Run IsoQuant
    iso = ISOQUANT_SC(reference_fa, genedb_gtf, rg.read_groups, hq_fastqs)

    // Extract transcripts
    tx = EXTRACT_TRANSCRIPTS(reference_fa, iso.isoquant_dir)
    
    // Call ORFs
    orf = ORF_CALL(tx.models_fa_1line, pipeline_bin)
    
    // Assembly quality assessment (optional, requires API key)
    assembly_qc = LLM_ASSEMBLY_QC(iso.isoquant_dir, pipeline_bin)

  emit:
    filtered_fastq  = filtered_fastq      // Pass through for downstream use
    read_groups     = rg.read_groups
    isoquant_out    = iso.isoquant_dir
    models_gtf      = tx.models_gtf
    models_fa       = tx.models_fa
    models_fa_1line = tx.models_fa_1line
    orf_fa          = orf.orf_fa
    assembly_qc_report = assembly_qc.assembly_qc_report
    assembly_qc_metrics = assembly_qc.assembly_qc_metrics
}
