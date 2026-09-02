// main.nf
nextflow.enable.dsl=2

include { InputAndPreprocess } from './modules/InputAndPreprocess.nf'
include { AlignmentAndTranscriptReconstruction } from './modules/AlignmentAndTranscriptReconstruction.nf'
include { PreliminarySingleCellAnalysis } from './modules/PreliminarySingleCellAnalysis.nf'
include { AlternativeSplicing } from './modules/AlternativeSplicing.nf'
include { IntegrativeAnalysis } from './modules/IntegrativeAnalysis.nf'

workflow {
  def llm_key_present = !(System.getenv('LLM_API_KEY') ?: '').trim().isEmpty()
  def llm_base_url = (params.llm_base_url ?: '').toString().trim()
  if (params.llm_enabled && !llm_key_present) {
    error "LLM is enabled but LLM_API_KEY is not set"
  }
  if (params.llm_enabled && !llm_base_url) {
    error "LLM is enabled but LLM_BASE_URL/params.llm_base_url is empty"
  }
  if (params.llm_enabled && !(llm_base_url ==~ /(?i)https?:\/\/.+/)) {
    error "LLM_BASE_URL/params.llm_base_url must be an absolute http(s) URL"
  }
  if (!(params.comparisons instanceof List) || params.comparisons.isEmpty()) {
    error "params.comparisons must contain at least one [group1, group2, name] entry"
  }
  params.comparisons.each { comparison ->
    if (!(comparison instanceof List) || comparison.size() != 3) {
      error "Each params.comparisons entry must be [group1, group2, name]: ${comparison}"
    }
    if (comparison.any { value -> value == null || !value.toString().trim() }) {
      error "Comparison values must be non-empty: ${comparison}"
    }
    if (comparison.any { value -> value.toString() != value.toString().trim() }) {
      error "Comparison values cannot have leading or trailing whitespace: ${comparison}"
    }
    if (comparison[0].toString().trim() == comparison[1].toString().trim()) {
      error "Comparison groups must be different: ${comparison}"
    }
    if (comparison[0..1].any { value -> value.toString().contains('\t') || value.toString().contains('\n') || value.toString().contains('\r') }) {
      error "Comparison group names cannot contain tabs or newlines: ${comparison}"
    }
    if (!(comparison[2].toString() ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
      error "Comparison names must contain only letters, numbers, '.', '_' or '-': ${comparison[2]}"
    }
  }
  def comparison_names = params.comparisons.collect { comparison -> comparison[2].toString() }
  if (comparison_names.toSet().size() != comparison_names.size()) {
    error "params.comparisons contains duplicate comparison names: ${comparison_names}"
  }
  if (!(params.sample_stage_mapping instanceof List) || params.sample_stage_mapping.isEmpty()) {
    error "params.sample_stage_mapping must contain at least one [sample, stage] entry"
  }
  params.sample_stage_mapping.each { mapping ->
    if (!(mapping instanceof List) || mapping.size() != 2 ||
        mapping.any { value -> value == null || !value.toString().trim() }) {
      error "Each sample-stage entry must contain two non-empty values: ${mapping}"
    }
    if (mapping.any { value -> value.toString().contains('\t') || value.toString().contains('\n') || value.toString().contains('\r') }) {
      error "Sample-stage values cannot contain tabs or newlines: ${mapping}"
    }
    if (mapping.any { value -> value.toString() != value.toString().trim() }) {
      error "Sample-stage values cannot have leading or trailing whitespace: ${mapping}"
    }
    if (!(mapping[0].toString() ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
      error "Mapped sample IDs must contain only letters, numbers, '.', '_' or '-': ${mapping[0]}"
    }
  }
  def mapped_samples = params.sample_stage_mapping.collect { mapping -> mapping[0].toString().trim() }
  if (mapped_samples.toSet().size() != mapped_samples.size()) {
    error "params.sample_stage_mapping contains duplicate samples: ${mapped_samples}"
  }
  def mapped_stages = params.sample_stage_mapping
    .collect { mapping -> mapping[1].toString().trim() }
    .toSet()
  params.comparisons.each { comparison ->
    def missing_groups = comparison[0..1]
      .collect { value -> value.toString().trim() }
      .findAll { group -> !mapped_stages.contains(group) }
    if (missing_groups) {
      error "Comparison ${comparison[2]} uses groups absent from sample_stage_mapping: ${missing_groups.unique()}"
    }
  }

  // Module (1): InputAndPreprocess
  pre = InputAndPreprocess(params.metadata, params.whitelist)

  // Module (2): AlignmentAndTranscriptReconstruction
  align = AlignmentAndTranscriptReconstruction(
    pre.filtered_fastq,
    params.reference_fa,
    params.genedb_gtf
  )

  // Module (3): PreliminarySingleCellAnalysis
  sc = PreliminarySingleCellAnalysis(align.isoquant_out)

  // Module (4): AlternativeSplicing
  splicing = AlternativeSplicing(
    sc.seurat_tr,
    sc.seurat_gene,
    align.models_gtf,         
    align.models_fa_1line,    
    align.isoquant_out        
  )
  
  // Module (5): IntegrativeAnalysis 
  IntegrativeAnalysis(
    sc.seurat_tr,           // seurat_tr
    sc.seurat_gene,         // seurat_gene
    align.models_gtf,       // gtf_file
    splicing.dominant_results, // dominant_results
    align.orf_fa            // orf_fasta
  )

}
