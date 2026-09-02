# scLRA-LLM

scLRA-LLM is a Nextflow DSL2 pipeline for single-cell long-read RNA-seq analysis. It connects read preprocessing, IsoQuant transcript reconstruction, gene- and transcript-level single-cell analysis, alternative-splicing analysis, ORF analysis, and optional LLM-assisted reporting in one workflow.

LLM steps use an OpenAI-compatible Chat Completions endpoint and are not tied to a specific provider. All biological analysis steps can run with LLM support disabled.

## Workflow

| Module | Main tasks |
| --- | --- |
| `InputAndPreprocess` | Metadata validation, BLAZE barcode processing, NanoFilt filtering, NanoPlot QC, optional LLM read-QC report |
| `AlignmentAndTranscriptReconstruction` | IsoQuant transcript reconstruction and quantification, transcript sequence extraction, ORF prediction, optional LLM assembly-QC report |
| `PreliminarySingleCellAnalysis` | Gene/transcript matrices, Seurat QC and clustering, marker detection, optional LLM cell-type annotation and QC report |
| `AlternativeSplicing` | SUPPA2 event quantification, dominant-isoform analysis, isoform-switch preparation and consequences, optional LLM interpretation |
| `IntegrativeAnalysis` | DDTU, DTE, DGE and DGU integration, ORF clustering, optional LLM interpretation |

The five modules are orchestrated from [`main.nf`](main.nf), with process definitions under [`modules/`](modules/).

## Requirements

- Linux
- Nextflow with DSL2 support; the current workflow was validated with Nextflow 25.10.2
- A Java version supported by the installed Nextflow release (Java 17 or newer for current Nextflow releases)
- Conda or Mamba for the process environments
- Sufficient CPU, memory and storage for IsoQuant and single-cell analysis

The pipeline provides Conda YAML files for its Python and command-line dependencies. The main environment currently does not fully lock the R/Bioconductor stack; see [R environment](#r-environment).

## Installation

```bash
git clone https://github.com/lishenglilab/scLRA-LLM.git
cd scLRA-LLM
unzip assets/3M-february-2018.zip -d assets
```

The archive contains the 10x 3M barcode whitelist referenced by the default configuration. If a different chemistry or whitelist is used, update `params.whitelist` instead. No API credentials are required when LLM support is disabled.

## Input configuration

Before running, edit [`nextflow.config`](nextflow.config) and check at least:

- `params.metadata`
- `params.whitelist`
- `params.reference_fa`
- `params.genedb_gtf`
- `params.sample_stage_mapping`
- `params.comparisons`
- Thread counts and analysis thresholds appropriate for the server

### Metadata

Metadata may be TSV or CSV and must contain `sample_id` and `fastq` columns. Absolute FASTQ paths are recommended.

```text
sample_id<TAB>fastq
AD1<TAB>/absolute/path/AD1.fastq
AD2<TAB>/absolute/path/AD2.fastq
N1<TAB>/absolute/path/N1.fastq
```

Sample IDs must be unique and may contain letters, numbers, `.`, `_` and `-`.

### Groups and comparisons

Each sample must map to a biological group:

```groovy
sample_stage_mapping = [
  ["AD1", "AD"],
  ["AD2", "AD"],
  ["N1",  "Normal"]
]
```

Comparisons use `[group1, group2, comparison_name]`:

```groovy
comparisons = [
  ["AD", "Normal", "AD_vs_Normal"]
]
```

Comparison names must be unique and may contain letters, numbers, `.`, `_` and `-`. Both groups must occur in `sample_stage_mapping`.

## Provider-neutral LLM configuration

The optional LLM processes use the OpenAI Python SDK against an OpenAI-compatible Chat Completions API. Configure credentials with environment variables:

```bash
export LLM_BASE_URL='https://provider.example/v1'
export LLM_API_KEY='your-key'
export LLM_MODEL='provider-model-name'
```

- `LLM_BASE_URL` and `LLM_API_KEY` identify and authenticate the endpoint.
- `LLM_MODEL` is optional for backward compatibility, but should normally be set for a non-DeepSeek endpoint. OpenAI-compatible APIs require a model string and cannot reliably infer one from a multi-model URL.
- If `LLM_API_KEY` is absent, LLM steps are disabled automatically and all non-LLM analyses continue.
- Use `--llm_enabled false` to disable LLM calls explicitly.
- Enabling LLM without a key or a valid absolute HTTP(S) base URL fails during workflow validation.

The API key is deliberately not stored in `nextflow.config`, Nextflow parameters, task command lines, metadata or committed `.env` files. Do not pass real credentials on the command line.

Some output directories and JSON filenames retain historical `deepseek_*` names so existing downstream consumers do not break. These names do not restrict the configured provider; provider-neutral CLI entry points live in [`bin/llm/`](bin/llm/).

## Running the pipeline

Select the compatible R interpreter when more than one R installation is available:

```bash
export RSCRIPT_BIN='/usr/bin/Rscript'
"$RSCRIPT_BIN" --version
```

Run the full workflow:

```bash
nextflow run main.nf -resume
```

Run without any LLM requests:

```bash
nextflow run main.nf -resume --llm_enabled false
```

Useful structure-only checks:

```bash
NXF_OFFLINE=true nextflow lint main.nf nextflow.config modules/*.nf
NXF_OFFLINE=true nextflow run main.nf -preview --llm_enabled false
```

Nextflow caching is enabled. Keep `work/` and use `-resume` to reuse successfully completed tasks.

## Output layout

The default output root is `results/`:

```text
results/
├── InputAndPreprocess/
├── AlignmentAndTranscriptReconstruction/
├── PreliminarySingleCellAnalysis/
├── AlternativeSplicing/
└── IntegrativeAnalysis/
```

Each module publishes its reports, tables, plots and reusable objects below its corresponding directory. Large intermediate files remain in the Nextflow `work/` directory.

## R environment

The R scripts use packages including Seurat, tidyverse, data.table, harmony, rtracklayer, IsoformSwitchAnalyzeR, DESeq2, muscat, BSgenome.Hsapiens.UCSC.hg38, Biostrings, ggVennDiagram and UpSetR.

The current [`envs/sclong.yml`](envs/sclong.yml) mainly locks Python and command-line packages. For production deployment, create a dedicated R/Bioconductor environment matched to the server's R version and set `RSCRIPT_BIN` to that environment's `Rscript`. Avoid installing packages dynamically during a production run.

## Validation and compatibility

The optimized version has been checked with:

- Nextflow lint and preview
- Python 3.11 and the workflow's Python 3.7 environment
- Parsing and CLI startup for all eight R scripts
- Conda YAML parsing
- Offline regression cases for malformed FASTQ and metadata, NanoStats and SQANTI parsing, cell-annotation parsing, empty differential results, assembly error outputs and ORF ambiguity handling

The optimization did not add new biological analysis stages. To preserve the existing analysis scope, isoform-switch preparation still uses the first entry in `params.comparisons`, while dominant-isoform and integrative differential analyses iterate over all comparisons.

## Security

- Never commit API keys, tokens or private endpoint credentials.
- Keep local secret files outside the repository; `.env*` files are ignored by default.
- If a credential has ever been committed, deleting it in a later commit is not sufficient. Revoke and rotate the credential immediately, then clean the Git history if required.

## License

This project is distributed under the terms in [`LICENSE`](LICENSE).
