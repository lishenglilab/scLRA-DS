#!/usr/bin/env python3
"""
Assembly QC analysis using an OpenAI-compatible API
Analyzes transcriptome assembly metrics and generates evaluation report
"""

import argparse
import json
import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from llm.client import add_llm_arguments, client_from_args, response_text

def parse_assembly_files(isoquant_dir):
    """
    Parse assembly metrics from IsoQuant output files
    
    Parameters:
    isoquant_dir (str): Path to IsoQuant output directory
    
    Returns:
    dict: Parsed assembly metrics
    """
    metrics = {}
    
    # Define file paths
    gene_counts_file = os.path.join(isoquant_dir, "OUT", "OUT.discovered_gene_counts.tsv")
    transcript_counts_file = os.path.join(isoquant_dir, "OUT", "OUT.discovered_transcript_counts.tsv")
    transcript_models_file = os.path.join(isoquant_dir, "OUT", "OUT.transcript_models.gtf")
    novel_vs_known_file = os.path.join(isoquant_dir, "OUT", "OUT.novel_vs_known.SQANTI-like.tsv")
    
    print(f"Analyzing assembly from: {isoquant_dir}")
    
    # Check if necessary files exist
    required_files = [gene_counts_file, transcript_counts_file, transcript_models_file]
    for f in required_files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Required IsoQuant output not found: {f}")
    
    # 1. Load gene counts - FIXED: don't skip header with comment='#'
    try:
        gene_counts = pd.read_csv(gene_counts_file, sep="\t")
        
        # Check if column names were properly read
        print(f"Gene counts columns: {gene_counts.columns.tolist()}")
        print(f"Gene counts first few rows:")
        print(gene_counts.head())
        
        # If the first column is '#feature_id', rename it to 'feature_id'
        if '#feature_id' in gene_counts.columns:
            gene_counts = gene_counts.rename(columns={'#feature_id': 'feature_id'})
        
        # Exclude rows starting with __ and count non-zero
        valid_genes = gene_counts[~gene_counts['feature_id'].str.startswith('__', na=False)]
        non_zero_genes = valid_genes[valid_genes['count'] > 0]
        total_genes = len(non_zero_genes)
        metrics['total_genes'] = int(total_genes)
        
        # Count novel genes
        novel_genes = sum(non_zero_genes['feature_id'].str.startswith('novel_gene_', na=False))
        metrics['novel_genes'] = int(novel_genes)
        
        print(f"  Total genes: {total_genes}")
        print(f"  Novel genes: {novel_genes}")
    except Exception as e:
        print(f"ERROR loading gene counts: {e}")
        print(f"Gene counts file content (first 10 lines):")
        with open(gene_counts_file, 'r') as f:
            for i, line in enumerate(f):
                if i < 10:
                    print(f"  Line {i}: {line.strip()}")
                else:
                    break
        metrics['total_genes'] = 0
        metrics['novel_genes'] = 0
    
    # 2. Load transcript counts - FIXED: don't skip header with comment='#'
    try:
        transcript_counts = pd.read_csv(transcript_counts_file, sep="\t")
        
        # Check if column names were properly read
        print(f"Transcript counts columns: {transcript_counts.columns.tolist()}")
        print(f"Transcript counts first few rows:")
        print(transcript_counts.head())
        
        # If the first column is '#feature_id', rename it to 'feature_id'
        if '#feature_id' in transcript_counts.columns:
            transcript_counts = transcript_counts.rename(columns={'#feature_id': 'feature_id'})
        
        valid_transcripts = transcript_counts[~transcript_counts['feature_id'].str.startswith('__', na=False)]
        non_zero_transcripts = valid_transcripts[valid_transcripts['count'] > 0]
        total_transcripts = len(non_zero_transcripts)
        metrics['total_transcripts'] = int(total_transcripts)
        
        # Count novel transcripts
        novel_transcripts = sum(non_zero_transcripts['feature_id'].str.startswith('transcript', na=False))
        metrics['novel_transcripts'] = int(novel_transcripts)
        
        print(f"  Total transcripts: {total_transcripts}")
        print(f"  Novel transcripts: {novel_transcripts}")
    except Exception as e:
        print(f"ERROR loading transcript counts: {e}")
        print(f"Transcript counts file content (first 10 lines):")
        with open(transcript_counts_file, 'r') as f:
            for i, line in enumerate(f):
                if i < 10:
                    print(f"  Line {i}: {line.strip()}")
                else:
                    break
        metrics['total_transcripts'] = 0
        metrics['novel_transcripts'] = 0
    
    # 3. Calculate transcript lengths from GTF
    transcript_lengths = []
    try:
        with open(transcript_models_file, 'r') as f:
            current_transcript = None
            current_length = 0
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                feature_type = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                
                # Parse attributes
                attrs = {}
                for attr in parts[8].split(';'):
                    attr = attr.strip()
                    if ' ' in attr:
                        key, value = attr.split(' ', 1)
                        attrs[key] = value.strip('"')
                
                if feature_type == 'transcript':
                    if current_transcript and current_length > 0:
                        transcript_lengths.append(current_length)
                    current_transcript = attrs.get('transcript_id', '')
                    current_length = 0
                elif feature_type == 'exon':
                    if current_transcript:
                        current_length += (end - start + 1)
            
            # Add last transcript
            if current_transcript and current_length > 0:
                transcript_lengths.append(current_length)
        
        if transcript_lengths:
            metrics['transcript_length_min'] = int(min(transcript_lengths))
            metrics['transcript_length_max'] = int(max(transcript_lengths))
            metrics['transcript_length_median'] = int(np.median(transcript_lengths))
            metrics['transcript_length_mean'] = float(np.mean(transcript_lengths))
            print(f"  Transcript length range: {metrics['transcript_length_min']} - {metrics['transcript_length_max']} bp")
            print(f"  Median transcript length: {metrics['transcript_length_median']} bp")
            print(f"  Number of transcripts with length calculated: {len(transcript_lengths)}")
        else:
            print("  WARNING: No transcript lengths calculated")
            metrics['transcript_length_min'] = 0
            metrics['transcript_length_max'] = 0
            metrics['transcript_length_median'] = 0
            metrics['transcript_length_mean'] = 0.0
    except Exception as e:
        print(f"ERROR calculating transcript lengths: {e}")
        metrics['transcript_length_min'] = 0
        metrics['transcript_length_max'] = 0
        metrics['transcript_length_median'] = 0
        metrics['transcript_length_mean'] = 0.0
    
    # 4. Analyze novel transcript types if file exists
    novel_type_counts = {}
    if os.path.exists(novel_vs_known_file):
        try:
            with open(novel_vs_known_file, "r") as handle:
                first_field = handle.readline().split("\t", 1)[0].lstrip("#").strip()

            # IsoQuant 3.7 outputs can be either headered or headerless. Only
            # columns 1 and 6 are required for this metric.
            has_header = first_field == "isoform"
            novel_data = pd.read_csv(
                novel_vs_known_file,
                sep="\t",
                header=0 if has_header else None,
                usecols=[0, 5],
                names=None if has_header else ["isoform", "structural_category"],
            )
            if has_header:
                novel_data.columns = [str(col).lstrip("#") for col in novel_data.columns]

            novel_transcripts_data = novel_data[
                novel_data["isoform"].astype(str).str.startswith("transcript", na=False)
            ]
            novel_type_counts = (
                novel_transcripts_data["structural_category"].value_counts().to_dict()
            )
            print(f"  Novel transcript types: {len(novel_type_counts)} categories")
        except Exception as e:
            novel_type_counts = {"error": f"Failed to parse novel types: {str(e)}"}
            print(f"  ERROR parsing novel types: {e}")
    else:
        print(f"  Novel vs known file not found, skipping type analysis")
    
    metrics['novel_transcript_types'] = novel_type_counts
    
    # 5. Calculate additional metrics
    if metrics['total_genes'] > 0:
        metrics['transcripts_per_gene'] = float(
            metrics['total_transcripts'] / metrics['total_genes']
        )
    else:
        metrics['transcripts_per_gene'] = 0.0
    # Keep the historical JSON key for existing result consumers. Its value
    # has always represented transcripts per gene despite the old name.
    metrics['genes_per_transcript'] = metrics['transcripts_per_gene']
    
    # Calculate percentage of novel features
    if metrics['total_genes'] > 0:
        metrics['novel_genes_percent'] = float(metrics['novel_genes'] / metrics['total_genes'] * 100)
    else:
        metrics['novel_genes_percent'] = 0.0
    
    if metrics['total_transcripts'] > 0:
        metrics['novel_transcripts_percent'] = float(metrics['novel_transcripts'] / metrics['total_transcripts'] * 100)
    else:
        metrics['novel_transcripts_percent'] = 0.0
    
    # Debug: print summary
    print("\n=== SUMMARY ===")
    print(f"Total genes: {metrics['total_genes']}")
    print(f"Total transcripts: {metrics['total_transcripts']}")
    print(f"Novel genes: {metrics['novel_genes']} ({metrics['novel_genes_percent']:.1f}%)")
    print(f"Novel transcripts: {metrics['novel_transcripts']} ({metrics['novel_transcripts_percent']:.1f}%)")
    
    return metrics

def generate_assembly_report(metrics, client, species="Human", tissue="Skin", model="deepseek-reasoner"):
    """
    Generate an assembly quality report from the parsed metrics
    
    Parameters:
    metrics (dict): Dictionary containing assembly metrics
    client: Configured OpenAI-compatible client
    species (str): The species of the sample (default is "Human")
    tissue (str): The tissue type of the sample (default is "Skin")
    model (str): Model name exposed by the configured endpoint
    
    Returns:
    str: The content of the generated QC report
    """
    # Format novel transcript types for display
    novel_types_str = "\n"
    if metrics['novel_transcript_types']:
        for type_name, count in metrics['novel_transcript_types'].items():
            novel_types_str += f"  - {type_name}: {count}\n"
    else:
        novel_types_str = "  No novel transcript type information available\n"
    
    # Construct the system prompt
    system_prompt = f"""
You are a bioinformatics expert analyzing single-cell long-read RNA-seq transcriptome assembly results.

Please evaluate the following transcriptome assembly metrics from {species} {tissue} tissue and provide:
1. An overall assessment of assembly quality
2. Insights into transcript discovery and novelty
3. Recommendations for downstream analysis
4. Potential issues or areas for improvement

ASSEMBLY METRICS:
- Total genes detected: {metrics.get('total_genes', 'N/A')}
- Total transcripts detected: {metrics.get('total_transcripts', 'N/A')}
- Novel genes discovered: {metrics.get('novel_genes', 'N/A')} ({metrics.get('novel_genes_percent', 0):.1f}% of total)
- Novel transcripts discovered: {metrics.get('novel_transcripts', 'N/A')} ({metrics.get('novel_transcripts_percent', 0):.1f}% of total)
- Transcript length range: {metrics.get('transcript_length_min', 'N/A')} - {metrics.get('transcript_length_max', 'N/A')} bp
- Median transcript length: {metrics.get('transcript_length_median', 'N/A')} bp
- Mean transcript length: {metrics.get('transcript_length_mean', 'N/A'):.1f} bp
- Transcripts per gene ratio: {metrics.get('transcripts_per_gene', metrics.get('genes_per_transcript', 0)):.2f}

NOVEL TRANSCRIPT TYPES:
{novel_types_str}

Please consider:
1. Are the gene and transcript counts reasonable for single-cell long-read data?
2. How does the novel transcript/gene discovery rate compare to expectations for this technology?
3. Are transcript lengths within expected ranges for Nanopore sequencing?
4. What do the novel transcript type distributions suggest about assembly quality (e.g., completeness, fragmentation)?
5. Any red flags or exceptional observations that need attention?

Provide a comprehensive but concise report suitable for researchers, including specific recommendations for quality improvement if needed.
"""

    # Send the request to the configured OpenAI-compatible API.
    messages = [{"role": "system", "content": system_prompt}]
    
    completion = client.chat.completions.create(
        model=model,
        messages=messages
    )

    # Get the response content from the model
    response_content = response_text(completion)
    
    return response_content

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Generate an assembly QC report using an OpenAI-compatible API")
    
    # Arguments
    parser.add_argument('--isoquant_dir', required=True, help='Path to IsoQuant output directory')
    parser.add_argument('--output_text_file', required=True, help='Path to save the output QC report')
    parser.add_argument('--species', default="Human", help='The species of the sample (default: Human)')
    parser.add_argument('--tissue', default="Skin", help='The tissue type of the sample (default: Skin)')
    parser.add_argument('--save_metrics', action='store_true', help='Save metrics to JSON file')
    add_llm_arguments(parser)

    # Parse the arguments
    args = parser.parse_args()
    
    try:
        client = client_from_args(args)
    except ValueError as exc:
        parser.error(str(exc))
    
    # Parse the assembly files. This is an optional reporting step, so keep an
    # unexpected IsoQuant layout from aborting the biological workflow.
    try:
        metrics = parse_assembly_files(args.isoquant_dir)
    except Exception as exc:
        metrics = {"status": "error", "error": str(exc)}
        report_content = (
            "Failed to parse the IsoQuant outputs for the assembly QC report.\n"
            f"Error: {exc}\n"
        )
        if args.save_metrics:
            metrics_file = args.output_text_file.replace('.txt', '_metrics.json')
            with open(metrics_file, 'w') as f:
                json.dump(metrics, f, indent=2)
        with open(args.output_text_file, 'w') as f:
            f.write(report_content)
        print(report_content, file=sys.stderr)
        return
    
    # Save metrics if requested
    if args.save_metrics:
        metrics_file = args.output_text_file.replace('.txt', '_metrics.json')
        with open(metrics_file, 'w') as f:
            json.dump(metrics, f, indent=2)
        print(f"Metrics saved to: {metrics_file}")
    
    # Keep this optional reporting step from aborting the biological workflow
    # when the configured provider is temporarily unavailable.
    try:
        report_content = generate_assembly_report(
            metrics,
            client=client,
            species=args.species,
            tissue=args.tissue,
            model=args.model,
        )
    except Exception as exc:
        report_content = (
            "Failed to generate the assembly QC report via the configured LLM API.\n"
            f"Error: {exc}\n"
        )
        print(report_content, file=sys.stderr)
    
    # Save the result to the output text file
    with open(args.output_text_file, 'w') as f:
        f.write(report_content)
    
    print(f"Assembly quality report saved to: {args.output_text_file}")

if __name__ == '__main__':
    main()
