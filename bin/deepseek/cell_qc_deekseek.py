#!/usr/bin/env python3
"""Single-cell QC reporting through an OpenAI-compatible LLM endpoint."""

import pandas as pd
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from llm.client import add_llm_arguments, client_from_args, response_text

def generate_qc_report(args):
    """
    Generate a quality control report for single-cell RNA-seq data.
    
    Parameters:
    args (Namespace): The arguments passed from the command line.
    
    Returns:
    str: The content of the generated QC report.
    """
    client = client_from_args(args)

    # Read the CSV/TSV file containing the QC data
    if args.csv_file_path.endswith('.tsv'):
        df = pd.read_csv(args.csv_file_path, sep='\t')
    else:
        df = pd.read_csv(args.csv_file_path)
    if df.empty:
        raise ValueError("QC metrics input does not contain any data rows")
    
    # Debug: print the columns to verify
    print(f"DataFrame columns: {df.columns.tolist()}")
    print(f"DataFrame shape: {df.shape}")
    print(df.head())
    
    # Map the columns from R script to the expected metric names
    # R script column names:
    # total_gene_umi_mean, total_gene_umi_min, total_gene_umi_max, total_gene_umi_std,
    # genes_detected_mean, genes_detected_min, genes_detected_max, genes_detected_std,
    # total_tr_umi_mean, total_tr_umi_min, total_tr_umi_max, total_tr_umi_std,
    # transcripts_detected_mean, transcripts_detected_min, transcripts_detected_max, transcripts_detected_sd,
    # mt_rna_mean, mt_rna_min, mt_rna_max, mt_rna_std,
    # gene_transcript_corr, doublet_rate
    
    # Fixed metrics_columns mapping
    metrics_columns = {
        "Total Gene UMI Counts": {
            'mean': 'total_gene_umi_mean',
            'min': 'total_gene_umi_min', 
            'max': 'total_gene_umi_max',
            'std': 'total_gene_umi_std'
        },
        "Genes Detected per Cell": {
            'mean': 'genes_detected_mean',
            'min': 'genes_detected_min',
            'max': 'genes_detected_max',
            'std': 'genes_detected_std'
        },
        "Total Transcript UMI Counts": {
            'mean': 'total_tr_umi_mean',
            'min': 'total_tr_umi_min',
            'max': 'total_tr_umi_max',
            'std': 'total_tr_umi_std'
        },
        "Transcripts Detected per Cell": {
            'mean': 'transcripts_detected_mean',
            'min': 'transcripts_detected_min',
            'max': 'transcripts_detected_max',
            'std': 'transcripts_detected_sd'
        },
        "Mitochondrial RNA Percentage": {
            'mean': 'mt_rna_mean',
            'min': 'mt_rna_min',
            'max': 'mt_rna_max',
            'std': 'mt_rna_std'
        },
        "Gene-Transcript Expression Correlation": {
            'value': 'gene_transcript_corr'
        },
        "Doublet Rate": {
            'value': 'doublet_rate'
        }
    }
    
    # Check if all required columns exist
    missing_columns = []
    for metric, columns in metrics_columns.items():
        for key, col_name in columns.items():
            if col_name not in df.columns:
                missing_columns.append(col_name)
    
    if missing_columns:
        print(f"WARNING: Missing columns: {missing_columns}")
        print("Available columns:", df.columns.tolist())
    
    # Build the system prompt for the API call
    system_prompt = f"""You are a bioinformatics expert analyzing single-cell RNA-seq data quality.

Please evaluate the following single-cell RNA-seq quality control metrics from {args.species} {args.tissue} tissue and provide:
1. An overall assessment of data quality
2. Analysis of each metric compared to expected values for long-read single-cell RNA-seq
3. Recommendations for data filtering or processing if needed
4. Suggestions for downstream analysis based on the data quality

SINGLE-CELL QC METRICS:

"""
    
    # Add each metric to the prompt
    for metric_name, columns in metrics_columns.items():
        if 'value' in columns:
            # Single value metric
            col_name = columns['value']
            if col_name in df.columns:
                value = df[col_name].iloc[0]
                if pd.notna(value):
                    system_prompt += f"{metric_name}: {value:.3f}\n"
                else:
                    system_prompt += f"{metric_name}: Not available\n"
            else:
                system_prompt += f"{metric_name}: Not available\n"
        else:
            # Multi-value metric (mean, min, max, std)
            values = []
            for stat, col_name in columns.items():
                if col_name in df.columns:
                    value = df[col_name].iloc[0]
                    if pd.notna(value):
                        if stat in ['mean', 'std']:
                            values.append(f"{stat}: {value:.2f}")
                        else:
                            values.append(f"{stat}: {value:.0f}")
            if values:
                system_prompt += f"{metric_name}: ({', '.join(values)})\n"
            else:
                system_prompt += f"{metric_name}: Not available\n"
    
    system_prompt += """

Please consider the following aspects in your evaluation:
1. **UMI Counts**: Are the total UMI counts per cell reasonable for long-read single-cell data?
2. **Genes/Transcripts Detected**: Is the number of detected genes/transcripts per cell adequate?
3. **Mitochondrial RNA**: Is the mitochondrial RNA percentage within acceptable limits (typically < 20%)?
4. **Doublet Rate**: Is the doublet rate reasonable for single-cell experiments?
5. **Gene-Transcript Correlation**: Is there good correlation between gene and transcript expression?
6. **Data Consistency**: Are all metrics consistent with high-quality single-cell data?

Provide a comprehensive but concise report suitable for researchers, including specific recommendations if data quality issues are detected.
"""

    print("Generated prompt preview (first 1000 chars):")
    print(system_prompt[:1000] + "...")

    # Create the messages to send to the API
    messages = [{"role": "system", "content": system_prompt}]
    
    # Call the configured OpenAI-compatible API to generate the report.
    try:
        completion = client.chat.completions.create(
            model=args.model,
            messages=messages
        )

        # Get the response content from the model
        response_content = response_text(completion)
        
        # Save the result to the output text file
        with open(args.output_text_file, 'w') as f:
            f.write(response_content)
        
        print(f"Single-cell quality control report saved to {args.output_text_file}")
        
        return response_content
    
    except Exception as e:
        print(f"Error calling LLM API: {e}")
        # Create a basic error report
        error_report = f"Failed to generate QC report via the configured LLM API.\nError: {str(e)}\n\nPlease check the endpoint, API key, and network connection."
        with open(args.output_text_file, 'w') as f:
            f.write(error_report)
        return error_report


def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Generate a QC report for single-cell RNA-seq data.")
    
    # Arguments
    parser.add_argument('--csv_file_path', required=True, help='Path to the CSV/TSV file containing QC data')
    parser.add_argument('--output_text_file', required=True, help='Path to save the output QC report')
    parser.add_argument('--species', default="Human", help='The species of the sample (default: Human)')
    parser.add_argument('--tissue', default="Skin", help='The tissue type of the sample (default: Skin)')
    add_llm_arguments(parser)
    
    # Parse the arguments
    args = parser.parse_args()

    # Check if input file exists
    import os
    if not os.path.exists(args.csv_file_path):
        print(f"ERROR: Input file not found: {args.csv_file_path}")
        sys.exit(1)
    
    # Generate the QC report
    try:
        generate_qc_report(args)
    except ValueError as exc:
        parser.error(str(exc))


if __name__ == '__main__':
    main()
