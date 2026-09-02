#!/usr/bin/env python3
"""LLM interpretation of genes shared by four differential analyses."""

import pandas as pd
import argparse
import sys
import os
import json
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from llm.client import (
    add_llm_arguments,
    client_from_args,
    response_text,
    response_total_tokens,
)

def parse_differential_results(result_dir, comparison_name):
    """
    Parse differential analysis results from the four TXT files
    
    Parameters:
    result_dir (str): Path to differential analysis results directory
    comparison_name (str): Comparison name (e.g., AD_vs_Normal)
    
    Returns:
    dict: Combined gene information from all four analyses
    """
    try:
        # Define file paths
        files = {
            'DDTU': os.path.join(result_dir, f"Differential_Dominant_Transcript_Usage_{comparison_name}.txt"),
            'DTE': os.path.join(result_dir, f"Differential_Transcript_Expression_{comparison_name}.txt"),
            'DGE': os.path.join(result_dir, f"Differential_Gene_Expression_{comparison_name}.txt"),
            'DGU': os.path.join(result_dir, f"Differential_Gene_Usage_{comparison_name}.txt")
        }
        
        print(f"Looking for files in: {result_dir}")
        for analysis, filepath in files.items():
            print(f"  {analysis}: {filepath}")
            print(f"    Exists: {os.path.exists(filepath)}")
        
        # Check if files exist
        missing_files = [k for k, v in files.items() if not os.path.exists(v)]
        if missing_files:
            print(f"ERROR: Missing files: {missing_files}")
            return None
        
        # Read each file
        results = {}
        for analysis, filepath in files.items():
            try:
                print(f"Reading {analysis} file...")
                df = pd.read_csv(filepath, sep='\t')
                print(f"  Success! Shape: {df.shape}, Columns: {list(df.columns)}")
                
                # Standardize column names for easier processing
                if 'gene_name' not in df.columns and 'gene' in df.columns:
                    df['gene_name'] = df['gene']
                
                results[analysis] = df
            except Exception as e:
                print(f"ERROR reading {analysis} file {filepath}: {str(e)}")
                return None
        
        return results
        
    except Exception as e:
        print(f"ERROR parsing differential results: {str(e)}")
        return None

def find_intersection_genes(differential_results):
    """
    Find genes that appear in all four analyses
    
    Parameters:
    differential_results (dict): Dictionary with results from four analyses
    
    Returns:
    list: Genes in the intersection of all four analyses
    """
    print("\nFinding genes in intersection of all four analyses...")
    
    # Extract significant genes from each analysis
    gene_sets = {}
    
    # DDTU
    ddtu_df = differential_results.get('DDTU')
    if ddtu_df is not None:
        print(f"DDTU columns: {list(ddtu_df.columns)}")
        # Look for type column or change column
        type_col = 'type' if 'type' in ddtu_df.columns else ('change' if 'change' in ddtu_df.columns else None)
        gene_col = 'gene_name' if 'gene_name' in ddtu_df.columns else 'gene'
        
        if type_col and gene_col in ddtu_df.columns:
            ddtu_sig = ddtu_df[ddtu_df[type_col] != 'notSig'][gene_col].dropna().unique()
            gene_sets['DDTU'] = set(ddtu_sig)
            print(f"  DDTU significant genes: {len(gene_sets['DDTU'])}")
    
    # DTE
    dte_df = differential_results.get('DTE')
    if dte_df is not None:
        print(f"DTE columns: {list(dte_df.columns)}")
        type_col = 'change' if 'change' in dte_df.columns else 'type'
        gene_col = 'gene_name' if 'gene_name' in dte_df.columns else 'gene'
        
        if type_col in dte_df.columns and gene_col in dte_df.columns:
            dte_sig = dte_df[dte_df[type_col] != 'Stable'][gene_col].dropna().unique()
            gene_sets['DTE'] = set(dte_sig)
            print(f"  DTE significant genes: {len(gene_sets['DTE'])}")
    
    # DGE
    dge_df = differential_results.get('DGE')
    if dge_df is not None:
        print(f"DGE columns: {list(dge_df.columns)}")
        type_col = 'change' if 'change' in dge_df.columns else 'type'
        gene_col = 'gene' if 'gene' in dge_df.columns else 'gene_name'
        
        if type_col in dge_df.columns and gene_col in dge_df.columns:
            dge_sig = dge_df[dge_df[type_col] != 'Stable'][gene_col].dropna().unique()
            gene_sets['DGE'] = set(dge_sig)
            print(f"  DGE significant genes: {len(gene_sets['DGE'])}")
    
    # DGU
    dgu_df = differential_results.get('DGU')
    if dgu_df is not None:
        print(f"DGU columns: {list(dgu_df.columns)}")
        type_col = 'type' if 'type' in dgu_df.columns else 'change'
        gene_col = 'gene_name' if 'gene_name' in dgu_df.columns else 'gene'
        
        if type_col in dgu_df.columns and gene_col in dgu_df.columns:
            dgu_sig = dgu_df[dgu_df[type_col] != 'notSig'][gene_col].dropna().unique()
            gene_sets['DGU'] = set(dgu_sig)
            print(f"  DGU significant genes: {len(gene_sets['DGU'])}")
    
    # Find intersection
    if len(gene_sets) == 4:
        intersection_genes = set.intersection(*gene_sets.values())
        print(f"\nGenes in intersection of all 4 analyses: {len(intersection_genes)}")
        return sorted(intersection_genes), gene_sets
    else:
        print(f"ERROR: Not all gene sets available. Found {len(gene_sets)} sets.")
        print(f"Available sets: {list(gene_sets.keys())}")
        return [], gene_sets

def get_gene_details(gene_name, differential_results):
    """
    Get detailed information about a gene from all four analyses
    
    Parameters:
    gene_name (str): Gene symbol
    differential_results (dict): Results from four analyses
    
    Returns:
    dict: Gene details from each analysis
    """
    details = {
        'gene_name': gene_name,
        'DDTU': None,
        'DTE': None,
        'DGE': None,
        'DGU': None
    }
    
    # DDTU details
    ddtu_df = differential_results.get('DDTU')
    if ddtu_df is not None:
        gene_col = 'gene_name' if 'gene_name' in ddtu_df.columns else 'gene'
        if gene_col in ddtu_df.columns:
            gene_data = ddtu_df[ddtu_df[gene_col] == gene_name]
            if not gene_data.empty:
                type_col = 'type' if 'type' in gene_data.columns else 'change'
                significant = gene_data[gene_data[type_col] != 'notSig']
                row = significant.iloc[0] if not significant.empty else gene_data.iloc[0]
                details['DDTU'] = {
                    'type': row[type_col],
                    'transcript_id': row['transcript_id'] if 'transcript_id' in gene_data.columns else 'Unknown'
                }
    
    # DTE details
    dte_df = differential_results.get('DTE')
    if dte_df is not None:
        gene_col = 'gene_name' if 'gene_name' in dte_df.columns else 'gene'
        if gene_col in dte_df.columns:
            gene_data = dte_df[dte_df[gene_col] == gene_name]
            if not gene_data.empty:
                type_col = 'change' if 'change' in gene_data.columns else 'type'
                significant = gene_data[gene_data[type_col] != 'Stable']
                row = significant.iloc[0] if not significant.empty else gene_data.iloc[0]
                details['DTE'] = {
                    'change': row[type_col],
                    'avg_log2FC': float(row['avg_log2FC']) if 'avg_log2FC' in gene_data.columns else 0,
                    'p_val_adj': float(row['p_val_adj']) if 'p_val_adj' in gene_data.columns else 1
                }
    
    # DGE details
    dge_df = differential_results.get('DGE')
    if dge_df is not None:
        gene_col = 'gene' if 'gene' in dge_df.columns else 'gene_name'
        if gene_col in dge_df.columns:
            gene_data = dge_df[dge_df[gene_col] == gene_name]
            if not gene_data.empty:
                type_col = 'change' if 'change' in gene_data.columns else 'type'
                significant = gene_data[gene_data[type_col] != 'Stable']
                row = significant.iloc[0] if not significant.empty else gene_data.iloc[0]
                details['DGE'] = {
                    'change': row[type_col],
                    'avg_log2FC': float(row['avg_log2FC']) if 'avg_log2FC' in gene_data.columns else 0,
                    'p_val_adj': float(row['p_val_adj']) if 'p_val_adj' in gene_data.columns else 1
                }
    
    # DGU details
    dgu_df = differential_results.get('DGU')
    if dgu_df is not None:
        gene_col = 'gene_name' if 'gene_name' in dgu_df.columns else 'gene'
        if gene_col in dgu_df.columns:
            gene_data = dgu_df[dgu_df[gene_col] == gene_name]
            if not gene_data.empty:
                type_col = 'type' if 'type' in gene_data.columns else 'change'
                significant = gene_data[gene_data[type_col] != 'notSig']
                row = significant.iloc[0] if not significant.empty else gene_data.iloc[0]
                details['DGU'] = {
                    'type': row[type_col],
                    'lgFC': float(row['lgFC']) if 'lgFC' in gene_data.columns else 0,
                    'ident1_pct': float(row['ident1_pct']) if 'ident1_pct' in gene_data.columns else 0,
                    'ident2_pct': float(row['ident2_pct']) if 'ident2_pct' in gene_data.columns else 0
                }
    
    return details

def generate_prompt(intersection_genes, gene_details, study_description, comparison_name):
    """
    Generate a prompt for the configured LLM API
    
    Parameters:
    intersection_genes (list): List of genes in the intersection of all four analyses
    gene_details (dict): Detailed information for each gene
    study_description (str): Study context
    comparison_name (str): Comparison name
    
    Returns:
    str: Formatted prompt
    """
    # Extract comparison groups
    if '_vs_' in comparison_name:
        group1, group2 = comparison_name.split('_vs_', 1)
    else:
        group1, group2 = "Condition1", "Condition2"
    
    # Format gene information
    genes_info = []
    for i, gene in enumerate(intersection_genes[:15], 1):  # Limit to first 15 genes for prompt length
        details = gene_details.get(gene, {})
        
        gene_info = f"{i}. {gene}\n"
        
        if details.get('DDTU'):
            ddtu = details['DDTU']
            gene_info += f"   • DDTU (Transcript Usage): {ddtu.get('type', 'Unknown')} change"
            if 'transcript_id' in ddtu and ddtu['transcript_id'] != 'Unknown':
                gene_info += f" (Transcript: {ddtu['transcript_id']})"
            gene_info += "\n"
        
        if details.get('DTE'):
            dte = details['DTE']
            change = dte.get('change', 'Unknown')
            log2fc = dte.get('avg_log2FC', 0)
            direction = "up" if change == 'Up' else "down" if change == 'Down' else ""
            gene_info += f"   • DTE (Transcript Expression): {direction}regulated (log2FC: {log2fc:.2f})\n"
        
        if details.get('DGE'):
            dge = details['DGE']
            change = dge.get('change', 'Unknown')
            log2fc = dge.get('avg_log2FC', 0)
            direction = "up" if change == 'Up' else "down" if change == 'Down' else ""
            gene_info += f"   • DGE (Gene Expression): {direction}regulated (log2FC: {log2fc:.2f})\n"
        
        if details.get('DGU'):
            dgu = details['DGU']
            gtype = dgu.get('type', 'Unknown')
            lgfc = dgu.get('lgFC', 0)
            gene_info += f"   • DGU (Gene Usage): {gtype}-biased (log2FC: {lgfc:.2f})\n"
        
        genes_info.append(gene_info)
    
    prompt = f"""You are a senior bioinformatics expert and molecular biologist. Please provide a comprehensive biological interpretation of genes that are significantly altered across all four analysis types in our single-cell transcriptomics study.

STUDY CONTEXT:
{study_description}

ANALYSIS DETAILS:
We are comparing {group1} vs {group2} ({comparison_name}).
The analysis includes four complementary approaches:
1. DDTU (Differential Dominant Transcript Usage): Changes in the main transcript isoform used
2. DTE (Differential Transcript Expression): Changes in individual transcript expression levels
3. DGE (Differential Gene Expression): Changes in overall gene expression
4. DGU (Differential Gene Usage): Changes in the proportion of cells expressing each gene

KEY FINDING:
We have identified {len(intersection_genes)} genes that are significantly altered in ALL FOUR analysis types. These genes represent the core molecular signature of the biological changes between {group1} and {group2}.

GENE-SPECIFIC INFORMATION (showing {min(15, len(intersection_genes))} of {len(intersection_genes)} genes):
{chr(10).join(genes_info)}

YOUR TASK:
Please provide a comprehensive biological interpretation covering:

1. OVERALL SIGNIFICANCE:
   - Why is it particularly meaningful that these genes appear in all four analysis types?
   - What does this multi-level convergence suggest about the robustness and biological importance of these changes?

2. BIOLOGICAL PATHWAYS AND FUNCTIONS:
   - Based on the specific genes identified, what biological pathways or cellular processes are most affected?
   - How do these pathways relate to the disease/condition being studied?

3. REGULATORY MECHANISMS:
   - What do the convergent changes at transcript usage, transcript expression, gene expression, and gene usage levels suggest about underlying regulatory mechanisms?
   - Are these changes likely driven by transcriptional regulation, post-transcriptional regulation, or both?

4. DISEASE/PHENOTYPE RELEVANCE:
   - How might these multi-level changes contribute to the phenotypic differences between {group1} and {group2}?
   - Do any of these genes have known relevance to the specific disease/condition?

5. THERAPEUTIC IMPLICATIONS:
   - Based on these findings, what might be promising therapeutic targets or intervention strategies?
   - Which of these genes or pathways might be most amenable to therapeutic modulation?

Please provide a detailed, well-organized response suitable for both biologists and bioinformaticians. Focus on biological insight rather than statistical description.
"""
    return prompt

def call_llm_api(prompt, client, model="deepseek-chat"):
    """
    Call an OpenAI-compatible API for biological interpretation
    
    Parameters:
    prompt (str): Prompt for the model
    client: Configured OpenAI-compatible client
    model (str): Model to use
    
    Returns:
    dict: API response and status
    """
    try:
        response = client.chat.completions.create(
            model=model,
            messages=[
                {
                    "role": "system", 
                    "content": "You are a senior bioinformatics expert and molecular biologist specializing in transcriptomics, gene regulation, and disease mechanisms."
                },
                {
                    "role": "user", 
                    "content": prompt
                }
            ],
            temperature=0.7,
            max_tokens=3000
        )
        
        interpretation = response_text(response)
        
        return {
            'status': 'success',
            'interpretation': interpretation,
            'model': model,
            'tokens_used': response_total_tokens(response)
        }
        
    except Exception as e:
        return {
            'status': 'error',
            'message': f'API call failed: {str(e)}'
        }

def main():
    parser = argparse.ArgumentParser(
        description='Generate biological interpretation for genes shared by four analyses using an OpenAI-compatible API'
    )
    
    # Required arguments
    parser.add_argument('--result_dir', required=True, 
                       help='Path to differential analysis results directory')
    parser.add_argument('--comparison_name', required=True,
                       help='Comparison name (e.g., AD_vs_Normal)')
    parser.add_argument('--output_report', required=True,
                       help='Path to save the biological interpretation report')
    parser.add_argument('--output_api', required=True,
                       help='Path to save the API response details')
    
    # Optional arguments
    parser.add_argument('--study_description', default='',
                       help='Study context description')
    add_llm_arguments(parser, default_model='deepseek-chat')
    
    args = parser.parse_args()
    
    try:
        client = client_from_args(args)
    except ValueError as exc:
        parser.error(str(exc))
    
    # Set default study description if not provided
    if not args.study_description:
        args.study_description = """This is a single-cell full-length transcriptome sequencing analysis comparing AD (atopic dermatitis) and Normal (normal skin cells) conditions. 
The study aims to identify differentially expressed transcripts and alternative splicing events in AD, 
exploring disease-associated molecular changes at multiple levels."""
    
    print(f"Analyzing genes in intersection of all four differential analyses for comparison: {args.comparison_name}")
    print(f"Results directory: {args.result_dir}")
    
    # Parse differential results
    print("\nParsing differential analysis results...")
    differential_results = parse_differential_results(args.result_dir, args.comparison_name)
    
    if differential_results is None:
        print("ERROR: Could not parse differential results")
        sys.exit(1)
    
    # Find genes in intersection of all four analyses
    intersection_genes, gene_sets = find_intersection_genes(differential_results)
    
    # Check if we found any genes
    if not intersection_genes:
        print("\nWARNING: No genes found in intersection of all four analyses.")
        print("This suggests that while there are significant changes at individual analysis levels,")
        print("no single gene shows convergent changes across all analysis types.")
        
        with open(args.output_report, 'w') as f:
            f.write("No genes were found in the intersection of all four analysis types.\n\n")
            f.write("Individual analysis statistics:\n")
            for analysis, gene_set in gene_sets.items():
                f.write(f"  {analysis}: {len(gene_set)} significant genes\n")
            f.write("\nInterpretation: The absence of genes in the intersection suggests that the molecular changes\n")
            f.write("between conditions are distributed across different regulatory layers rather than concentrated\n")
            f.write("in a few genes. This could indicate complex, multi-layered regulatory mechanisms at play.\n")
        
        with open(args.output_api, 'w') as f:
            f.write(json.dumps({
                'status': 'no_genes',
                'message': 'No genes found in intersection of all four analyses',
                'individual_counts': {k: len(v) for k, v in gene_sets.items()},
                'timestamp': datetime.now().isoformat()
            }, indent=2))
        summary_file = args.output_report.replace('.txt', '_summary.txt')
        with open(summary_file, 'w') as f:
            f.write("No genes were shared by all four differential analyses.\n")
            for analysis, gene_set in sorted(gene_sets.items()):
                f.write(f"{analysis}: {len(gene_set)} significant genes\n")
        sys.exit(0)
    
    print(f"\nFound {len(intersection_genes)} genes in intersection of all four analyses:")
    for i, gene in enumerate(sorted(intersection_genes), 1):
        print(f"  {i}. {gene}")
    
    # Get gene details
    print("\nExtracting gene details...")
    gene_details = {}
    for gene in intersection_genes:
        gene_details[gene] = get_gene_details(gene, differential_results)
    
    # Generate prompt
    print("Generating prompt for LLM API...")
    prompt = generate_prompt(intersection_genes, gene_details, args.study_description, args.comparison_name)
    
    # Save prompt for debugging
    prompt_file = args.output_report.replace('.txt', '_prompt.txt')
    with open(prompt_file, 'w') as f:
        f.write(prompt)
    print(f"Prompt saved to: {prompt_file}")
    
    # Call the configured LLM API.
    print("\nCalling LLM API...")
    api_result = call_llm_api(prompt, client, args.model)
    
    if api_result['status'] == 'success':
        interpretation = api_result['interpretation']
        
        # Save interpretation report
        with open(args.output_report, 'w') as f:
            f.write(interpretation)
        
        # Save API response details
        api_details = {
            'timestamp': datetime.now().isoformat(),
            'status': 'success',
            'comparison': args.comparison_name,
            'intersection_genes_count': len(intersection_genes),
            'intersection_genes': sorted(intersection_genes),
            'individual_counts': {k: len(v) for k, v in gene_sets.items()},
            'model': api_result['model'],
            'tokens_used': api_result['tokens_used'],
            'prompt_length': len(prompt),
            'response_length': len(interpretation)
        }
        
        with open(args.output_api, 'w') as f:
            json.dump(api_details, f, indent=2)
        
        print(f"\nBiological interpretation generated for {len(intersection_genes)} genes in intersection")
        print(f"Report saved to: {args.output_report}")
        
    else:
        print(f"\nAPI call failed: {api_result['message']}")
        
        with open(args.output_report, 'w') as f:
            f.write(f"API call failed: {api_result['message']}\n")
            f.write("\nNo biological interpretation was generated due to API error.\n")
        
        with open(args.output_api, 'w') as f:
            f.write(json.dumps({
                'status': 'api_error',
                'message': api_result['message'],
                'timestamp': datetime.now().isoformat()
            }, indent=2))
    
    # Create gene summary file
    summary_file = args.output_report.replace('.txt', '_summary.txt')
    
    summary = f"""Differential Analysis Interpretation Summary
===============================================
Comparison: {args.comparison_name}
Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

Individual Analysis Statistics:
{'-' * 50}"""
    
    for analysis, gene_set in gene_sets.items():
        summary += f"\n{analysis}: {len(gene_set)} significant genes"
    
    summary += f"\n\nGenes in Intersection of All Four Analyses:"
    summary += f"\n{'-' * 50}"
    summary += f"\nTotal genes in intersection: {len(intersection_genes)}\n"
    
    for i, gene in enumerate(sorted(intersection_genes), 1):
        details = gene_details.get(gene, {})
        
        summary += f"\n{i}. {gene}\n"
        
        if details.get('DDTU'):
            ddtu_type = details['DDTU'].get('type', 'Unknown')
            summary += f"   • DDTU: {ddtu_type}\n"
        
        if details.get('DTE'):
            dte_change = details['DTE'].get('change', 'Unknown')
            dte_fc = details['DTE'].get('avg_log2FC', 0)
            summary += f"   • DTE: {dte_change} (log2FC: {dte_fc:.2f})\n"
        
        if details.get('DGE'):
            dge_change = details['DGE'].get('change', 'Unknown')
            dge_fc = details['DGE'].get('avg_log2FC', 0)
            summary += f"   • DGE: {dge_change} (log2FC: {dge_fc:.2f})\n"
        
        if details.get('DGU'):
            dgu_type = details['DGU'].get('type', 'Unknown')
            dgu_fc = details['DGU'].get('lgFC', 0)
            summary += f"   • DGU: {dgu_type} (log2FC: {dgu_fc:.2f})\n"
    
    summary += f"\n{'=' * 50}\n"
    summary += "Note: These genes show significant changes in all four analysis types:\n"
    summary += "1. Differential Dominant Transcript Usage (DDTU)\n"
    summary += "2. Differential Transcript Expression (DTE)\n"
    summary += "3. Differential Gene Expression (DGE)\n"
    summary += "4. Differential Gene Usage (DGU)\n"
    summary += "\nThis convergence suggests they are core regulators of the observed biological changes."
    
    with open(summary_file, 'w') as f:
        f.write(summary)
    
    print(f"Gene summary saved to: {summary_file}")
    print("\nDifferential analysis interpretation completed successfully!")

if __name__ == '__main__':
    main()
