import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from llm.client import add_llm_arguments, client_from_args, response_text


NUMBER_PATTERN = r"[-+]?(?:\d[\d,]*)(?:\.\d+)?"


def parse_number(value):
    """Parse NanoPlot numbers that may contain commas and decimals."""
    return float(value.replace(",", ""))

def parse_nanostats_file(file_path):
    """
    Parse the NanoStats file and extract key sequencing metrics.
    
    Parameters:
    file_path (str): Path to the NanoStats text file.
    
    Returns:
    dict: A dictionary containing the parsed sequencing metrics.
    """
    metrics = {}

    with open(file_path, 'r') as file:
        content = file.read()

    # Regular expressions to match relevant metrics
    def safe_search(pattern, content, default=None):
        match = re.search(pattern, content)
        if match:
            return parse_number(match.group(1))
        return default

    # Extract the data using safe_search
    metrics['mean_read_length'] = safe_search(rf"Mean read length:\s*({NUMBER_PATTERN})", content)
    metrics['mean_read_quality'] = safe_search(rf"Mean read quality:\s*({NUMBER_PATTERN})", content)
    metrics['median_read_length'] = safe_search(rf"Median read length:\s*({NUMBER_PATTERN})", content)
    metrics['median_read_quality'] = safe_search(rf"Median read quality:\s*({NUMBER_PATTERN})", content)
    metrics['number_of_reads'] = int(safe_search(rf"Number of reads:\s*({NUMBER_PATTERN})", content, default=0))
    metrics['read_length_n50'] = safe_search(rf"Read length N50:\s*({NUMBER_PATTERN})", content)
    metrics['stdev_read_length'] = safe_search(rf"STDEV read length:\s*({NUMBER_PATTERN})", content)
    metrics['total_bases'] = safe_search(rf"Total bases:\s*({NUMBER_PATTERN})", content)

    # Parsing quality cutoff information
    quality_cutoffs = re.findall(
        rf">Q(\d+):\s*({NUMBER_PATTERN})\s*\(({NUMBER_PATTERN})%\)\s*({NUMBER_PATTERN})Mb",
        content,
    )
    metrics['quality_cutoffs'] = {}
    for cutoff in quality_cutoffs:
        quality_score, count, percentage, mb = cutoff
        metrics['quality_cutoffs'][f'Q{quality_score}'] = {
            'count': int(parse_number(count)),
            'percentage': parse_number(percentage),
            'megabases': parse_number(mb)
        }

    # Parsing top reads information (mean basecall quality)
    top_reads = re.findall(
        rf"(\d+):\s*({NUMBER_PATTERN})\s*\(({NUMBER_PATTERN})\)",
        content.split("Top 5 longest reads", 1)[0],
    )
    metrics['top_reads_quality'] = [
        (int(top[0]), parse_number(top[1]), int(parse_number(top[2])))
        for top in top_reads
    ]

    return metrics


def generate_qc_report(metrics, client, species="Human", tissue="Skin", model="deepseek-reasoner", output_text_file='qc_report.txt'):
    """
    Generate a quality control report from the parsed NanoStats data.
    
    Parameters:
    metrics (dict): Dictionary containing parsed sequencing metrics.
    client: Configured OpenAI-compatible client.
    species (str): The species of the sample (default is "Human").
    tissue (str): The tissue type of the sample (default is "Skin").
    model (str): Model name exposed by the configured endpoint.
    output_text_file (str): Path to save the generated QC report.
    
    Returns:
    str: The content of the generated QC report.
    """
    # Construct the system prompt based on the parsed data and provided species/tissue
    system_prompt = (
        f"Generate a quality control report for the following single-cell Nanopore long-read RNA-seq data from {species} {tissue} tissue. "
        "The report should evaluate each of the following metrics and provide a conclusion about the data quality, "
        "with suggestions for improving the quality if necessary:\n\n"
        f"Mean Read Length: {metrics['mean_read_length']} bp\n"
        f"Mean Read Quality: {metrics['mean_read_quality']}\n"
        f"Median Read Length: {metrics['median_read_length']} bp\n"
        f"Median Read Quality: {metrics['median_read_quality']}\n"
        f"Number of Reads: {metrics['number_of_reads']}\n"
        f"Read Length N50: {metrics['read_length_n50']} bp\n"
        f"STDEV Read Length: {metrics['stdev_read_length']} bp\n"
        f"Total Bases: {metrics['total_bases']} bases\n\n"
        "Quality Cutoffs:\n"
    )

    for cutoff, data in metrics['quality_cutoffs'].items():
        system_prompt += f"> {cutoff}: {data['count']} reads ({data['percentage']}%) {data['megabases']}Mb\n"

    system_prompt += "\nTop 5 Highest Mean Basecall Quality Scores and Their Read Lengths:\n"
    for i, (read_num, quality, length) in enumerate(metrics['top_reads_quality']):
        system_prompt += f"{i + 1}: {read_num} (Quality: {quality}, Length: {length} bp)\n"

    # Send the request to the configured OpenAI-compatible API.
    messages = [{"role": "system", "content": system_prompt}]
    
    completion = client.chat.completions.create(
        model=model,
        messages=messages
    )

    # Get the response content from the model
    response_content = response_text(completion)
    
    # Save the result to the output text file
    with open(output_text_file, 'w') as f:
        f.write(response_content)
    
    print(f"Quality control report saved to {output_text_file}")

    return response_content


def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Generate a QC report for sequencing data.")
    
    # Arguments
    parser.add_argument('--file_path', required=True, help='Path to the NanoStats text file')
    parser.add_argument('--output_text_file', required=True, help='Path to save the output QC report')
    parser.add_argument('--species', default="Human", help='The species of the sample (default: Human)')
    parser.add_argument('--tissue', default="Skin", help='The tissue type of the sample (default: Skin)')
    add_llm_arguments(parser)

    # Parse the arguments
    args = parser.parse_args()

    try:
        client = client_from_args(args)
    except ValueError as exc:
        parser.error(str(exc))

    # Parse the NanoStats file and extract metrics
    metrics = parse_nanostats_file(args.file_path)

    # Keep this optional reporting step from aborting the biological workflow
    # when the configured provider is temporarily unavailable.
    try:
        generate_qc_report(
            metrics,
            client=client,
            species=args.species,
            tissue=args.tissue,
            model=args.model,
            output_text_file=args.output_text_file,
        )
    except Exception as exc:
        error_report = (
            "Failed to generate the read QC report via the configured LLM API.\n"
            f"Error: {exc}\n"
        )
        with open(args.output_text_file, "w") as handle:
            handle.write(error_report)
        print(error_report, file=sys.stderr)


if __name__ == '__main__':
    main()
