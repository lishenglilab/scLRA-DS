#!/usr/bin/env python3
import argparse
import gzip
import os
import sys
from pathlib import Path


def infer_sample_id_from_path(p: str) -> str:
    """
    Infer sample_id from a FASTQ path.

    Priority:
      1) Filename prefix before common tokens (e.g., AD1_matched_reads.fastq.gz -> AD1)
      2) Filename without extensions
      3) Parent directory name as a fallback
    """
    name = Path(p).name

    tokens = [
        "_matched_reads",
        ".matched_reads",
        ".hq",
        "_hq",
        ".fastq",
        ".fq",
    ]
    for t in tokens:
        if t in name:
            return name.split(t)[0]

    # Strip common suffixes
    for suf in [".fastq.gz", ".fq.gz", ".fastq", ".fq", ".gz"]:
        if name.endswith(suf):
            name = name[: -len(suf)]
            break

    # If still composite, use prefix before first underscore
    if "_" in name:
        return name.split("_")[0]

    return Path(p).parent.name


def iter_read_groups(fastq_gz: str, sample_id: str):
    opener = gzip.open if fastq_gz.endswith(".gz") else open
    with opener(fastq_gz, "rt") as f:
        record_number = 0
        while True:
            raw_header = f.readline()
            if raw_header == "":
                break
            record_number += 1
            header = raw_header.rstrip("\r\n")

            # FASTQ 4-line record: header, seq, '+', qual
            raw_sequence = f.readline()
            raw_separator = f.readline()
            raw_quality = f.readline()
            if not raw_sequence or not raw_separator or not raw_quality:
                raise ValueError(f"Incomplete FASTQ record after header: {header}")

            if not header.startswith("@"):
                raise ValueError(
                    f"FASTQ record {record_number} has an invalid header: {header}"
                )
            separator = raw_separator.rstrip("\r\n")
            if not separator.startswith("+"):
                raise ValueError(
                    f"FASTQ record {record_number} has an invalid separator line"
                )
            sequence = raw_sequence.rstrip("\r\n")
            quality = raw_quality.rstrip("\r\n")
            if len(sequence) != len(quality):
                raise ValueError(
                    f"FASTQ record {record_number} has sequence/quality length mismatch"
                )

            read_id = header[1:].split()[0]
            cell_barcode = read_id.split("_")[0] if "_" in read_id else read_id
            yield f"{read_id}\t{sample_id}_{cell_barcode}\n"


def main():
    ap = argparse.ArgumentParser(
        description="Generate IsoQuant --read_group TSV from one or more gzipped FASTQ files."
    )
    ap.add_argument(
        "--fastq",
        nargs="+",
        required=True,
        help="Input FASTQ(.gz) files (gzip-compressed is recommended).",
    )
    ap.add_argument(
        "--out",
        default="all_read_groups.tsv",
        help="Output TSV path (default: all_read_groups.tsv).",
    )
    ap.add_argument(
        "--sample-id",
        nargs="*",
        help="Optional sample IDs, same length and order as --fastq. "
             "If omitted, sample IDs are inferred from paths.",
    )
    args = ap.parse_args()

    fastqs = args.fastq
    outp = args.out

    sample_ids = args.sample_id
    if sample_ids is not None and len(sample_ids) > 0:
        if len(sample_ids) != len(fastqs):
            raise SystemExit(
                f"--sample-id length ({len(sample_ids)}) must equal --fastq length ({len(fastqs)})"
            )
    else:
        sample_ids = [infer_sample_id_from_path(p) for p in fastqs]

    missing_fastqs = [fq for fq in fastqs if not os.path.isfile(fq)]
    if missing_fastqs:
        raise SystemExit("Missing input FASTQ files: " + ", ".join(missing_fastqs))

    invalid_sample_ids = [
        sid for sid in sample_ids
        if not str(sid).strip() or any(char in str(sid) for char in "\t\r\n")
    ]
    if invalid_sample_ids:
        raise SystemExit("Sample IDs must be non-empty and cannot contain tabs/newlines")

    n_files = 0
    n_rows = 0

    with open(outp, "w") as out:
        for fq, sid in zip(fastqs, sample_ids):
            n_files += 1
            print(f"Processing {fq} as sample_id={sid}", file=sys.stderr)

            for line in iter_read_groups(fq, sid):
                out.write(line)
                n_rows += 1

    if n_files == 0:
        raise SystemExit("No input FASTQ files exist; nothing written.")
    if n_rows == 0:
        raise SystemExit("Input FASTQ files contained no read records; nothing written.")

    print(
        f"Done. Wrote {n_rows} read_group rows from {n_files} FASTQ files to {outp}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
