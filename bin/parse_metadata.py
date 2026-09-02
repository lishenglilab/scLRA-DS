#!/usr/bin/env python3
import csv
import argparse
import re
import sys
from collections import Counter

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True, help="metadata CSV/TSV with header")
    ap.add_argument("--out", default="samples.tsv", help="output tsv path")
    args = ap.parse_args()

    meta = args.metadata
    outp = args.out

    # auto-detect delimiter by first line
    with open(meta, "r", newline="") as f:
        head = f.readline()
    delim = "\t" if "\t" in head else ","

    with open(meta, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        if reader.fieldnames is None:
            raise SystemExit("Metadata appears empty or missing header row")

        fields = [c.strip() for c in reader.fieldnames]
        reader.fieldnames = fields
        need = {"sample_id", "fastq"}
        if not need.issubset(set(fields)):
            raise SystemExit(f"Metadata must contain columns: sample_id, fastq. Got: {fields}")

        rows = []
        for line_number, row in enumerate(reader, start=2):
            sid = (row.get("sample_id") or "").strip()
            fq = (row.get("fastq") or "").strip()
            if sid and fq:
                if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", sid):
                    raise SystemExit(
                        f"Invalid sample_id on line {line_number}: {sid!r}. "
                        "Use only letters, numbers, '.', '_' and '-'."
                    )
                rows.append((sid, fq))

    if not rows:
        raise SystemExit("No valid rows found in metadata (need non-empty sample_id and fastq).")

    sample_ids = [sid for sid, _ in rows]
    duplicate_ids = sorted(sid for sid, count in Counter(sample_ids).items() if count > 1)
    if duplicate_ids:
        raise SystemExit(
            "Metadata sample_id values must be unique. Duplicates: "
            + ", ".join(duplicate_ids)
        )

    with open(outp, "w", newline="") as out:
        out.write("sample_id\tfastq\n")
        for sid, fq in rows:
            out.write(f"{sid}\t{fq}\n")

    print(f"Wrote {len(rows)} samples to {outp}", file=sys.stderr)

if __name__ == "__main__":
    main()
