#!/usr/bin/env python3
"""Check whether the forward-strand percentage of the requested R1 reads falls inside the mixed range.

The input is the cmsearch tblout of the first N R1 reads. The reads must be in a FASTA file with
the record number as the ID (1, 2, 3, ...), so that repeated FASTQ names count independently.
"""

import argparse
import csv
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path


def count_reads(fasta):
    """Count the FASTA records that were searched."""
    with fasta.open() as reads:
        return sum(1 for line in reads if line.startswith(">"))


def classify_reads(hits, tested, args):
    """Calculate the forward percentage using the requested read count."""
    best = defaultdict(dict)
    for line in hits:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split()
        if len(fields) < 17 or fields[9] not in ("+", "-"):
            raise ValueError(f"Unexpected cmsearch tblout row: {line.rstrip()}")
        if fields[16] != "!":  # Only included hits count.
            continue

        read_id, strand, score = int(fields[0]), fields[9], float(fields[14])
        if not 1 <= read_id <= tested:
            raise ValueError(f"Unexpected target ID in cmsearch output: {read_id}")
        scores = best[read_id]
        if strand not in scores or score > scores[strand]:
            scores[strand] = score

    forward = reverse = 0
    for scores in best.values():
        # A missing strand ranks below every reported bit score.
        plus = scores.get("+", float("-inf"))
        minus = scores.get("-", float("-inf"))
        if plus > minus:
            forward += 1
        elif minus > plus:
            reverse += 1

    resolved = forward + reverse
    forward_percent = Fraction(100 * forward, args.reads)
    if resolved:
        mixed = 100 - args.cutoff <= forward_percent <= args.cutoff
        classification = "mixed" if mixed else "not_mixed"
    else:
        classification = "no_resolved_hits" if best else "no_hits"

    return {
        "sample": args.sample,
        "cm_database": args.cm_database,
        "reads_requested": args.reads,
        "reads_tested": tested,
        "resolved": resolved,
        "forward": forward,
        "reverse": reverse,
        "no_hit": tested - len(best),
        "forward_pct_of_requested": f"{float(forward_percent):.3f}",
        "cutoff_pct": args.cutoff_pct,
        "classification": classification,
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--sample", required=True, help="Sample ID for the output TSV")
    parser.add_argument(
        "--fasta",
        required=True,
        type=Path,
        help="FASTA of the reads searched by cmsearch, with record numbers as IDs",
    )
    parser.add_argument(
        "--tblout",
        required=True,
        type=Path,
        help="cmsearch --tblout table for the FASTA (searched with --cut_ga)",
    )
    parser.add_argument(
        "--cm-database",
        required=True,
        help="Name of the covariance-model database, for the output TSV",
    )
    parser.add_argument(
        "--reads",
        type=int,
        default=1000,
        help="Number of R1 reads requested and denominator for the forward percentage (default: 1000)",
    )
    parser.add_argument(
        "--cutoff-pct",
        default="60",
        help="Upper bound of the inclusive forward percentage range (default: 60, which gives 40 to 60)",
    )
    args = parser.parse_args()
    try:
        args.cutoff = Fraction(args.cutoff_pct)
    except (ValueError, ZeroDivisionError):
        parser.error("--cutoff-pct must be a numeric percentage")
    if not 50 <= args.cutoff <= 100:
        parser.error("--cutoff-pct must be between 50 and 100")
    if args.reads < 1:
        parser.error("--reads must be positive")
    return args


def main():
    args = parse_args()
    tested = count_reads(args.fasta)
    with args.tblout.open() as hits:
        row = classify_reads(hits, tested, args)

    writer = csv.DictWriter(
        sys.stdout, fieldnames=row, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerow(row)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        sys.exit(f"Error: {error}")
