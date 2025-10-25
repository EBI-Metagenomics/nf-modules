#!/usr/bin/env python
import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Compare predicted proteins (transcripts) between BRAKER and MetaEuk GFF files using exon/CDS overlap."
    )
    parser.add_argument('--fasta_aa', type=str, required=True, help="BRAKER protein fasta")
    parser.add_argument('--fasta_fn', type=str, required=True, help="BRAKER nucleotide fasta")
    parser.add_argument('--output', type=str, required=True, help="Output prefix (sample name)")
    return parser.parse_args()


def rename_fasta_ids(fasta_file, output, suffix):
    '''
    Rename BRAKER protein fasta IDs to append '|braker_only' suffix
    '''
    records = SeqIO.parse(fasta_file, "fasta")

    for rec in records:
        new_id = f"{rec.id}|braker_only"
        rec.id = new_id

    SeqIO.write(records, output + suffix, "fasta" )

def main():
    args = parse_args()
    print(f"Renaming {args.output}...")

    rename_fasta_ids(args.fast_aa, args.output, "_braker_unique.faa")
    rename_fasta_ids(args.fasta_fn, args.output, "_braker_unique.ffn")

    print("Renaming done")


if __name__ == "__main__":
    main()

