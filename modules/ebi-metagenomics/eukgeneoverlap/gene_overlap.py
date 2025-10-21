#!/usr/bin/env python
import argparse
from collections import defaultdict
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Compare predicted proteins (transcripts) between BRAKER and MetaEuk GFF files using exon/CDS overlap."
    )
    parser.add_argument('--gff_file1', type=str, required=True, help="BRAKER GFF file")
    parser.add_argument('--gff_file2', type=str, required=True, help="MetaEuk GFF file")
    parser.add_argument('--fasta1', type=str, required=False, help="BRAKER FASTA file (proteins)")
    parser.add_argument('--fasta2', type=str, required=False, help="MetaEuk FASTA file (proteins)")
    parser.add_argument('--output', type=str, required=True, help="Output prefix (sample name)")
    parser.add_argument('--threshold', type=float, default=0.9, help="Reciprocal exon/CDS overlap threshold")
    return parser.parse_args()



def read_cds_coords(gff_file, source_type):
    """
    Return {contig: {gene_id: [(start1,end1), (start2,end2), ...]}} for all CDS/exon coordinates.
    """
    genes = defaultdict(dict)
    gene_count = 0

    with open(gff_file, 'r') as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip().split("\t")
            if len(cols) != 9:
                continue
            contig, source, feature, start, end, score, strand, phase, attr = cols
            if feature not in ("CDS", "exon"):
                continue

            # Determine gene ID
            if source_type == "metaeuk":
                gene_id = attr.split(';')[0].split('=')[1]
            else:  # braker
                gene_id = [x.split('=')[1] for x in attr.split(';') if x.startswith("Parent=")][0]

            if gene_id not in genes[contig]:
                genes[contig][gene_id] = []
                gene_count += 1
            genes[contig][gene_id] = {"strand": strand, "exons": []}
            genes[contig][gene_id]["exons"].append((int(start), int(end)))
    
    # Sort coordinates per gene
    for contig in genes:
        for gene_id in genes[contig]:
            genes[contig][gene_id]["exons"] = sorted(
                genes[contig][gene_id]["exons"], key=lambda x: x[0]
        )

    print(f"Read {gene_count} genes with CDS/exons from {gff_file} ({source_type})")
    return genes


def exon_overlap(exons1, exons2, threshold=0.9):
    """
    exons1, exons2: list of (start, end) tuples
    Calculate total overlapping bases between two sets of exons.
    Returns True if reciprocal coverage ≥ threshold.
    """
    # sum of the exon lengths per gene
    total_length1 = sum(e[1]-e[0]+1 for e in exons1)
    total_length2 = sum(e[1]-e[0]+1 for e in exons2)

    # Compute total overlap
    ov = 0
    i, j = 0, 0
    while i < len(exons1) and j < len(exons2):
        a1, a2 = exons1[i]
        b1, b2 = exons2[j]
        start = max(a1, b1)
        end = min(a2, b2)
        if start <= end:
            ov += end - start + 1
        if a2 < b2:
            i += 1
        else:
            j += 1

    cov1 = ov / total_length1
    cov2 = ov / total_length2
    return cov1 >= threshold and cov2 >= threshold, ov, cov1, cov2


def parse_fasta(file_name, source_type):
    """
    Return {gene_id: SeqRecord} dictionary.
    """
    records = {}
    for rec in SeqIO.parse(file_name, "fasta"):
        header = rec.id
        if source_type == "metaeuk":
            gene_id = header.split('|')[0]  # take only Target_ID part
        else:  # braker
            gene_id = header  # remove .t1, .t2, etc.
        records[gene_id] = rec
    return records


def compare_genes(gff1_coords, gff2_coords, fasta1=None, fasta2=None, output="output", threshold=0.9):
    overlap = []
    uniq_1 = []
    uniq_2 = []
    matched_gff_1 = set()
    matched_gff_2 = set()

    # Compare genes per contig
    for contig in set(gff1_coords.keys()) & set(gff2_coords.keys()):
        for gene1, data1 in gff1_coords[contig].items():
            for gene2, data2 in gff2_coords[contig].items():
                # only look for overlaps if predictions are on the same strand
                if data1["strand"] != data2["strand"]:
                    continue
                is_same, ov, cov1, cov2 = exon_overlap(data1["exons"], data2["exons"], threshold)
                if is_same:
                    overlap.append([contig, gene1, gene2])
                    matched_gff_1.add((contig, gene1))
                    matched_gff_2.add((contig, gene2))

    # Unique genes
    for contig in gff1_coords:
        for gene1 in gff1_coords[contig]:
            if (contig, gene1) not in matched_gff_1:
                uniq_1.append([contig, gene1])
    for contig in gff2_coords:
        for gene2 in gff2_coords[contig]:
            if (contig, gene2) not in matched_gff_2:
                uniq_2.append([contig, gene2])

    # Write text outputs
    with open(output + "_overlap.tsv", "w") as out:
        out.write("contig\tbraker_gene\tmetaeuk_gene\n")
        for row in overlap:
            out.write("\t".join(row) + "\n")

    with open(output + "_braker_unique.tsv", "w") as out:
        out.write("contig\tbraker_gene\n")
        for row in uniq_1:
            out.write("\t".join(row) + "\n")

    with open(output + "_metaeuk_unique.tsv", "w") as out:
        out.write("contig\tmetaeuk_gene\n")
        for row in uniq_2:
            out.write("\t".join(row) + "\n")

    # Handle FASTA outputs
    if fasta1 and fasta2:
        braker_records = parse_fasta(fasta1, "braker")
        metaeuk_records = parse_fasta(fasta2, "metaeuk")

        overlap_records = []
        braker_unique_records = []
        metaeuk_unique_records = []

        added_overlap = set()
        added_braker_unique = set()
        added_metaeuk_unique = set()

        # Overlaps
        for contig, g1, g2 in overlap:
            if g1 in braker_records and g1 not in added_overlap:
                rec = braker_records[g1]
                rec.id = f"{rec.id}|overlap"
                rec.description = ""
                overlap_records.append(rec)
                added_overlap.add(g1)
            if g2 in metaeuk_records and g2 not in added_overlap:
                rec = metaeuk_records[g2]
                rec.id = f"{rec.id}|overlap"
                rec.description = ""
                overlap_records.append(rec)
                added_overlap.add(g2)

        # Unique BRAKER
        for contig, g1 in uniq_1:
            if g1 in braker_records and g1 not in added_braker_unique:
                rec = braker_records[g1]
                rec.id = f"{rec.id}|braker_only"
                rec.description = ""
                braker_unique_records.append(rec)
                added_braker_unique.add(g1)

        # Unique MetaEuk
        for contig, g2 in uniq_2:
            if g2 in metaeuk_records and g2 not in added_metaeuk_unique:
                rec = metaeuk_records[g2]
                rec.id = f"{rec.id}|metaeuk_only"
                rec.description = ""
                metaeuk_unique_records.append(rec)
                added_metaeuk_unique.add(g2)

        # Write FASTA outputs
        SeqIO.write(overlap_records, output + "_overlap.faa", "fasta")
        SeqIO.write(braker_unique_records, output + "_braker_unique.faa", "fasta")
        SeqIO.write(metaeuk_unique_records, output + "_metaeuk_unique.faa", "fasta")
        print("FASTA outputs written.")


def main():
    args = parse_args()
    print(f"Analyzing {args.output}...")

    gff1_coords = read_cds_coords(args.gff_file1, "braker")
    gff2_coords = read_cds_coords(args.gff_file2, "metaeuk")

    compare_genes(gff1_coords, gff2_coords, args.fasta1, args.fasta2, args.output, args.threshold)
    print("Comparison done!")


if __name__ == "__main__":
    main()
