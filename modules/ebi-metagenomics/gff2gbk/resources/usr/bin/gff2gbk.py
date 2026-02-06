#!/usr/bin/env python

import argparse
import gzip
import logging
import os

from BCBio import GFF
from Bio import SeqIO

logger = logging.getLogger(__name__)


def _open_text_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _file_is_readable(path):
    return path and os.path.exists(path) and os.path.getsize(path) > 0


# -----------------------------
# FASTA loaders
# -----------------------------


def read_fasta_dict(fasta_path):
    records = {}
    with _open_text_maybe_gzip(fasta_path) as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            records[rec.id] = rec
    return records


# -----------------------------
# Translation helpers
# -----------------------------


def _first_qualifier(feature, keys):
    for k in keys:
        if k in feature.qualifiers:
            v = feature.qualifiers[k]
            if isinstance(v, list):
                return v[0]
            return v
    return None


def _safe_feature_id(feature):
    return _first_qualifier(feature, ["ID", "locus_tag", "Name"])


def _tag_provenance(feature):
    src = _first_qualifier(feature, ["source"])
    if not src:
        return

    note = feature.qualifiers.get("note", [])
    if isinstance(note, str):
        note = [note]

    note.append(f"gene_caller={src}")
    feature.qualifiers["note"] = note


# -----------------------------
# Translation injection
# -----------------------------


def add_translations(records, protein_dict=None):
    added = 0
    missing = 0

    for rec in records:
        for feat in rec.features:
            if feat.type != "CDS":
                continue

            _tag_provenance(feat)

            cds_id = _safe_feature_id(feat)

            if not cds_id:
                missing += 1
                continue

            # Prefer protein FASTA
            if protein_dict and cds_id in protein_dict:
                prot = str(protein_dict[cds_id].seq)
                feat.qualifiers["translation"] = [prot]
                added += 1
                continue

            # Fallback translation
            nuc = feat.extract(rec.seq)
            nuc = nuc[: (len(nuc) // 3) * 3]

            if len(nuc) == 0:
                missing += 1
                continue

            prot = str(nuc.translate(table=11))
            feat.qualifiers["translation"] = [prot]
            added += 1

    logger.info("Translations added: %s | missing CDS: %s", added, missing)


# -----------------------------
# Main conversion
# -----------------------------


def convert_to_gbk(contigs, gff, output_gbk, proteins=None):
    contigs_dict = read_fasta_dict(contigs)

    protein_dict = None
    if proteins:
        logger.info("Loading protein FASTA")
        protein_dict = read_fasta_dict(proteins)

    records = []

    with _open_text_maybe_gzip(gff) as handle:
        for rec in GFF.parse(handle, base_dict=contigs_dict):
            if "molecule_type" not in rec.annotations:
                rec.annotations["molecule_type"] = "DNA"

            records.append(rec)

    add_translations(records, protein_dict)

    with open(output_gbk, "w") as out:
        SeqIO.write(records, out, "genbank")

    logger.info("GBK written to %s", output_gbk)


# -----------------------------
# CLI
# -----------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Convert MGnify combined GFF + contigs + optional protein FASTA into GBK"
    )

    parser.add_argument("--gff", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--proteins", required=False)
    parser.add_argument("--output_gbk", required=True)

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    convert_to_gbk(
        contigs=args.contigs,
        gff=args.gff,
        output_gbk=args.output_gbk,
        proteins=args.proteins,
    )


if __name__ == "__main__":
    main()
