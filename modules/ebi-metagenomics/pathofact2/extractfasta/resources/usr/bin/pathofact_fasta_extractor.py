#!/usr/bin/env python

# Copyright 2024-2026 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import gzip
import logging
import os

from Bio import SeqIO

logger = logging.getLogger(__name__)


def _open_text_maybe_gzip(path):
    """
    Return a text-mode file handle for plain or gzipped files.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


class BlastHit:
    """Selected best hit for a query."""

    __slots__ = ("sseqid", "evalue", "bitscore")

    def __init__(self, sseqid, evalue, bitscore):
        self.sseqid = sseqid
        self.evalue = evalue
        self.bitscore = bitscore


def _file_is_readable(path):
    """Return True if path exists and is non-empty."""
    if not path:
        return False
    if not os.path.exists(path):
        logger.warning("File not found: %s", path)
        return False
    if os.path.getsize(path) == 0:
        logger.info("Skipping empty file: %s", path)
        return False
    return True


def read_fasta_to_dict(fasta_path):
    """
    Read a FASTA file (optionally .gz) and return a dict of {sequence_id: SeqRecord}.
    """
    sequences = {}
    if not _file_is_readable(fasta_path):
        return sequences

    handle = _open_text_maybe_gzip(fasta_path)
    for record in SeqIO.parse(handle, "fasta"):
        sequences[record.id] = record
    handle.close()

    logger.info("Loaded %s sequences from FASTA: %s", len(sequences), fasta_path)
    return sequences


def parse_blastp_best_hits(blastp_out):
    """
    Parse DIAMOND blastp format 6 output with columns:
    qseqid sseqid pident length qlen slen evalue bitscore

    Select best hit per query:
      - prefer highest bitscore
      - then lowest evalue

    Store: {qseqid: BlastHit}
    """
    best = {}
    if not _file_is_readable(blastp_out):
        return best

    handle = _open_text_maybe_gzip(blastp_out)
    for line in handle:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        fields = line.split()
        if len(fields) < 8:
            continue

        qseqid = fields[0]
        sseqid = fields[1].split("(")[0]
        evalue = float(fields[6])
        bitscore = float(fields[7])

        candidate = BlastHit(sseqid=sseqid, evalue=evalue, bitscore=bitscore)
        current = best.get(qseqid)

        if current is None:
            best[qseqid] = candidate
            continue

        if candidate.bitscore > current.bitscore:
            best[qseqid] = candidate
        elif (
            candidate.bitscore == current.bitscore and candidate.evalue < current.evalue
        ):
            best[qseqid] = candidate

    logger.info(
        "Selected best BLASTP hits for %s queries from: %s", len(best), blastp_out
    )
    return best


def parse_pathofact2_predictions_tsv(tsv_path):
    """
    Parse PathoFact2 TSV output already filtered by threshold.
    Header: Sequence	Prediction	Probability

    Store: {sequence_id: probability}
    """
    preds = {}
    if not _file_is_readable(tsv_path):
        return preds

    handle = _open_text_maybe_gzip(tsv_path)
    # Read header line (file is non-empty by _file_is_readable)
    header = next(handle).rstrip("\n").split("\t")
    if len(header) < 3 or header[0] != "Sequence":
        logger.warning("Unexpected header in %s: %s", tsv_path, header)

    for line in handle:
        line = line.strip()
        if not line:
            continue

        fields = line.split("\t")
        if len(fields) < 3:
            continue

        seq_id = fields[0]
        probability = float(fields[2])
        preds[seq_id] = probability

    logger.info("Loaded %s PathoFact2 predictions from: %s", len(preds), tsv_path)
    return preds


def collect_detected_sequence_ids(blast_hits, tox_preds, vf_preds):
    """Union of all sequence IDs detected by any method."""
    ids = set()
    ids.update(blast_hits.keys())
    ids.update(tox_preds.keys())
    ids.update(vf_preds.keys())
    return ids


def write_detected_fasta(sequences, detected_ids, output_fasta):
    """
    Write a FASTA containing only sequences whose IDs are in detected_ids.
    """
    records = []
    missing = 0

    for seq_id in sorted(detected_ids):
        rec = sequences.get(seq_id)
        if rec is None:
            missing += 1
            continue
        records.append(rec)

    SeqIO.write(records, output_fasta, "fasta")
    logger.info("Wrote %s sequences to: %s", len(records), output_fasta)

    if missing > 0:
        logger.warning(
            "%s detected IDs were not found in FASTA and were skipped.", missing
        )


def write_support_table(blast_hits, tox_preds, vf_preds, output_tsv):
    """
    Write support TSV with header:
    sequence_id    detection_method    support_value_type    support_value

    detection_method:
      - blastp
      - pathofact2_tox
      - pathofact2_vf

    support_value_type:
      - blastp -> evalue
      - pathofact2_tox / pathofact2_vf -> probability
    """
    header = (
        "sequence_id\tdetection_method\tsupport_value_type\tsupport_value\tvfdb_hit\n"
    )

    rows = [header]

    for seq_id in sorted(blast_hits.keys()):
        hit = blast_hits[seq_id]
        rows.append(f"{seq_id}\tblastp\tevalue\t{hit.evalue}\t{hit.sseqid}\n")

    for seq_id in sorted(tox_preds.keys()):
        rows.append(f"{seq_id}\tpathofact2_tox\tprobability\t{tox_preds[seq_id]}\t\n")

    for seq_id in sorted(vf_preds.keys()):
        rows.append(f"{seq_id}\tpathofact2_vf\tprobability\t{vf_preds[seq_id]}\t\n")

    with open(output_tsv, "w") as out:
        out.writelines(rows)

    logger.info("Wrote support table with %s rows to: %s", len(rows) - 1, output_tsv)


def main():
    parser = argparse.ArgumentParser(
        description="Extract PathoFact2 candidate proteins from FASTA and report support from BLASTP and PathoFact2."
    )
    parser.add_argument(
        "-f", "--fasta", required=True, help="Protein sequences FASTA (optionally .gz)"
    )
    parser.add_argument(
        "-b", "--blastp_out", required=True, help="DIAMOND blastp tabular output"
    )
    parser.add_argument(
        "-t", "--pathofact2_tox", required=True, help="PathoFact2 toxins TSV (filtered)"
    )
    parser.add_argument(
        "-v", "--pathofact2_vf", required=True, help="PathoFact2 VF TSV (filtered)"
    )
    parser.add_argument("-o", "--output_prefix", required=True, help="Output prefix")
    parser.add_argument(
        "--verbose", action="store_true", default=False, help="Enable DEBUG logging"
    )
    args = parser.parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler()],
    )

    logger.info("Reading FASTA")
    seqs = read_fasta_to_dict(args.fasta)

    logger.info("Parsing BLASTP output (best hit per query)")
    blast_hits = parse_blastp_best_hits(args.blastp_out)

    logger.info("Parsing PathoFact2 toxin predictions")
    tox_preds = parse_pathofact2_predictions_tsv(args.pathofact2_tox)

    logger.info("Parsing PathoFact2 VF predictions")
    vf_preds = parse_pathofact2_predictions_tsv(args.pathofact2_vf)

    detected_ids = collect_detected_sequence_ids(blast_hits, tox_preds, vf_preds)
    logger.info("Total detected sequence IDs (union of methods): %s", len(detected_ids))

    out_fasta = f"{args.output_prefix}_pathofact2.fasta"
    out_tsv = f"{args.output_prefix}_support.tsv"

    logger.info("Writing outputs")
    write_detected_fasta(seqs, detected_ids, out_fasta)
    write_support_table(blast_hits, tox_preds, vf_preds, out_tsv)

    logger.info("Done")


if __name__ == "__main__":
    main()
