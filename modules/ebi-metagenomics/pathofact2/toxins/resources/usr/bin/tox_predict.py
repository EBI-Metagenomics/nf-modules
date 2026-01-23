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
# AUTHOR: Luis F. Delgado (lfdelzam@hotmail.com)

import argparse
import os

import joblib
from joblib import Parallel, delayed
from scipy.sparse import hstack


# Function to read sequences and their IDs from a FASTA file
def read_fasta_sequences(file_path):
    sequences = []
    with open(file_path) as fin:
        seq_id = None
        seq = ""
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
                if seq_id and seq:
                    sequences.append((seq_id, seq))
                seq_id = line.split()[0][1:]  # Remove the '>' character to get the ID
                seq = ""
            else:
                seq += line
        if seq_id and seq:
            sequences.append((seq_id, seq))  # Append the last sequence
    return sequences


# Function to preprocess sequences into k-mer features using loaded vectorizers
def preprocess_sequences(sequences, kmer_sizes, vectorizer_dir, num_cpus):
    def process_kmer(k):
        vectorizer_path = os.path.join(vectorizer_dir, f"vectorizer_{k}.pkl")
        vectorizer = joblib.load(vectorizer_path)  # Load the saved vectorizer
        seqs = [
            seq for _, seq in sequences
        ]  # Extract only the sequences for vectorization
        x_kmer = vectorizer.transform(seqs)
        return x_kmer

    # Process k-mer features in parallel
    x_kmers = Parallel(n_jobs=num_cpus)(delayed(process_kmer)(k) for k in kmer_sizes)
    x_combined = hstack(x_kmers)
    return x_combined


# Main function to handle argument parsing and workflow execution
def main():
    parser = argparse.ArgumentParser(
        description="Predict sequences using a pre-trained model."
    )
    parser.add_argument(
        "-m",
        "--model",
        type=str,
        required=True,
        help="Path to the saved model (final_model.joblib)",
    )
    parser.add_argument(
        "-s",
        "--sequences",
        type=str,
        required=True,
        help="Path to the FASTA file containing new sequences",
    )
    parser.add_argument(
        "-k",
        "--kmers",
        type=int,
        nargs="+",
        default=[5],
        help="List of k-mer sizes used during training",
    )
    parser.add_argument(
        "-v",
        "--vectorizer_dir",
        type=str,
        required=True,
        help="Directory containing the vectorizer files",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path to the output file to save predictions",
    )
    parser.add_argument(
        "--cut_off",
        type=str,
        default="0",
        help="probaility cut_off to be cosidered as positive classification",
    )
    parser.add_argument(
        "-cpu",
        "--cpus",
        type=int,
        default=-1,
        help="Number of CPUs to use for parallel processing. Default is -1 (use all available CPUs)",
    )

    args = parser.parse_args()

    num_cpus = (
        args.cpus if args.cpus > 0 else -1
    )  # Use all available CPUs if non-positive value is provided

    # Step 1: Load the saved model
    best_model = joblib.load(args.model)

    # Step 2: Load new sequences for prediction
    sequences_with_ids = read_fasta_sequences(args.sequences)

    # Step 3: Preprocess new sequences using the saved vectorizers
    x_new = preprocess_sequences(
        sequences_with_ids, args.kmers, args.vectorizer_dir, num_cpus
    )

    # print(best_model.get_params())
    # Step 4: Predict using the loaded model
    predictions = best_model.predict(x_new)
    probabilities = best_model.predict_proba(x_new)

    # Step 5: Output the predictions to a file
    with open(args.output, "w") as fout:
        fout.write("Sequence\tPrediction\tProbability\n")
        for i, (seq_id, _) in enumerate(sequences_with_ids):
            prob = probabilities[i][1]
            pred = predictions[i]
            if prob >= float(args.cut_off):
                fout.write(f"{seq_id}\t{pred}\t{prob:.2f}\n")


if __name__ == "__main__":
    main()
