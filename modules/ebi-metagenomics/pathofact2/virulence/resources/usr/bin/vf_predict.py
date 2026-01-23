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
import math
import random
import re
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

# import pickle
from functools import reduce

import joblib
import numpy as np
import pandas as pd

# import time
# import psutil

# process = psutil.Process()
# Memory_usage_before = process.memory_info().rss / 1024 ** 2
# start_time = time.time()


def readfasta(file):
    myfasta = []
    name = None
    seq_list = []

    with open(file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    # Join sequence list into a single string and append
                    myfasta.append([name, "".join(seq_list)])
                name = line.split()[0][1:]
                seq_list = []  # Reset sequence list
            else:
                # Use re.sub directly, but append lines instead of using +=
                seq_list.append(re.sub("[^ARNDCQEGHILKMFPSTWYV-]", "-", line.upper()))

        # Add the last sequence
        if name is not None:
            myfasta.append([name, "".join(seq_list)])

    return myfasta


def aac(fastas):
    aa_string = "ACDEFGHIKLMNPQRSTVWY"
    encodings = []
    header = ["#"]
    for i in aa_string:
        header.append(i)
    encodings.append(header)
    for i in fastas:
        name, sequence = i[0], re.sub("-", "", i[1])
        count = Counter(sequence)
        for key in count:
            count[key] = count[key] / len(sequence)
        code = [name]
        for aa in aa_string:
            code.append(count[aa])
        encodings.append(code)
    df = pd.DataFrame(encodings[1:], columns=encodings[0])
    return df


def count_occurrences(seq1, sequence):
    return sum(sequence.count(aa) for aa in seq1)


def count2(aa_set, sequence):
    indices = [i + 1 for i, aa in enumerate(sequence) if aa in aa_set]
    length = len(sequence)
    number = len(indices)
    cutoff_nums = [
        1,
        math.floor(0.25 * number),
        math.floor(0.50 * number),
        math.floor(0.75 * number),
        number,
    ]
    cutoff_nums = [max(i, 1) for i in cutoff_nums]
    return [
        (indices[cutoff - 1] / length) * 100 if cutoff <= number else 0
        for cutoff in cutoff_nums
    ]


def ctdc_for_sequence(name, sequence, group1, group2, group3, properties):
    sequence = sequence.replace("-", "")
    encodings = []
    seq_len = len(sequence)
    for prop in properties:
        c1 = count_occurrences(group1[prop], sequence) / seq_len
        c2 = count_occurrences(group2[prop], sequence) / seq_len
        c3 = 1 - c1 - c2
        encodings.extend([c1, c2, c3])
    return [name] + encodings


def ctd_d_for_sequence(name, sequence, group1, group2, group3, properties):
    sequence = sequence.replace("-", "")
    encodings = []
    for prop in properties:
        encodings.extend(
            count2(group1[prop], sequence)
            + count2(group2[prop], sequence)
            + count2(group3[prop], sequence)
        )
    return [name] + encodings


def ctd_t_for_sequence(name, sequence, group1, group2, group3, properties):
    sequence = sequence.replace("-", "")
    aa_pairs = [sequence[i : i + 2] for i in range(len(sequence) - 1)]
    encodings = []
    for prop in properties:
        c1221 = sum(
            1
            for pair in aa_pairs
            if (pair[0] in group1[prop] and pair[1] in group2[prop])
            or (pair[0] in group2[prop] and pair[1] in group1[prop])
        )
        c1331 = sum(
            1
            for pair in aa_pairs
            if (pair[0] in group1[prop] and pair[1] in group3[prop])
            or (pair[0] in group3[prop] and pair[1] in group1[prop])
        )
        c2332 = sum(
            1
            for pair in aa_pairs
            if (pair[0] in group2[prop] and pair[1] in group3[prop])
            or (pair[0] in group3[prop] and pair[1] in group2[prop])
        )
        total_pairs = len(aa_pairs)
        encodings.extend(
            [c1221 / total_pairs, c1331 / total_pairs, c2332 / total_pairs]
        )
    return [name] + encodings


def process_sequence(name, sequence, aadict):
    sequence = sequence.replace("-", "")
    code = [name]
    tmp_code = np.zeros(400)
    for j in range(len(sequence) - 1):
        tmp_code[aadict[sequence[j]] * 20 + aadict[sequence[j + 1]]] += 1
    total = tmp_code.sum()
    if total > 0:
        tmp_code /= total
    code.extend(tmp_code)
    return code


def dpc(fastas, num_cpus):
    aa_string = "ACDEFGHIKLMNPQRSTVWY"
    di_peptides = [aa1 + aa2 for aa1 in aa_string for aa2 in aa_string]
    aadict = {aa: i for i, aa in enumerate(aa_string)}
    encodings = []
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        future_to_sequence = {
            executor.submit(process_sequence, name, sequence, aadict): (name, sequence)
            for name, sequence in fastas
        }
        for future in as_completed(future_to_sequence):
            encodings.append(future.result())
    df = pd.DataFrame(encodings, columns=["#"] + di_peptides)
    return df


def ctdc(fastas, group1, group2, group3, properties, num_cpus):
    encodings = []
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        future_to_sequence = {
            executor.submit(
                ctdc_for_sequence, name, sequence, group1, group2, group3, properties
            ): (name, sequence)
            for name, sequence in fastas
        }
        for future in as_completed(future_to_sequence):
            encodings.append(future.result())
    headers = ["#"] + [f"{p}.G{g}" for p in properties for g in range(1, 4)]
    df = pd.DataFrame(encodings, columns=headers)
    return df


def ctd_d(fastas, group1, group2, group3, properties, num_cpus):
    encodings = []
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        future_to_sequence = {
            executor.submit(
                ctd_d_for_sequence, name, sequence, group1, group2, group3, properties
            ): (name, sequence)
            for name, sequence in fastas
        }
        for future in as_completed(future_to_sequence):
            encodings.append(future.result())
    headers = ["#"] + [
        f"{p}.{g}.residue{d}"
        for p in properties
        for g in ("1", "2", "3")
        for d in ["0", "25", "50", "75", "100"]
    ]
    df = pd.DataFrame(encodings, columns=headers)
    return df


def ctd_t(fastas, group1, group2, group3, properties, num_cpus):
    encodings = []
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        future_to_sequence = {
            executor.submit(
                ctd_t_for_sequence, name, sequence, group1, group2, group3, properties
            ): (name, sequence)
            for name, sequence in fastas
        }
        for future in as_completed(future_to_sequence):
            encodings.append(future.result())
    headers = ["#"] + [
        f"{p}.{tr}" for p in properties for tr in ("Tr1221", "Tr1331", "Tr2332")
    ]
    df = pd.DataFrame(encodings, columns=headers)
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process fasta input with random forest virulence prediction model"
    )
    parser.add_argument("--file", help="Input fasta file", required=True)
    parser.add_argument(
        "--model", metavar="model", help="The model file (SAV format)", required=True
    )
    parser.add_argument(
        "--cpus", metavar="cpus", type=int, help="number of CPUs", default=1
    )
    parser.add_argument(
        "--outfile",
        metavar="outfile",
        help="The output file (TSV format)",
        required=True,
    )
    args = parser.parse_args()

    group1 = {
        "hydrophobicity_PRAM900101": "RKEDQN",
        "hydrophobicity_ARGP820101": "QSTNGDE",
        "hydrophobicity_ZIMJ680101": "QNGSWTDERA",
        "hydrophobicity_PONP930101": "KPDESNQT",
        "hydrophobicity_CASG920101": "KDEQPSRNTG",
        "hydrophobicity_ENGD860101": "RDKENQHYP",
        "hydrophobicity_FASG890101": "KERSQD",
        "normwaalsvolume": "GASTPDC",
        "polarity": "LIFWCMVY",
        "polarizability": "GASDT",
        "charge": "KR",
        "secondarystruct": "EALMQKRH",
        "solventaccess": "ALFCGIVW",
    }

    group2 = {
        "hydrophobicity_PRAM900101": "GASTPHY",
        "hydrophobicity_ARGP820101": "RAHCKMV",
        "hydrophobicity_ZIMJ680101": "HMCKV",
        "hydrophobicity_PONP930101": "GRHA",
        "hydrophobicity_CASG920101": "AHYMLV",
        "hydrophobicity_ENGD860101": "SGTAW",
        "hydrophobicity_FASG890101": "NTPG",
        "normwaalsvolume": "NVEQIL",
        "polarity": "PATGS",
        "polarizability": "CPNVEQIL",
        "charge": "ANCQGHILMFPSTWYV",
        "secondarystruct": "VIYCWFT",
        "solventaccess": "RKQEND",
    }

    group3 = {
        "hydrophobicity_PRAM900101": "CLVIMFW",
        "hydrophobicity_ARGP820101": "LYPFIW",
        "hydrophobicity_ZIMJ680101": "LPFYI",
        "hydrophobicity_PONP930101": "YMFWLCVI",
        "hydrophobicity_CASG920101": "FIWC",
        "hydrophobicity_ENGD860101": "CVLIMF",
        "hydrophobicity_FASG890101": "AYHWVMFLIC",
        "normwaalsvolume": "MHKFRYW",
        "polarity": "HQRKNED",
        "polarizability": "KMHFRYW",
        "charge": "DE",
        "secondarystruct": "GNPSD",
        "solventaccess": "MSPTHY",
    }

    properties = (
        "hydrophobicity_PRAM900101",
        "hydrophobicity_ARGP820101",
        "hydrophobicity_ZIMJ680101",
        "hydrophobicity_PONP930101",
        "hydrophobicity_CASG920101",
        "hydrophobicity_ENGD860101",
        "hydrophobicity_FASG890101",
        "normwaalsvolume",
        "polarity",
        "polarizability",
        "charge",
        "secondarystruct",
        "solventaccess",
    )

    random.seed(1234)
    fastas = readfasta(args.file)
    encoding1 = aac(fastas)
    encoding2 = dpc(fastas, args.cpus)
    encoding3 = ctdc(fastas, group1, group2, group3, properties, args.cpus)
    encoding4 = ctd_t(fastas, group1, group2, group3, properties, args.cpus)
    encoding5 = ctd_d(fastas, group1, group2, group3, properties, args.cpus)

    data_frames = [encoding1, encoding2, encoding3, encoding4, encoding5]
    df = reduce(
        lambda left, right: pd.merge(left, right, on=["#"], how="outer"), data_frames
    )

    X = df.drop("#", axis=1)
    rfc = joblib.load(open(args.model, "rb"))
    rfc_predict = rfc.predict(X)
    rfc_probs = rfc.predict_proba(X)[:, 1]

    prediction = pd.DataFrame(
        {"Sequence": df["#"], "Prediction": rfc_predict, "Probability": rfc_probs}
    )
    prediction.to_csv(args.outfile, sep="\t", index=False)


# end_time = time.time()
# print(f"Execution time: {end_time - start_time} seconds")
# print(f"Memory usage: {(process.memory_info().rss / 1024 ** 2)-Memory_usage_before} MB")
