#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# Copyright 2024 EMBL - European Bioinformatics Institute
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

# Load necessary packages
library(dada2)

args = commandArgs(trailingOnly=TRUE)  # Expects at most three arguments: prefix, forward fastq, and reverse fastq (optional)
prefix = args[1]  # Prefix
path_f = args[2]  # Forward fastq
path_r = ifelse(length(args) > 2, args[3], NA)  # Reverse fastq, if provided

# Learn error model
err_f = learnErrors(path_f, multithread=TRUE)
if (!is.na(path_r)) {
    err_r = learnErrors(path_r, multithread=TRUE)
}

# Dereplicate sequences
drp_f = derepFastq(path_f)
if (!is.na(path_r)) {
    drp_r = derepFastq(path_r)
}

# Generate stranded ASVs
dada_f = dada(drp_f, err=err_f, multithread=TRUE)
if (!is.na(path_r)) {
    dada_r = dada(drp_r, err=err_r, multithread=TRUE)
}

# Merge stranded ASVs
if (!is.na(path_r)) {
    merged = mergePairs(dada_f, drp_f, dada_r, drp_r, verbose=TRUE)
} else {
    merged = dada_f
}

if (length(merged$sequence) == 0) {
    print("No ASVs - stopping script early.")
    quit(status=1)
} else {
    # Make ASV count table
    seqtab = makeSequenceTable(merged)

    # Remove chimeras
    seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

    # Save ASV count table
    write.table(seqtab.nochim, file = paste0("./", prefix, "_asv_counts.tsv"), sep = "\t", row.names=FALSE)

    # Save ASV sequences to FASTA file
    seqtab.length = length(seqtab.nochim)
    unqs = getUniques(seqtab.nochim)
    id_list = paste("seq", 1:length(unqs), sep="_")
    uniquesToFasta(unqs, paste0("./", prefix, "_asvs.fasta"), id_list)
}
