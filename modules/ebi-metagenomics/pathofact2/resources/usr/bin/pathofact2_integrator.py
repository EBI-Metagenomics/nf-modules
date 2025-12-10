#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2024-2025 EMBL - European Bioinformatics Institute
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
import fileinput
import logging
import os
from collections import defaultdict
from typing import Dict, List, Optional, Union, DefaultDict

# Set up logger
logger = logging.getLogger(__name__)


def validate_inputs(optional_inputs: Dict[str, Optional[str]]) -> List[str]:
    """
    Validate that input files exist and are not empty or header-only.

    Args:
        optional_inputs: Dictionary mapping input names to file paths (can be None)

    Returns:
        List of valid input names whose files should be parsed.
    """
    valid_inputs: List[str] = []

    for name, path in optional_inputs.items():
        if not path:
            continue

        # Check existence
        if not os.path.exists(path):
            logger.warning(f"File not found for '{name}' → {path}")
            continue

        # Check for empty files
        if os.path.getsize(path) == 0:
            logger.info(f"Skipping empty file for '{name}' → {path}")
            continue

        # Read first few lines using hook_compressed (handles .gz automatically)
        with fileinput.hook_compressed(path, "r", encoding="utf-8", errors="ignore") as input_table:
            lines = [line for _, line in zip(range(3), input_table)]  # read up to 3 lines

        if len(lines) == 0:
            logger.info(f"Skipping empty file for '{name}' → {path}")
            continue
        elif len(lines) == 1:
            logger.info(f"Skipping header-only file for '{name}' → {path}")
            continue

        # If passed all checks
        valid_inputs.append(name)

    return valid_inputs


def parse_pathofact_tsv(
    pathofact2_annotation: DefaultDict[str, Dict[str, Dict[str, Union[List[str], str]]]], input_file: str
) -> DefaultDict[str, Dict[str, Dict[str, Union[List[str], str]]]]:
    """
    Parse Pathofact2 annotation files filtered in tsv format

    Args:
        pathofact2_annotation: Nested defaultdict to store Pathofact2 annotations
        input_file: Path to the tsv input file
                    This file has 3 columns:
                    Sequence	Prediction	Probability

    Returns:
        Updated pathofact2_annotation defaultdict
    """
    with fileinput.hook_compressed(input_file, "r", encoding="utf-8") as input_table:
        # Skipping the header line
        next(input_table)
        for line in input_table:
            line_l: List[str] = line.rstrip().split("\t")
            protein_id: str = line_l[0]
            probability: int = line_l[2]

            pathofact2_annotation[protein_id][analysis_software_name] = {"drug_class": drug_class_list, "seq_identity": seq_identity}

    return pathofact2_annotation



def parse_annotation_dict(pathofact2_annotation: DefaultDict[str, Dict[str, Dict[str, Union[List[str], str]]]]) -> Dict[str, List[str]]:
    """
    Parse Pathofact2 annotation dictionary and format protein attributes.

    Args:
        pathofact2_annotation: Nested defaultdict containing pathofact2 annotations

    Returns:
        Dictionary mapping protein IDs to lists of formatted attribute strings
    """
    protein_attributes: Dict[str, List[str]] = {}

    for protein, tools in pathofact2_annotation.items():
        all_drugs: List[str] = []
        tool_names: List[str] = []
        tool_identities: List[str] = []

        for tool, info in tools.items():
            # Collect drug classes
            drug_classes: Union[List[str], str] = info.get("drug_class", [])
            if isinstance(drug_classes, list):
                all_drugs.extend(drug_classes)
            else:
                all_drugs.append(drug_classes)

            # Collect tool name and identity value
            tool_names.append(tool)
            seq_identity: Union[List[str], str] = info.get("seq_identity", "NA")
            if isinstance(seq_identity, str):
                tool_identities.append(seq_identity)
            else:
                tool_identities.append("NA")

        # Clean up and deduplicate drug classes
        unique_drugs: List[str] = list(dict.fromkeys(all_drugs))
        unique_drugs = [item.replace(" ", "_") for item in unique_drugs]

        # Build attribute strings
        drugs_string: str = "drug_class=" + ",".join(unique_drugs)
        tools_string: str = "amr_tool=" + ",".join(tool_names)
        idents_string: str = "amr_tool_ident=" + ",".join(tool_identities)

        # Append to dictionary
        protein_attributes[protein] = [drugs_string, tools_string, idents_string]

    return protein_attributes


def parse_gff(cds_gff: str, output_file: str, protein_attributes: Dict[str, List[str]]) -> None:
    """
    Parse GFF file and integrate Pathofact2 annotations into output GFF.

    Args:
        cds_gff: Path to input GFF file containing CDS coordinates
        output_file: Path to output GFF file
        protein_attributes: Dictionary mapping protein IDs to attribute lists
    """
    with (
        fileinput.hook_compressed(cds_gff, "r", encoding="utf-8") as input_table,
        open(output_file, "w") as output_gff,
    ):
        output_gff.write("##gff-version 3\n")
        for line in input_table:
            line = line.rstrip()
            line_l: List[str] = line.split("\t")
            # Annotation lines have exactly 9 columns
            if len(line_l) == 9:
                (
                    contig,
                    seq_source,
                    seq_type,
                    start,
                    end,
                    score,
                    strand,
                    phase,
                    attr,
                ) = line.rstrip().split("\t")
                if seq_type == "CDS":
                    features_list: List[str] = attr.split(";")
                    feature_id: str = features_list[0].split("=")[1]
                    if feature_id in protein_attributes:
                        attr = attr.rstrip(";")
                        new_attribute: str = attr + ";" + ";".join(protein_attributes[feature_id])
                        line_l.pop(-1)
                        line_l.append(new_attribute)
                        to_print: str = "\t".join(line_l)
                        output_gff.write(to_print + "\n")


def main() -> None:
    """
    Main function to orchestrate Pathofact2 results
    """
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description="Integration of toxins and virulence factors predicted and annotated using Pathofact2 into a single gff3 file"
    )
    parser.add_argument("-t", "--toxins", dest="toxins_tsv", help="Result of Pathofact2 tool for toxins prediction after filtering", required=False, default=None)
    parser.add_argument("-r", "--virulence", dest="virulence_tsv", help="Result of Pathofact2 tool for virulence factors after filtering", required=False, default=None)
    parser.add_argument(
        "-c", "--cds_gff", dest="cds_gff", help="GFF file containing the coordinates of the CDSs used for Pathofact2 annotation", required=True, default=None
    )
    parser.add_argument("-o", "--output", dest="output", help="Name of the output file", required=False, default="integrated_result.gff")
    parser.add_argument("-v", "--verbose", dest="verbose", help="Enable verbose logging (DEBUG level)", action="store_true", default=False)
    args: argparse.Namespace = parser.parse_args()

    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.StreamHandler(),  # Console output
        ],
    )

    optional_inputs: Dict[str, Optional[str]] = {
        "toxins": args.toxins_tsv,
        "virulence": args.virulence_tsv,
    }

    valid_inputs: List[str] = validate_inputs(optional_inputs)
    if valid_inputs:
        pathofact2_annotation: DefaultDict[str, Dict[str, Dict[str, Union[List[str], str]]]] = defaultdict(dict)
        logger.info(f"The valid inputs provided for integration are: {', '.join(valid_inputs)}")
        logger.info("Parsing valid inputs")

        for name in valid_inputs:
            logger.debug(f"Processing {name} input file: {optional_inputs[name]}")
                pathofact2_annotation = parse_pathofact_tsv(pathofact2_annotation, optional_inputs[name])

        logger.info("Parsing the annotation information per protein")
        protein_attributes: Dict[str, List[str]] = parse_annotation_dict(pathofact2_annotation)
        logger.debug(f"Found Pathofact2 annotations for {len(protein_attributes)} proteins")

        logger.info("Parsing gff file and writing output file")
        parse_gff(args.cds_gff, args.output, protein_attributes)
        logger.info(f"Output written to: {args.output}")

    else:
        logger.info("No inputs to parse, the output file is not generated.")


if __name__ == "__main__":
    main()
