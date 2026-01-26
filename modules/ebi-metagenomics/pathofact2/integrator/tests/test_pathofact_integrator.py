#!/usr/bin/env python

import sys
from pathlib import Path

import pytest

# Import functions from the script
from pathofact2_integrator import (
    parse_cdd,
    parse_gff,
    parse_ips,
    parse_pathofact_support,
    validate_inputs,
)

# Add the script directory to the path
script_dir = Path(__file__).parent.parent / "resources" / "usr" / "bin"
sys.path.insert(0, str(script_dir))


@pytest.fixture
def tmp_files(tmp_path):
    """Create temporary test files with realistic data."""

    # PathOFact2 support file data
    pred_support_data = {
        "header": [
            "sequence_id",
            "detection_method",
            "support_value_type",
            "support_value",
            "vfdb_hit",
        ],
        "rows": [
            ["WP_000031783.1", "blastp", "evalue", "6.5e-248", "VFG046459"],
            ["WP_000242755.1", "blastp", "evalue", "2.58e-94", "VFG042734"],
            ["WP_001273183.1", "blastp", "evalue", "0.0", "VFG049144"],
            ["WP_000242755.1", "pathofact2_tox", "probability", "1.0", ""],
            ["WP_001273183.1", "pathofact2_tox", "probability", "1.0", ""],
            ["WP_000163771.1", "pathofact2_vf", "probability", "1.0", ""],
            ["WP_001273183.1", "pathofact2_vf", "probability", "0.9977883", ""],
        ],
    }

    # CDD annotation file data
    cdd_annot_data = {
        "header": [
            "query",
            "hit_type",
            "pssm_id",
            "from",
            "to",
            "evalue",
            "bitscore",
            "accession",
            "short_name",
            "incomplete",
            "superfamily_pssm_id",
        ],
        "rows": [
            [
                "WP_000031783.1",
                "Specific",
                "206671",
                "11",
                "203",
                "4.85688e-133",
                "377.694",
                "cd01884",
                "EF_Tu",
                "-",
                "476819",
            ],
            [
                "WP_000031783.1",
                "Specific",
                "294006",
                "300",
                "389",
                "1.02934e-65",
                "202.357",
                "cd03707",
                "EFTU_III",
                "-",
                "470675",
            ],
            [
                "WP_000031783.1",
                "Specific",
                "293898",
                "211",
                "297",
                "1.80452e-52",
                "168.466",
                "cd03697",
                "EFTU_II",
                "-",
                "445922",
            ],
            [
                "WP_000242755.1",
                "Specific",
                "237999",
                "8",
                "117",
                "3.94612e-35",
                "118.968",
                "cd00038",
                "CAP_ED",
                "-",
                "469590",
            ],
            [
                "WP_000242755.1",
                "Specific",
                "238044",
                "140",
                "207",
                "1.57469e-15",
                "66.9205",
                "cd00092",
                "HTH_CRP",
                "-",
                "481199",
            ],
        ],
    }

    # InterProScan annotation file data (no header)
    ips_annot_data = {
        "rows": [
            [
                "WP_000242755.1",
                "9a58754a5c561c0784b75360a962f090",
                "210",
                "CDD",
                "cd00092",
                "HTH_CRP",
                "140",
                "207",
                "1.14133E-15",
                "T",
                "15-01-2026",
                "-",
                "-",
                "-",
                "-",
            ],
            [
                "WP_000242755.1",
                "9a58754a5c561c0784b75360a962f090",
                "210",
                "CDD",
                "cd00038",
                "CAP_ED",
                "8",
                "117",
                "2.86014E-35",
                "T",
                "15-01-2026",
                "IPR000595",
                "Cyclic nucleotide-binding domain",
                "-",
                "Reactome:R-BTA-1296072",
            ],
            [
                "WP_000031783.1",
                "3ba5ec971eda5dc5f5c1c093a5df4e11",
                "394",
                "CDD",
                "cd03697",
                "EFTU_II",
                "211",
                "297",
                "1.30792E-52",
                "T",
                "15-01-2026",
                "IPR033720",
                "Elongation factor Tu, domain 2",
                "-",
                "Reactome:R-HSA-5389840",
            ],
            [
                "WP_000031783.1",
                "3ba5ec971eda5dc5f5c1c093a5df4e11",
                "394",
                "CDD",
                "cd01884",
                "EF_Tu",
                "11",
                "203",
                "0.0",
                "T",
                "15-01-2026",
                "IPR041709",
                "Elongation factor Tu (EF-Tu), GTP-binding domain",
                "-",
                "Reactome:R-HSA-5389840",
            ],
            [
                "WP_000031783.1",
                "3ba5ec971eda5dc5f5c1c093a5df4e11",
                "394",
                "CDD",
                "cd03707",
                "EFTU_III",
                "300",
                "389",
                "7.46066E-66",
                "T",
                "15-01-2026",
                "-",
                "-",
                "-",
                "-",
            ],
        ]
    }

    # GFF file data
    gff_data = {
        "header": [
            "##gff-version 3",
            "##sequence-region Contig001 1 100000",
            "##sequence-region Contig002 1 80000",
        ],
        "rows": [
            [
                "Contig001",
                "Prodigal",
                "CDS",
                "51000",
                "51837",
                ".",
                "+",
                "0",
                "ID=WP_000259031.1;Name=WP_000259031.1;product=hypothetical protein",
            ],
            [
                "Contig001",
                "Prodigal",
                "CDS",
                "87000",
                "87630",
                ".",
                "+",
                "0",
                "ID=WP_000242755.1;Name=WP_000242755.1;product=hypothetical protein",
            ],
            [
                "Contig002",
                "Prodigal",
                "CDS",
                "3000",
                "3744",
                ".",
                "+",
                "0",
                "ID=WP_000163771.1;Name=WP_000163771.1;product=hypothetical protein",
            ],
            [
                "Contig002",
                "Prodigal",
                "CDS",
                "38000",
                "39182",
                ".",
                "-",
                "0",
                "ID=WP_000031783.1;Name=WP_000031783.1;product=hypothetical protein",
            ],
            [
                "Contig003",
                "Prodigal",
                "CDS",
                "12000",
                "15102",
                ".",
                "+",
                "0",
                "ID=WP_001273183.1;Name=WP_001273183.1;product=hypothetical protein",
            ],
        ],
    }

    # Create pred_support.tsv
    pred_support = tmp_path / "pred_support.tsv"
    pred_support_content = "\t".join(pred_support_data["header"]) + "\n"
    pred_support_content += "\n".join(
        ["\t".join(row) for row in pred_support_data["rows"]]
    )
    pred_support.write_text(pred_support_content + "\n")

    # Create CDD annotation file
    cdd_annot = tmp_path / "cdd_annot.tsv"
    cdd_annot_content = "\t".join(cdd_annot_data["header"]) + "\n"
    cdd_annot_content += "\n".join(["\t".join(row) for row in cdd_annot_data["rows"]])
    cdd_annot.write_text(cdd_annot_content + "\n")

    # Create IPS annotation file (no header)
    ips_annot = tmp_path / "ips_annot.tsv"
    ips_annot_content = "\n".join(["\t".join(row) for row in ips_annot_data["rows"]])
    ips_annot.write_text(ips_annot_content + "\n")

    # Create GFF file
    gff_file = tmp_path / "input.gff"
    gff_content = "\n".join(gff_data["header"]) + "\n"
    gff_content += "\n".join(["\t".join(row) for row in gff_data["rows"]])
    gff_file.write_text(gff_content + "\n")

    # Create empty file
    empty_file = tmp_path / "empty.tsv"
    empty_file.write_text("")

    # Create header-only file
    header_only = tmp_path / "header_only.tsv"
    header_only.write_text("sequence_id\tdetection_method\tsupport_value_type\n")

    return {
        "pred_support": str(pred_support),
        "cdd_annot": str(cdd_annot),
        "ips_annot": str(ips_annot),
        "gff_file": str(gff_file),
        "empty_file": str(empty_file),
        "header_only": str(header_only),
        "tmp_path": tmp_path,
    }


class TestValidateInputs:
    """Test input validation function."""

    def test_valid_inputs(self, tmp_files):
        """Test validation passes for valid files."""
        to_validate = {
            "gff": tmp_files["gff_file"],
            "annotation": tmp_files["cdd_annot"],
        }
        result = validate_inputs(to_validate)
        assert len(result) == 0, "Valid files should pass validation"

    def test_missing_file(self, tmp_files):
        """Test validation fails for non-existent files."""
        to_validate = {
            "gff": "/nonexistent/file.gff",
        }
        result = validate_inputs(to_validate)
        assert len(result) == 1
        assert "File not found" in result["gff"]

    def test_empty_file(self, tmp_files):
        """Test validation fails for empty files."""
        to_validate = {
            "annotation": tmp_files["empty_file"],
        }
        result = validate_inputs(to_validate)
        assert len(result) == 1
        assert result["annotation"] == "Empty file"

    def test_header_only_file(self, tmp_files):
        """Test validation fails for header-only files."""
        to_validate = {
            "support": tmp_files["header_only"],
        }
        result = validate_inputs(to_validate)
        assert len(result) == 1
        assert result["support"] == "Header only"

    def test_none_input(self):
        """Test validation handles None inputs."""
        to_validate = {
            "optional_file": None,
        }
        result = validate_inputs(to_validate)
        assert len(result) == 1
        assert result["optional_file"] == "No input provided"


class TestParsePathofactSupport:
    """Test parsing of PathOFact2 support file."""

    def test_parse_support_file(self, tmp_files):
        """Test parsing support file with mixed detection methods."""
        result = parse_pathofact_support(tmp_files["pred_support"])

        # Check all proteins are present
        assert "WP_000031783.1" in result
        assert "WP_000242755.1" in result
        assert "WP_001273183.1" in result
        assert "WP_000163771.1" in result

        # Check blastp annotation format
        assert "vfdb=VFG046459" in result["WP_000031783.1"]
        assert "blastp_eval=6.5e-248" in result["WP_000031783.1"]

        # Check pathofact2 annotation format
        assert "pathofact2_tox_prob=1.0" in result["WP_000242755.1"]
        assert "pathofact2_vf_prob=1.0" in result["WP_000163771.1"]

        # Check multiple annotations are concatenated
        assert "vfdb=VFG042734" in result["WP_000242755.1"]
        assert "pathofact2_tox_prob=1.0" in result["WP_000242755.1"]

        # WP_001273183.1 has all three: blastp, tox, vf
        wp_001273183_annot = result["WP_001273183.1"]
        assert "vfdb=VFG049144" in wp_001273183_annot
        assert "pathofact2_tox_prob=1.0" in wp_001273183_annot
        assert "pathofact2_vf_prob=0.9977883" in wp_001273183_annot


class TestParseCDD:
    """Test parsing of CDD annotation file."""

    def test_parse_cdd_basic(self, tmp_files):
        """Test basic CDD parsing."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        result = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])

        # Check proteins with both pathofact and CDD annotations
        assert "WP_000031783.1" in result
        assert "WP_000242755.1" in result

        # Check CDD format: accession:short_name
        assert "cd01884:EF_Tu" in result["WP_000031783.1"]
        assert "cd00038:CAP_ED" in result["WP_000242755.1"]
        assert "cd00092:HTH_CRP" in result["WP_000242755.1"]

    def test_parse_cdd_multiple_hits(self, tmp_files):
        """Test CDD with multiple hits per protein."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        result = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])

        # WP_000031783.1 has 3 CDD hits
        wp_031783_annot = result["WP_000031783.1"]
        assert "cd01884:EF_Tu" in wp_031783_annot
        assert "cd03707:EFTU_III" in wp_031783_annot
        assert "cd03697:EFTU_II" in wp_031783_annot

        # Check they're comma-separated
        assert wp_031783_annot.count(",") == 2

    def test_parse_cdd_preserves_pathofact(self, tmp_files):
        """Test that CDD parsing preserves PathOFact annotations."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        result = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])

        # Check pathofact annotations are still present
        assert "vfdb=VFG046459" in result["WP_000031783.1"]
        assert "blastp_eval=6.5e-248" in result["WP_000031783.1"]
        assert ";cdd=" in result["WP_000031783.1"]

    def test_parse_cdd_only_annotates_pathofact_proteins(self, tmp_files):
        """Test that only proteins in pathofact_attrs get CDD annotations."""
        # Create minimal pathofact_attrs with only one protein
        pathofact_attrs = {"WP_000242755.1": "vfdb=VFG042734;blastp_eval=2.58e-94"}
        result = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])

        # Only WP_000242755.1 should be in result
        assert "WP_000242755.1" in result
        assert "WP_000031783.1" not in result


class TestParseIPS:
    """Test parsing of InterProScan annotation file."""

    def test_parse_ips_basic(self, tmp_files):
        """Test basic IPS parsing."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        result = parse_ips(pathofact_attrs, tmp_files["ips_annot"])

        # Check proteins with both pathofact and IPS CDD annotations
        assert "WP_000031783.1" in result
        assert "WP_000242755.1" in result

    def test_parse_ips_filters_cdd_only(self, tmp_files):
        """Test that IPS parser only extracts CDD annotations."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        result = parse_ips(pathofact_attrs, tmp_files["ips_annot"])

        # Check CDD annotations are present
        assert "cd00092:HTH_CRP" in result["WP_000242755.1"]
        assert "cd00038:CAP_ED" in result["WP_000242755.1"]

        # Check format matches: accession:description
        assert "cd03697:EFTU_II" in result["WP_000031783.1"]
        assert "cd01884:EF_Tu" in result["WP_000031783.1"]
        assert "cd03707:EFTU_III" in result["WP_000031783.1"]

    def test_parse_ips_multiple_hits(self, tmp_files):
        """Test IPS with multiple CDD hits per protein."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        result = parse_ips(pathofact_attrs, tmp_files["ips_annot"])

        # WP_000031783.1 has 3 CDD annotations
        wp_031783_annot = result["WP_000031783.1"]
        assert wp_031783_annot.count(",") == 2  # comma-separated

        # WP_000242755.1 has 2 CDD annotations
        wp_000242755_annot = result["WP_000242755.1"]
        assert wp_000242755_annot.count(",") == 1


class TestParseGFF:
    """Test GFF parsing and annotation integration."""

    def test_parse_gff_integration(self, tmp_files):
        """Test GFF parsing integrates annotations correctly."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        proteins_annotation = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])

        output_file = tmp_files["tmp_path"] / "output.gff"
        parse_gff(tmp_files["gff_file"], str(output_file), proteins_annotation)

        # Read output
        content = output_file.read_text()
        lines = content.strip().split("\n")

        # Check header
        assert lines[0] == "##gff-version 3"

        # Check annotated lines contain pathofact + CDD annotations
        annotated_lines = [line for line in lines if "WP_000242755.1" in line]
        assert len(annotated_lines) == 1
        assert "vfdb=VFG042734" in annotated_lines[0]
        assert "cd00038:CAP_ED" in annotated_lines[0]
        assert "cd00092:HTH_CRP" in annotated_lines[0]

    def test_parse_gff_only_annotated_proteins(self, tmp_files):
        """Test that only proteins with annotations are in output."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        proteins_annotation = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])

        output_file = tmp_files["tmp_path"] / "output.gff"
        parse_gff(tmp_files["gff_file"], str(output_file), proteins_annotation)

        content = output_file.read_text()

        # Proteins with annotations should be present
        assert "WP_000242755.1" in content
        assert "WP_000031783.1" in content

        # Protein without PathOFact prediction should NOT be in output
        assert "WP_000259031.1" not in content

    def test_parse_gff_preserves_coordinates(self, tmp_files):
        """Test that GFF coordinates are preserved correctly."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        proteins_annotation = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])

        output_file = tmp_files["tmp_path"] / "output.gff"
        parse_gff(tmp_files["gff_file"], str(output_file), proteins_annotation)

        content = output_file.read_text()
        lines = [line for line in content.split("\n") if "WP_000242755.1" in line]

        assert len(lines) == 1
        fields = lines[0].split("\t")

        # Check coordinates match input
        assert fields[0] == "Contig001"
        assert fields[3] == "87000"
        assert fields[4] == "87630"
        assert fields[6] == "+"

    def test_parse_gff_attribute_format(self, tmp_files):
        """Test that attributes are formatted correctly."""
        pathofact_attrs = {
            "WP_000242755.1": "vfdb=VFG042734;blastp_eval=2.58e-94;cdd=cd00038:CAP_ED"
        }

        output_file = tmp_files["tmp_path"] / "output.gff"
        parse_gff(tmp_files["gff_file"], str(output_file), pathofact_attrs)

        content = output_file.read_text()
        lines = [line for line in content.split("\n") if "WP_000242755.1" in line]

        assert len(lines) == 1
        attributes = lines[0].split("\t")[8]

        # Check attribute starts with ID
        assert attributes.startswith("ID=WP_000242755.1;")

        # Check semicolon-separated format
        assert "vfdb=VFG042734" in attributes
        assert "blastp_eval=2.58e-94" in attributes
        assert "cdd=cd00038:CAP_ED" in attributes


class TestIntegration:
    """Integration tests for complete workflow."""

    def test_full_workflow_cdd(self, tmp_files):
        """Test complete workflow with CDD annotations."""
        # Step 1: Parse PathOFact support
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        assert len(pathofact_attrs) == 4

        # Step 2: Parse CDD annotations
        proteins_annotation = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])
        assert len(proteins_annotation) == 2  # Only 2 proteins have CDD hits

        # Step 3: Parse GFF and write output
        output_file = tmp_files["tmp_path"] / "final_output.gff"
        parse_gff(tmp_files["gff_file"], str(output_file), proteins_annotation)

        # Verify output
        content = output_file.read_text()
        assert "##gff-version 3" in content
        assert "WP_000242755.1" in content
        assert "WP_000031783.1" in content
        assert "vfdb=" in content
        assert "cdd=" in content

    def test_full_workflow_ips(self, tmp_files):
        """Test complete workflow with IPS annotations."""
        # Step 1: Parse PathOFact support
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])

        # Step 2: Parse IPS annotations
        proteins_annotation = parse_ips(pathofact_attrs, tmp_files["ips_annot"])
        assert len(proteins_annotation) == 2

        # Step 3: Parse GFF and write output
        output_file = tmp_files["tmp_path"] / "final_output_ips.gff"
        parse_gff(tmp_files["gff_file"], str(output_file), proteins_annotation)

        # Verify output
        content = output_file.read_text()
        assert "##gff-version 3" in content
        assert "WP_000242755.1" in content
        assert "WP_000031783.1" in content


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
