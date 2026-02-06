from __future__ import annotations

import sys
from pathlib import Path

import pytest
from Bio import SeqIO

# Add the script directory to the path
script_dir = Path(__file__).parent.parent / "resources" / "usr" / "bin"
sys.path.insert(0, str(script_dir))

from gff2gbk import (  # noqa: E402
    convert_to_gbk,
)


def _write_fasta(path: Path, records: dict[str, str]) -> None:
    """
    Write a minimal FASTA file from a dict: {id: sequence}.
    """
    with path.open("w", encoding="utf-8") as fh:
        for rec_id, seq in records.items():
            fh.write(f">{rec_id}\n")
            # write as single line to keep snapshots stable
            fh.write(f"{seq}\n")


def _write_gff3(path: Path, rows: dict[str, dict[str, str]]) -> None:
    """
    Write a minimal GFF3 file from a dict-of-dicts.

    rows example:
      {
        "cds1": {"seqid":"ERZ...", "source":"Pyrodigal", "type":"CDS", "start":"1", "end":"90", "strand":"-", "phase":".", "id":"ERZ..._1"},
        ...
      }

    Keys are just for readability; ordering preserved by insertion order (Py3.7+).
    """
    with path.open("w", encoding="utf-8") as fh:
        fh.write("##gff-version 3\n")
        for _, r in rows.items():
            seqid = r["seqid"]
            source = r["source"]
            ftype = r["type"]
            start = r["start"]
            end = r["end"]
            score = r.get("score", ".")
            strand = r["strand"]
            phase = r.get("phase", ".")
            attributes = f"ID={r['id']}"
            fh.write(
                "\t".join(
                    [seqid, source, ftype, start, end, score, strand, phase, attributes]
                )
                + "\n"
            )


def _extract_cds_translations(gbk_path: Path) -> dict[str, str]:
    """
    Read GBK and return {CDS_ID: translation_string}.
    """
    out: dict[str, str] = {}
    for rec in SeqIO.parse(str(gbk_path), "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            cds_id = feat.qualifiers.get("ID", [None])[0]
            translation = feat.qualifiers.get("translation", [None])[0]
            if cds_id is None:
                raise AssertionError("CDS feature missing ID qualifier")
            if translation is None:
                raise AssertionError(f"CDS {cds_id} missing translation qualifier")
            out[str(cds_id)] = str(translation)
    return out


def _extract_cds_notes(gbk_path: Path) -> dict[str, str]:
    """
    Return {CDS_ID: note_string}.
    """
    out: dict[str, str] = {}
    for rec in SeqIO.parse(str(gbk_path), "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            cds_id = feat.qualifiers.get("ID", [None])[0]
            note = feat.qualifiers.get("note", [None])[0]
            if cds_id is None:
                raise AssertionError("CDS feature missing ID qualifier")
            if note is None:
                raise AssertionError(f"CDS {cds_id} missing note qualifier")
            out[str(cds_id)] = str(note)
    return out


@pytest.fixture()
def dummy_inputs(tmp_path: Path) -> dict[str, Path]:
    """
    Builds dummy contigs, GFF, and proteins as dicts and writes them to tmp_path.

    Returns paths as a dict.
    """
    contigs = {
        "ERZ1290913_80": (
            "CATCGTTGCACCCACCCCGGCAGCAATCCGCCGCCGTGGGTTCACACTCGGCTTTTCTGTACTTGTGATGCTGTCGAACTCGATGGTCATGCCAATCCCCTGTCGCTTTCCCTGTGCAGG"
            "CCCCGTCGTTTGCCCTCAGCTTGTCAGACGCTTGGACACCTCGCTCGGTTCCCTGATCGACACATTTTCCTGCCAAATTTCTTGCCACCGTTATCCGGTTCGTGCCACACCGGCCGTCGCG"
            "AGCGGTCGTCGGCTGAACAGGCCCGAAACTCAGGGTCGTTGACCGCCCAAGTCCCGGGTGGAGCTCGGTACGGTTTCTTTCATGGACACGACGAGCACGGCGAAGCAGATCGATCTTTACTG"
            "GCGGCCCGGTTGTGGCTTCTGCTCGAGCCTCATGCGTGGCCTCGACAAGCTTGGCGTCGAACGAGTTGAGCACAATATTTGGGACAACAAAGCTGATGCCGCCATCGTGCGCAAGCACGCGA"
            "ACGGAAACGAAGTCGTGCCCACCGTCGTCATCGGCGACAAGGGCCTCGTCAACCCATCTGCCGGTGCGCTCCTGGTATTTCTCGCCGAGAACGCACCCCATCTGCTCCCCGAAGGCGTCGAAG"
            "CACCCCAGCCTGGCAAGGTCGGTCGGTTCGTCGGTCGGGTCCTCGGAAACTGAATAACTCTCGCCGCCGTGGCGTGCCTAGTTGTCGTTCGCTGACCTGAGCCCCGCCGGGATGACGACTCAC"
            "CTTGTATTGTTCGCGCCATGACTGACACATCGCTTCCCGCTTGGGACGCCGACGGCATGGACATCGCCGAATCGATTCTTGACCTGATCGGCAACACGCCGATGGTGCGCATGAAACGCGTCA"
            "GTGAGGAACACGGCATCCGCTGCACCCTGGCGATGAAACATGAGGTCACGAACCCAGGTGGTTCCTCGAAAGACCGACCGGCACTCGAGATGATCCTCGCAGCCGAAGCAAGTGGCGAACTTC"
            "AGCCCGGCGGAACGATCGTCGAGCCCACGTCTGGCAACACCGGTGTCGGACTCGCCATCGTGGCTGCCCAGCGGGGCTACAAGTGCGTGTTCGTCATGACCAACAAGGCCGGCAAGGAAAAGGT"
            "CGACCTGCTCCGTGCCTACGGCGCCGAAGTCGTCGTCTGTGATGTCGCCGTCGCACCGGAAGATCCAAACAGCTACTACTCCGTGGCCGAGCGACTCACCCAGGAACGCGGCGCGTTTCGTCCG"
            "AATCAGTACGCCAACCCGCACAACCCGGCAGCCCACACGAAGACCACCGGCCCCGAAATTTGGGAGCAAACCAAGGGCCGCATCACGCACTTCATTGCGGGCGCGGGCACATGCGGTTCGCTCAC"
            "CGGTACTGGCCGCTACCTCAAGTCGATGAATCCGGACATCAAGATCATCGCTGCTGACCCCGAGAAGTCGGTGTTCTCTGGCGGTTCGGGTCGC"
        )
    }

    gff_rows = {
        "cds1": {
            "seqid": "ERZ1290913_80",
            "source": "Pyrodigal",
            "type": "CDS",
            "start": "1",
            "end": "90",
            "strand": "-",
            "phase": ".",
            "id": "ERZ1290913_80_1",
        },
        "cds2": {
            "seqid": "ERZ1290913_80",
            "source": "Pyrodigal",
            "type": "CDS",
            "start": "323",
            "end": "661",
            "strand": "+",
            "phase": ".",
            "id": "ERZ1290913_80_2",
        },
        "cds3": {
            "seqid": "ERZ1290913_80",
            "source": "Pyrodigal",
            "type": "CDS",
            "start": "749",
            "end": "1444",
            "strand": "+",
            "phase": ".",
            "id": "ERZ1290913_80_3",
        },
        "cds4": {
            "seqid": "ERZ1290913_80",
            "source": "FragGeneScanRS",
            "type": "CDS",
            "start": "2",
            "end": "283",
            "strand": "+",
            "phase": ".",
            "id": "ERZ1290913_80_2_283_+",
        },
    }

    proteins = {
        "ERZ1290913_80_1": "MTIEFDSITSTEKPSVNPRRRIAAGVGATM",
        "ERZ1290913_80_2": (
            "MDTTSTAKQIDLYWRPGCGFCSSLMRGLDKLGVERVEHNIWDNKADAAIVRKHANGNEVVPTVVIGDKGLVNPSAGALLVFLAENAPHLLPEGVEAPQPGKVGRFVGRVLGN"
        ),
        "ERZ1290913_80_2_283_+": (
            "IVAPTPAAIRRRGFTLGFSVLVMLSNSMVMPIPCRFPCAGPVVCPQLVRRLDTSLGSLIDTFSCQISCHRYPVRATPAVASGRRLNRPETQGR"
        ),
        "ERZ1290913_80_3": (
            "MTDTSLPAWDADGMDIAESILDLIGNTPMVRMKRVSEEHGIRCTLAMKHEVTNPGGSSKDRPALEMILAAEASGELQPGGTIVEPTSGNTGVGLAIVAAQRGYKCVFVMTNKAGKEKVDLLRAYGAEVVVCDVAVAPEDPNSYYSVAERLTQERGAFRPNQYANPHNPAAHTKTTGPEIWEQTKGRITHFIAGAGTCGSLTGTGRYLKSMNPDIKIIAADPEKSVFSGGSGR"
        ),
    }

    contigs_path = tmp_path / "contigs.fna"
    gff_path = tmp_path / "calls.gff"
    proteins_path = tmp_path / "proteins.faa"

    _write_fasta(contigs_path, contigs)
    _write_gff3(gff_path, gff_rows)
    _write_fasta(proteins_path, proteins)

    return {
        "contigs": contigs_path,
        "gff": gff_path,
        "proteins": proteins_path,
        "proteins_dict": proteins,  # for assertions
    }


@pytest.mark.parametrize("with_proteins", [True, False])
def test_convert_to_gbk_translations_and_provenance(
    tmp_path: Path, dummy_inputs: dict[str, Path], with_proteins: bool
):
    out_gbk = tmp_path / (
        "out_with_proteins.gbk" if with_proteins else "out_no_proteins.gbk"
    )

    proteins_path: Path | None = dummy_inputs["proteins"] if with_proteins else None

    convert_to_gbk(
        contigs=str(dummy_inputs["contigs"]),
        gff=str(dummy_inputs["gff"]),
        output_gbk=str(out_gbk),
        proteins=str(proteins_path) if proteins_path else None,
    )

    assert out_gbk.exists()
    translations = _extract_cds_translations(out_gbk)
    notes = _extract_cds_notes(out_gbk)

    # Check we have all 4 CDS IDs
    expected_ids = {
        "ERZ1290913_80_1",
        "ERZ1290913_80_2",
        "ERZ1290913_80_3",
        "ERZ1290913_80_2_283_+",
    }
    assert set(translations.keys()) == expected_ids
    assert set(notes.keys()) == expected_ids

    # Provenance notes must include the caller source
    assert "gene_caller=Pyrodigal" in notes["ERZ1290913_80_1"]
    assert "gene_caller=Pyrodigal" in notes["ERZ1290913_80_2"]
    assert "gene_caller=Pyrodigal" in notes["ERZ1290913_80_3"]
    assert "gene_caller=FragGeneScanRS" in notes["ERZ1290913_80_2_283_+"]

    if with_proteins:
        # When proteins are supplied, translation MUST exactly match provided FASTA (no added '*')
        prot = dummy_inputs["proteins_dict"]
        assert translations["ERZ1290913_80_1"] == prot["ERZ1290913_80_1"]
        assert translations["ERZ1290913_80_2"] == prot["ERZ1290913_80_2"]
        assert translations["ERZ1290913_80_3"] == prot["ERZ1290913_80_3"]
        assert translations["ERZ1290913_80_2_283_+"] == prot["ERZ1290913_80_2_283_+"]

        assert not translations["ERZ1290913_80_2"].endswith("*")
        assert not translations["ERZ1290913_80_2_283_+"].endswith("*")
    else:
        # Without proteins, translation is derived from nucleotide extraction; in your expected output:
        # ERZ1290913_80_2 ends with '*' and ERZ1290913_80_2_283_+ ends with '*'
        assert translations["ERZ1290913_80_2"].endswith("*")
        assert translations["ERZ1290913_80_2_283_+"].endswith("*")

        # Others should still be non-empty and look like proteins
        assert len(translations["ERZ1290913_80_1"]) > 0
        assert len(translations["ERZ1290913_80_3"]) > 0
