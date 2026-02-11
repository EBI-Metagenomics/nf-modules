from __future__ import annotations

import sys
from pathlib import Path

import pytest
from Bio import SeqIO

# Add the script directory to the path
script_dir = Path(__file__).parent.parent / "resources" / "usr" / "bin"
sys.path.insert(0, str(script_dir))

from gbk_generator import (  # noqa: E402
    build_records_from_gff,
    build_records_from_prodigal_faa,
    validate_prodigal_faa_headers,
    write_genbank,
)


# -------------------------------------------------------------------------
# Small file writers
# -------------------------------------------------------------------------
def _write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        for rid, seq in records.items():
            fh.write(f">{rid}\n{seq}\n")


def _write_faa_with_headers(path: Path, header_to_seq: dict[str, str]) -> None:
    """
    header_to_seq keys are the *full header line without leading '>'*.
    """
    with path.open("w", encoding="utf-8") as fh:
        for header, seq in header_to_seq.items():
            fh.write(f">{header}\n{seq}\n")


def _write_gff3(path: Path, rows: list[dict[str, str]]) -> None:
    """
    rows expects dicts with keys:
      seqid, source, type, start, end, strand, id
    Optionally: product, locus_tag
    """
    with path.open("w", encoding="utf-8") as fh:
        fh.write("##gff-version 3\n")
        for r in rows:
            attrs = [f"ID={r['id']}"]
            if "locus_tag" in r and r["locus_tag"]:
                attrs.append(f"locus_tag={r['locus_tag']}")
            if "product" in r and r["product"]:
                attrs.append(f"product={r['product']}")
            attr_s = ";".join(attrs)
            fh.write(
                "\t".join(
                    [
                        r["seqid"],
                        r["source"],
                        r.get("type", "CDS"),
                        str(r["start"]),
                        str(r["end"]),
                        ".",
                        r["strand"],
                        ".",
                        attr_s,
                    ]
                )
                + "\n"
            )


# -------------------------------------------------------------------------
# GenBank helpers
# -------------------------------------------------------------------------
def _load_gbk_features(gbk_path: Path) -> list:
    recs = list(SeqIO.parse(str(gbk_path), "genbank"))
    assert recs, "No GenBank records parsed"
    # single contig test data => one record, but keep generic
    feats = []
    for rec in recs:
        feats.extend(rec.features)
    return feats


def _cds_by_protein_id(gbk_path: Path) -> dict[str, dict[str, list[str]]]:
    """
    Return mapping: protein_id -> qualifiers dict (lists).
    """
    out: dict[str, dict[str, list[str]]] = {}
    for rec in SeqIO.parse(str(gbk_path), "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            pid = feat.qualifiers.get("protein_id", [None])[0]
            assert pid, "CDS missing protein_id qualifier"
            out[str(pid)] = feat.qualifiers
    return out


# -------------------------------------------------------------------------
# Prodigal header fixtures (good vs bad)
# -------------------------------------------------------------------------
def _good_prodigal_headers() -> dict[str, str]:
    """
    Well-formed headers based on your examples:
      <gene_id> # <start> # <end> # <strand> # <extra...>
    """
    return {
        # Canonical example
        "NC_000913_4 # 3734 # 5020 # 1 # ID=1_4;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-10bp;gc_cont=0.528": "MKT",
        # MGnify combined gene caller style (Pyrodigal)
        "ERZ1290913_1_2 # 61 # 4599 # -1 # ID=3_2;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.600": "MKT",
        "ERZ1290913_107_19 # 24433 # 28155 # 1 # ID=132_19;partial=00;start_type=GTG;rbs_motif=None;rbs_spacer=None;gc_cont=0.628": "MKT",
    }


def _bad_headers_like_cgc() -> dict[str, str]:
    """
    Badly-formed headers inspired by common CGC / non-Prodigal problems.
    These should fail parse_prodigal_header() in gbk_generator.py.
    """
    return {
        # 1) Missing the required " # " fields entirely
        "ERZ1290913_105_9309_9384_-": "MKT",
        # 2) Has separators but missing strand field (only 3 parts after split)
        "ERZ1290913_82_1 # 10 # 90 # ID=oops_missing_strand": "MKT",
        # 3) Strand not 1 or -1
        "ERZ1290913_116_2 # 100 # 200 # 0 # ID=bad_strand": "MKT",
        # 4) Non-integer coordinates
        "ERZ1290913_126_3 # start # 200 # 1 # ID=bad_coord": "MKT",
        # 5) Coordinates < 1
        "ERZ1290913_126_4 # 0 # 10 # 1 # ID=coord_lt_1": "MKT",
    }


# -------------------------------------------------------------------------
# Tests
# -------------------------------------------------------------------------
def test_validate_prodigal_headers_counts_bad(tmp_path: Path) -> None:
    faa = tmp_path / "calls.faa"
    headers = {}
    headers.update(_good_prodigal_headers())
    headers.update(_bad_headers_like_cgc())
    _write_faa_with_headers(faa, headers)

    total, bad = validate_prodigal_faa_headers(str(faa))
    assert total == len(headers)
    assert bad == len(_bad_headers_like_cgc())


def test_prodigal_mode_strict_headers_raises(tmp_path: Path) -> None:
    # Minimal contigs: only needed so build_records_from_prodigal_faa can proceed
    contigs = tmp_path / "contigs.fna"
    _write_fasta(contigs, {"NC_000913": "ATG" * 1000})

    faa = tmp_path / "calls.faa"
    headers = {}
    headers.update(_good_prodigal_headers())
    headers.update(_bad_headers_like_cgc())
    _write_faa_with_headers(faa, headers)

    with pytest.raises(
        ValueError, match=r"FAA records do not follow Prodigal header convention"
    ):
        build_records_from_prodigal_faa(
            contigs_path=str(contigs),
            faa_path=str(faa),
            proteins_path=None,
            prefix="test",
            default_product="hypothetical protein",
            locus_tag_prefix="",
            gene_from_locus_tag=False,
            skip_missing_contigs=True,
            require_translation=False,
            transl_table=11,
            drop_terminal_stop=True,
            strict_headers=True,
        )


def test_gff_mode_emits_source_feature_and_conventional_qualifiers(
    tmp_path: Path,
) -> None:
    contigs = tmp_path / "contigs.fna"
    gff = tmp_path / "calls.gff"
    prots = tmp_path / "proteins.faa"
    out_gbk = tmp_path / "out.gbk"

    # One contig with enough length
    _write_fasta(contigs, {"ERZ1290913_80": "A" * 2000})

    # Two CDS, two different sources, with products
    _write_gff3(
        gff,
        rows=[
            {
                "seqid": "ERZ1290913_80",
                "source": "Pyrodigal",
                "type": "CDS",
                "start": "1",
                "end": "90",
                "strand": "+",
                "id": "ERZ1290913_80_1",
                "product": "foo enzyme",
            },
            {
                "seqid": "ERZ1290913_80",
                "source": "FragGeneScanRS",
                "type": "CDS",
                "start": "100",
                "end": "189",
                "strand": "+",
                "id": "ERZ1290913_80_2_100_189_+",
                "product": "bar enzyme",
            },
        ],
    )

    # Protein FASTA provides translations by ID
    _write_faa_with_headers(
        prots,
        {
            "ERZ1290913_80_1": "MKT",
            "ERZ1290913_80_2_100_189_+": "VVV",
        },
    )

    records = build_records_from_gff(
        contigs_path=str(contigs),
        gff_path=str(gff),
        proteins_path=str(prots),
        prefix="genbank",
        default_product="hypothetical protein",
        locus_tag_prefix="LT_",
        gene_from_locus_tag=True,
        include_sources=None,
        exclude_sources=None,
        transl_table=11,
        drop_terminal_stop=True,
    )
    write_genbank(records, str(out_gbk))

    feats = _load_gbk_features(out_gbk)
    assert any(f.type == "source" for f in feats), "Expected contig-wide source feature"

    cds_map = _cds_by_protein_id(out_gbk)
    assert "ERZ1290913_80_1" in cds_map
    assert "ERZ1290913_80_2_100_189_+" in cds_map

    q1 = cds_map["ERZ1290913_80_1"]
    assert q1["protein_id"][0] == "ERZ1290913_80_1"
    assert q1["locus_tag"][0] == "LT_ERZ1290913_80_1"
    assert q1["product"][0] == "foo enzyme"
    assert q1["gene"][0] == "LT_ERZ1290913_80_1"
    assert q1["translation"][0] == "MKT"
    assert any("gene_caller=Pyrodigal" in n for n in q1.get("note", []))

    q2 = cds_map["ERZ1290913_80_2_100_189_+"]
    assert q2["translation"][0] == "VVV"
    assert any("gene_caller=FragGeneScanRS" in n for n in q2.get("note", []))


def test_gff_mode_include_sources_filters_by_column2(tmp_path: Path) -> None:
    contigs = tmp_path / "contigs.fna"
    gff = tmp_path / "calls.gff"
    prots = tmp_path / "proteins.faa"
    out_gbk = tmp_path / "out.gbk"

    _write_fasta(contigs, {"ERZ1290913_80": "A" * 2000})

    _write_gff3(
        gff,
        rows=[
            {
                "seqid": "ERZ1290913_80",
                "source": "Pyrodigal",
                "type": "CDS",
                "start": "1",
                "end": "90",
                "strand": "+",
                "id": "ERZ1290913_80_1",
            },
            {
                "seqid": "ERZ1290913_80",
                "source": "FragGeneScanRS",
                "type": "CDS",
                "start": "100",
                "end": "189",
                "strand": "+",
                "id": "ERZ1290913_80_2_100_189_+",
            },
        ],
    )

    _write_faa_with_headers(
        prots,
        {
            "ERZ1290913_80_1": "MKT",
            "ERZ1290913_80_2_100_189_+": "VVV",
        },
    )

    records = build_records_from_gff(
        contigs_path=str(contigs),
        gff_path=str(gff),
        proteins_path=str(prots),
        prefix="genbank",
        default_product="hypothetical protein",
        locus_tag_prefix="",
        gene_from_locus_tag=False,
        include_sources={"Pyrodigal"},
        exclude_sources=None,
        transl_table=11,
        drop_terminal_stop=True,
    )
    write_genbank(records, str(out_gbk))

    cds_map = _cds_by_protein_id(out_gbk)
    assert set(cds_map.keys()) == {"ERZ1290913_80_1"}


def test_gff_mode_fallback_translation_drops_terminal_stop(tmp_path: Path) -> None:
    """
    Ensure nucleotide fallback translation removes trailing '*' when the CDS ends with a stop codon.
    """
    contigs = tmp_path / "contigs.fna"
    gff = tmp_path / "calls.gff"
    out_gbk = tmp_path / "out.gbk"

    # Build a contig containing a CDS: ATG (M) + AAA (K) + TAA (stop)
    # Put it at positions 1..9, then pad.
    contig_seq = "ATGAAATAA" + ("A" * 200)
    _write_fasta(contigs, {"ERZ_STOP": contig_seq})

    _write_gff3(
        gff,
        rows=[
            {
                "seqid": "ERZ_STOP",
                "source": "Pyrodigal",
                "type": "CDS",
                "start": "1",
                "end": "9",
                "strand": "+",
                "id": "ERZ_STOP_1",
            }
        ],
    )

    # No proteins file provided => must fallback translate from contigs
    records = build_records_from_gff(
        contigs_path=str(contigs),
        gff_path=str(gff),
        proteins_path=None,
        prefix="genbank",
        default_product="hypothetical protein",
        locus_tag_prefix="",
        gene_from_locus_tag=False,
        include_sources=None,
        exclude_sources=None,
        transl_table=11,
        drop_terminal_stop=True,  # required behavior
    )
    write_genbank(records, str(out_gbk))

    cds_map = _cds_by_protein_id(out_gbk)
    tr = cds_map["ERZ_STOP_1"]["translation"][0]
    assert tr == "MK", f"Expected stop codon removed; got {tr!r}"
