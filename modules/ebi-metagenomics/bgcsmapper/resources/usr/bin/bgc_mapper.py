#!/usr/bin/env python3
"""
bgc_integrator.py  (table-based GFF parsing; NO BCBio)

Integrates BGC predictions from (optional) GECCO / antiSMASH / SanntiS outputs
into a single GFF3, using the *mandatory* base GFF.

Core rules:
  - Base CDS lines are preserved byte-for-byte for columns 1–9, except we append/replace ONLY:
        bgc_score
        bgc_tools
        + selected tool metadata keys
    in column 9.
  - Only CDS fully covered by at least one (merged) BGC region are written (i.e. only CDS with bgc_score).
  - Merge overlapping BGC regions per contig (data-portal rule).
  - BGC feature rows:
      * col3 type = bgc_region
      * ID first: ID=contig_id|bgc:start-end
      * if merged (members>1): col2 source = bgc_merged
      * if not merged (members==1): col2 source = original tool
      * ALWAYS add bgc_tools=<comma-separated tools> to bgc_region rows (merged and non-merged)
  - CDS rows:
      * ALWAYS add bgc_tools=<comma-separated tools that cover this CDS> when bgc_score is present.

Tool metadata propagation rules:
  - antiSMASH:
      * CDS entries: add ONLY the corresponding gene-level annotation for that CDS
        (matched by CDS ID == antiSMASH gene ID) plus the parent region product:
          antismash_gene_function      ← from antiSMASH gene_functions
          antismash_as_type            ← from as_type
          antismash_as_gene_clusters   ← from as_gene_clusters
          antismash_product            ← from parent region product
      * BGC entries:
          - antiSMASH contributes ONLY antismash_product to bgc_region rows
          - merged regions: antiSMASH contributes only product values (no gene-level aggregation)
  - GECCO:
      gecco_bgc_type               ← from GECCO Type (BGC + CDS)
  - SanntiS:
      nearest_MiBIG                ← from SanntiS nearest_MiBIG (BGC + CDS)
      nearest_MiBIG_class          ← from SanntiS nearest_MiBIG_class (BGC + CDS)
"""

from __future__ import annotations

import argparse
import fileinput
import logging
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

log = logging.getLogger(__name__)

# Keys that may be added to CDS/BGC attributes (order used when appending to CDS attributes)
META_KEYS_ORDER = [
    "antismash_gene_function",
    "antismash_product",
    "antismash_as_type",
    "antismash_as_gene_clusters",
    "gecco_bgc_type",
    "nearest_MiBIG",
    "nearest_MiBIG_class",
]


# ──────────────────────────────────────────────────────────────────────────────
# Data structures
# ──────────────────────────────────────────────────────────────────────────────


@dataclass(slots=True)
class CDSRec:
    contig: str
    start: int
    end: int
    line: str  # original base GFF line (no trailing newline)


@dataclass(slots=True)
class BGCRegion:
    contig: str
    start: int
    end: int
    tool: str  # gecco|sanntis|antismash
    attrs: dict[str, str] = field(
        default_factory=dict
    )  # region metadata (antiSMASH region keeps ONLY product)


@dataclass(slots=True)
class MergedRegion:
    contig: str
    start: int
    end: int
    members: list[BGCRegion]


# ──────────────────────────────────────────────────────────────────────────────
# CLI + validation
# ──────────────────────────────────────────────────────────────────────────────


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="bgc_integrator.py",
        description="Integrate optional GECCO/antiSMASH/SanntiS BGC calls into a base GFF3 (table parsing; no BCBio).",
    )
    p.add_argument(
        "--base_gff",
        required=True,
        type=Path,
        help="Mandatory base GFF (may be compressed).",
    )
    p.add_argument(
        "--gecco_gff",
        type=Path,
        default=None,
        help="Optional GECCO output GFF (may be compressed).",
    )
    p.add_argument(
        "--antismash_gff",
        type=Path,
        default=None,
        help="Optional antiSMASH output GFF (may be compressed).",
    )
    p.add_argument(
        "--sanntis_gff",
        type=Path,
        default=None,
        help="Optional SanntiS output GFF (may be compressed).",
    )
    p.add_argument(
        "--output_gff", required=True, type=Path, help="Output integrated GFF3."
    )
    p.add_argument(
        "--log_level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )
    return p.parse_args(argv)


def validate_inputs(args: argparse.Namespace) -> list[tuple[str, Path]]:
    if not args.base_gff.exists():
        raise FileNotFoundError(f"Base GFF not found: {args.base_gff}")

    optional: list[tuple[str, Path]] = []
    if args.gecco_gff:
        optional.append(("gecco", args.gecco_gff))
    if args.antismash_gff:
        optional.append(("antismash", args.antismash_gff))
    if args.sanntis_gff:
        optional.append(("sanntis", args.sanntis_gff))

    if not optional:
        raise ValueError(
            "At least one optional predictor GFF must be provided (gecco/antismash/sanntis)."
        )

    for tool, path in optional:
        if not path.exists():
            raise FileNotFoundError(f"{tool} GFF not found: {path}")

    return optional


# ──────────────────────────────────────────────────────────────────────────────
# GFF parsing utilities (table-based)
# ──────────────────────────────────────────────────────────────────────────────


def iter_gff_rows(path: Path) -> Iterable[tuple[str, list[str]]]:
    """Yield (raw_line_no_newline, cols) for non-comment GFF rows with 9 columns."""
    with fileinput.input(
        files=[str(path)], openhook=fileinput.hook_compressed
    ) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) != 9:
                log.debug("Skipping non-9col line in %s: %s", path, line)
                continue
            yield line, cols


def parse_attr_str(attr: str) -> dict[str, str]:
    if attr == "." or attr.strip() == "":
        return {}
    out: dict[str, str] = {}
    for part in attr.split(";"):
        if not part or "=" not in part:
            continue
        k, v = part.split("=", 1)
        out[k] = v
    return out


def _sanitize_attr_value(v: str) -> str:
    # Avoid raw ';' breaking attribute fields.
    return str(v).replace(";", ",")


def attrs_to_str_with_id_first(attrs: dict[str, str]) -> str:
    """Render attributes ensuring ID is first (used only for synthetic BGC region rows)."""
    if not attrs:
        return "."
    parts: list[str] = []
    if "ID" in attrs:
        parts.append(f"ID={_sanitize_attr_value(attrs['ID'])}")
    for k in sorted(attrs.keys()):
        if k == "ID":
            continue
        parts.append(f"{k}={_sanitize_attr_value(attrs[k])}")
    return ";".join(parts) if parts else "."


def _join_unique(values: Iterable[str]) -> str:
    uniq = sorted(
        {
            _sanitize_attr_value(v)
            for v in values
            if v is not None and str(v).strip() != ""
        }
    )
    return ",".join(uniq)


def extract_id_from_attr(attr: str) -> str | None:
    m = re.search(r"(?:(?<=;)|^)ID=([^;]+)", attr)
    return m.group(1) if m else None


def replace_or_append_attr(attr: str, key: str, value: str) -> str:
    """
    Replace key=... if present, else append ;key=value.
    Does NOT reformat or reorder anything else in the attribute string.
    """
    value = _sanitize_attr_value(value)
    if attr == "." or attr.strip() == "":
        return f"{key}={value}"

    if f"{key}=" in attr:
        pattern = rf"(?:(?<=;)|^){re.escape(key)}=[^;]*"
        return re.sub(pattern, f"{key}={value}", attr)

    return attr + f";{key}={value}"


# ──────────────────────────────────────────────────────────────────────────────
# Load base CDS lines (preserve exact original line)
# ──────────────────────────────────────────────────────────────────────────────


def load_base_cds(base_gff: Path) -> dict[str, list[CDSRec]]:
    contig_to_cds: dict[str, list[CDSRec]] = {}

    with fileinput.input(
        files=[str(base_gff)], openhook=fileinput.hook_compressed
    ) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) != 9:
                continue
            if cols[2] != "CDS":
                continue
            contig = cols[0]
            start = int(cols[3])
            end = int(cols[4])
            contig_to_cds.setdefault(contig, []).append(
                CDSRec(contig=contig, start=start, end=end, line=line)
            )

    for contig in contig_to_cds:
        contig_to_cds[contig].sort(key=lambda x: x.start)

    log.info(
        "Loaded base CDS: contigs=%d cds=%d",
        len(contig_to_cds),
        sum(len(v) for v in contig_to_cds.values()),
    )
    return contig_to_cds


# ──────────────────────────────────────────────────────────────────────────────
# Parse predictor outputs into BGC regions (table-based)
# ──────────────────────────────────────────────────────────────────────────────


def parse_gecco_regions(path: Path) -> list[BGCRegion]:
    regs: list[BGCRegion] = []
    for _, cols in iter_gff_rows(path):
        contig, _source, _ftype, start_s, end_s, *_rest, attr_s = cols
        attrs = parse_attr_str(attr_s)
        if "Type" in attrs:
            regs.append(
                BGCRegion(
                    contig=contig,
                    start=int(start_s),
                    end=int(end_s),
                    tool="gecco",
                    attrs={"gecco_bgc_type": attrs["Type"]},
                )
            )
    log.info("Parsed GECCO regions: %d", len(regs))
    return regs


def parse_sanntis_regions(path: Path) -> list[BGCRegion]:
    regs: list[BGCRegion] = []
    for _, cols in iter_gff_rows(path):
        contig, _source, _ftype, start_s, end_s, *_rest, attr_s = cols
        attrs = parse_attr_str(attr_s)
        if "nearest_MiBIG" in attrs or "nearest_MiBIG_class" in attrs:
            ra: dict[str, str] = {}
            if "nearest_MiBIG" in attrs:
                ra["nearest_MiBIG"] = attrs["nearest_MiBIG"]
            if "nearest_MiBIG_class" in attrs:
                ra["nearest_MiBIG_class"] = attrs["nearest_MiBIG_class"]
            regs.append(
                BGCRegion(
                    contig=contig,
                    start=int(start_s),
                    end=int(end_s),
                    tool="sanntis",
                    attrs=ra,
                )
            )
    log.info("Parsed SanntiS regions: %d", len(regs))
    return regs


def parse_antismash_regions_and_genes(
    path: Path,
) -> tuple[list[BGCRegion], dict[str, dict[str, str]]]:
    """
    Parse antiSMASH and return:
      1) regions: BGCRegion objects (one per antiSMASH region interval)
         - IMPORTANT: region attrs include ONLY antismash_product (and optional antismash_region_id)
      2) gene_ann_by_id: gene-level annotations keyed by gene ID (for CDS entries only)
    """
    regions_by_id: dict[str, tuple[str, int, int, str | None]] = {}
    gene_parent_by_id: dict[str, str] = {}
    gene_ann_by_id: dict[str, dict[str, str]] = {}

    for _, cols in iter_gff_rows(path):
        contig, _source, ftype, start_s, end_s, *_rest, attr_s = cols
        attrs = parse_attr_str(attr_s)

        if ftype in ("region", "biosynthetic-gene-cluster"):
            rid = attrs.get("ID")
            if not rid:
                continue
            product = attrs.get("product")
            regions_by_id[rid] = (contig, int(start_s), int(end_s), product)
            continue

        if ftype == "gene":
            gid = attrs.get("ID")
            parent = attrs.get("Parent")
            if not gid or not parent:
                continue

            gene_parent_by_id[gid] = parent
            g: dict[str, str] = {}

            if "gene_functions" in attrs:
                g["antismash_gene_function"] = attrs["gene_functions"]
            if "as_type" in attrs:
                g["antismash_as_type"] = attrs["as_type"]
            if "as_gene_clusters" in attrs:
                g["antismash_as_gene_clusters"] = attrs["as_gene_clusters"]

            if g:
                gene_ann_by_id[gid] = g

    regions: list[BGCRegion] = []
    for rid, (contig, start, end, product) in regions_by_id.items():
        ra: dict[str, str] = {"antismash_region_id": rid}
        if product:
            ra["antismash_product"] = product
        regions.append(
            BGCRegion(contig=contig, start=start, end=end, tool="antismash", attrs=ra)
        )

    for gid, parent in gene_parent_by_id.items():
        if parent in regions_by_id:
            _c, _s, _e, product = regions_by_id[parent]
            if product:
                gene_ann_by_id.setdefault(gid, {})
                gene_ann_by_id[gid]["antismash_product"] = product

    log.info("Parsed antiSMASH regions: %d", len(regions))
    return regions, gene_ann_by_id


# ──────────────────────────────────────────────────────────────────────────────
# Merge overlaps + scoring + metadata propagation to CDS
# ──────────────────────────────────────────────────────────────────────────────


def merge_overlaps(regions: Sequence[BGCRegion]) -> list[MergedRegion]:
    """Merge any overlapping intervals per contig (data-portal rule)."""
    by_contig: dict[str, list[BGCRegion]] = {}
    for r in regions:
        by_contig.setdefault(r.contig, []).append(r)

    merged: list[MergedRegion] = []
    for contig, regs in by_contig.items():
        if not regs:
            continue
        regs_sorted = sorted(regs, key=lambda r: r.start)

        cur_start, cur_end = regs_sorted[0].start, regs_sorted[0].end
        cur_members = [regs_sorted[0]]

        for r in regs_sorted[1:]:
            if r.start <= cur_end:
                cur_end = max(cur_end, r.end)
                cur_members.append(r)
            else:
                merged.append(
                    MergedRegion(
                        contig=contig, start=cur_start, end=cur_end, members=cur_members
                    )
                )
                cur_start, cur_end, cur_members = r.start, r.end, [r]

        merged.append(
            MergedRegion(
                contig=contig, start=cur_start, end=cur_end, members=cur_members
            )
        )

    return merged


def cds_within_region(start: int, end: int, rstart: int, rend: int) -> bool:
    return start >= rstart and end <= rend


def _tools_covering_cds(
    members: list[BGCRegion], cds_start: int, cds_end: int
) -> list[str]:
    """Return sorted unique tool names for member predictions that fully cover this CDS."""
    tools = sorted(
        {
            m.tool
            for m in members
            if cds_within_region(cds_start, cds_end, m.start, m.end)
        }
    )
    return tools


def _collect_member_meta_for_cds(
    members: list[BGCRegion], cds_start: int, cds_end: int
) -> dict[str, str]:
    """
    Collect metadata keys for a CDS from member predictions that cover the CDS.
    Values are unioned (unique values joined by commas).

    NOTE: antiSMASH BGC members only carry antismash_product at region-level (by design).
          antiSMASH gene-level keys are applied separately via CDS ID matching.
    """
    collected: dict[str, list[str]] = {k: [] for k in META_KEYS_ORDER}
    for m in members:
        if not cds_within_region(cds_start, cds_end, m.start, m.end):
            continue
        for k in META_KEYS_ORDER:
            if k in m.attrs and m.attrs[k]:
                collected[k].append(m.attrs[k])

    out: dict[str, str] = {}
    for k, vals in collected.items():
        if vals:
            out[k] = _join_unique(vals)
    return out


def score_and_filter_cds(
    contig_to_cds: dict[str, list[CDSRec]],
    merged_regions: Sequence[MergedRegion],
    antismash_gene_ann_by_id: dict[str, dict[str, str]],
) -> dict[str, list[str]]:
    """
    For each CDS fully inside a merged region:
      - compute bgc_score = (# member BGCs covering CDS) / (# member BGCs in region)
      - add bgc_tools = comma-separated list of tools covering the CDS
      - add tool metadata keys (from covering member regions)
      - apply antiSMASH gene-level keys ONLY for the matching CDS ID
      - output only those CDS (filtering out non-BGC CDS)
    """
    out: dict[str, list[str]] = {}

    for mr in merged_regions:
        cds_list = contig_to_cds.get(mr.contig, [])
        if not cds_list or not mr.members:
            continue

        denom = len(mr.members)

        for cds in cds_list:
            if cds.start > mr.end:
                break
            if not cds_within_region(cds.start, cds.end, mr.start, mr.end):
                continue

            hits = 0
            for m in mr.members:
                if cds_within_region(cds.start, cds.end, m.start, m.end):
                    hits += 1
            score = hits / denom

            cols = cds.line.split("\t")
            attr = cols[8]

            # 1) bgc_tools: tools that cover this CDS
            tools = _tools_covering_cds(mr.members, cds.start, cds.end)
            if tools:
                attr = replace_or_append_attr(attr, "bgc_tools", ",".join(tools))

            # 2) Metadata from region-members that cover the CDS
            meta = _collect_member_meta_for_cds(mr.members, cds.start, cds.end)

            # 3) antiSMASH gene-level metadata ONLY for this CDS (ID match)
            cds_id = extract_id_from_attr(attr)
            if cds_id and cds_id in antismash_gene_ann_by_id:
                gene_meta = antismash_gene_ann_by_id[cds_id]
                for k in (
                    "antismash_gene_function",
                    "antismash_as_type",
                    "antismash_as_gene_clusters",
                    "antismash_product",
                ):
                    if k in gene_meta and gene_meta[k]:
                        meta[k] = gene_meta[k]

            for k in META_KEYS_ORDER:
                if k in meta:
                    attr = replace_or_append_attr(attr, k, meta[k])

            # bgc_score last
            attr = replace_or_append_attr(attr, "bgc_score", f"{score:.2f}")

            cols[8] = attr
            out.setdefault(mr.contig, []).append("\t".join(cols))

    for contig in out:
        out[contig].sort(key=lambda ln: int(ln.split("\t")[3]))
    return out


# ──────────────────────────────────────────────────────────────────────────────
# Output region lines (BGC rows)
# ──────────────────────────────────────────────────────────────────────────────


def _collect_member_meta_for_region(members: list[BGCRegion]) -> dict[str, str]:
    """
    Union metadata keys across member regions (unique values joined by commas).
    antiSMASH contributes only antismash_product at region-level.
    """
    collected: dict[str, list[str]] = {k: [] for k in META_KEYS_ORDER}
    for m in members:
        for k in META_KEYS_ORDER:
            if k in m.attrs and m.attrs[k]:
                collected[k].append(m.attrs[k])

    out: dict[str, str] = {}
    for k, vals in collected.items():
        if vals:
            out[k] = _join_unique(vals)
    return out


def build_region_lines(
    merged_regions: Sequence[MergedRegion],
) -> list[tuple[str, int, str]]:
    """
    Emit region features:
      - members > 1 => source 'bgc_merged'
      - members == 1 => source original tool
      - type 'bgc_region'
      - ID first: contig|bgc:start-end
      - ALWAYS add bgc_tools=<comma-separated tools in this region> (merged and non-merged)
    """
    lines: list[tuple[str, int, str]] = []

    for mr in merged_regions:
        if not mr.members:
            continue

        tools = sorted({m.tool for m in mr.members})

        attrs: dict[str, str] = {
            "ID": f"{mr.contig}|bgc:{mr.start}-{mr.end}",
            "bgc_tools": ",".join(tools),
        }

        if len(mr.members) > 1:
            source = "bgc_merged"
            attrs["member_bgcs"] = str(len(mr.members))
            attrs.update(_collect_member_meta_for_region(mr.members))
        else:
            source = mr.members[0].tool
            attrs.update(mr.members[0].attrs)

        line = "\t".join(
            [
                mr.contig,
                source,
                "bgc_region",
                str(mr.start),
                str(mr.end),
                ".",
                ".",
                ".",
                attrs_to_str_with_id_first(attrs),
            ]
        )
        lines.append((mr.contig, mr.start, line))

    return lines


def write_gff(
    output: Path,
    region_lines: list[tuple[str, int, str]],
    cds_lines_by_contig: dict[str, list[str]],
) -> None:
    sortable: list[tuple[str, int, str]] = []
    sortable.extend(region_lines)

    for contig, cds_lines in cds_lines_by_contig.items():
        for ln in cds_lines:
            start = int(ln.split("\t")[3])
            sortable.append((contig, start, ln))

    sortable.sort(key=lambda x: (x[0], x[1]))

    with output.open("w") as out:
        out.write("##gff-version 3\n")
        for _, _, ln in sortable:
            out.write(ln + "\n")


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    try:
        optional_inputs = validate_inputs(args)
    except Exception as e:
        log.error(str(e))
        return 2

    contig_to_cds = load_base_cds(args.base_gff)

    regions: list[BGCRegion] = []
    antismash_gene_ann_by_id: dict[str, dict[str, str]] = {}

    for tool, path in optional_inputs:
        if tool == "gecco":
            regions.extend(parse_gecco_regions(path))
        elif tool == "sanntis":
            regions.extend(parse_sanntis_regions(path))
        elif tool == "antismash":
            a_regs, a_gene_ann = parse_antismash_regions_and_genes(path)
            regions.extend(a_regs)
            antismash_gene_ann_by_id.update(a_gene_ann)
        else:
            log.warning("Unknown tool %s (skipping)", tool)

    if not regions:
        log.warning("No BGC regions parsed from optional inputs; writing header only.")
        write_gff(args.output_gff, [], {})
        return 0

    merged_regions = merge_overlaps(regions)

    cds_lines_by_contig = score_and_filter_cds(
        contig_to_cds=contig_to_cds,
        merged_regions=merged_regions,
        antismash_gene_ann_by_id=antismash_gene_ann_by_id,
    )

    region_lines = build_region_lines(merged_regions)
    write_gff(args.output_gff, region_lines, cds_lines_by_contig)

    log.info("Wrote integrated GFF: %s", args.output_gff)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
