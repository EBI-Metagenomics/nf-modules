#!/usr/bin/env python3
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


"""
bgc_mapper.py

Integrates BGC predictions from (optional) GECCO / antiSMASH / SanntiS outputs
into a single GFF3, using the mandatory base GFF.

Core rules:
  - Base CDS lines are preserved byte-for-byte for columns 1–9, except we append/replace ONLY:
        bgc_support
        bgc_tools
        + selected tool metadata keys
    in column 9.
  - Only CDS fully covered by at least one (merged) BGC region are written (i.e. only CDS with bgc_support).
  - Merge overlapping BGC regions per contig (this is based on the MGnify BGCs data-portal rule).
  - BGC feature rows:
      * col3 type = bgc_region
      * ID first attribute: ID=contig_id|bgc:start-end
      * if merged (members>1): col2 source = bgc_merged
      * if not merged (members==1): col2 source = EXACT original source string from predictor GFF (may include version)
      * Always add bgc_tools=<comma-separated tools> to bgc_region rows (merged and non-merged)
  - CDS rows:
      * Always add bgc_tools=<comma-separated tools that cover this CDS>.

Tool metadata propagation rules:
  - antiSMASH:
      * CDS entries: add the corresponding gene-level annotation for that CDS
        (matched by CDS ID == antiSMASH gene ID) plus the parent region product:
          antismash_gene_function      ← from antiSMASH gene_functions
          antismash_as_type            ← from as_type
          antismash_as_gene_clusters   ← from as_gene_clusters
          antismash_product            ← from parent region product
      * BGC entries:
          - antiSMASH contributes antismash_product to bgc_region rows
          - merged regions: antiSMASH contributes only product values
  - GECCO:
      gecco_bgc_type ← from GECCO Type (BGC + CDS)
  - SanntiS:
      nearest_MiBIG ← from SanntiS nearest_MiBIG (BGC + CDS)
      nearest_MiBIG_class ← from SanntiS nearest_MiBIG_class (BGC + CDS)

Support/scoring rule (FIXED):
  - Let N = number of tools provided as input (min=1, max=3).
  - For each CDS covered by at leats one bgc prediction, compute:
        bgc_support = (# DISTINCT tools whose predictions cover the CDS) / N
    This yields the following valid scores:
      - N=3 → 0.33, 0.67, 1.00
      - N=2 → 0.50, 1.00
      - N=1 → 1.00

Optional json output:
  - Write antiSMASH sideloader JSON (records/subregions/details only; no protoclusters),
    with same basename as the output GFF.
  - Validate JSON against the official antiSMASH schemas (vendored locally).
"""

from __future__ import annotations

import argparse
import fileinput
import json
import logging
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from datetime import datetime, timezone
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
    tool: str  # gecco|sanntis|antismash (normalized tool key)
    source: str  # EXACT column 2 source string from the original predictor GFF (may include version)
    attrs: dict[str, str] = field(default_factory=dict)


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
        prog="bgc_mapper.py",
        description="Integrate optional GECCO/antiSMASH/SanntiS BGC calls into a base GFF3 file",
    )
    p.add_argument(
        "--base_gff",
        required=True,
        type=Path,
        help="Mandatory base GFF containing the coordinates of the original CDS used for BGC prediction (may be compressed).",
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
        "--validate_json",
        action="store_true",
        default=False,
        help="Validate the produced antiSMASH sideloader JSON against the local schema copies (dev/debug). Disabled by default.",
    )
    p.add_argument(
        "--output_gff", required=True, type=Path, help="Output integrated GFF3."
    )
    p.add_argument(
        "--verbose", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )
    return p.parse_args(argv)


def validate_inputs(args: argparse.Namespace) -> list[tuple[str, Path]]:
    print("validating inputs")

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
# GFF parsing utilities
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
    print("parsing base gff file")

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
# Parse predictor outputs into BGC regions
# ──────────────────────────────────────────────────────────────────────────────


def parse_gecco_regions(path: Path) -> list[BGCRegion]:
    print("Parsing gecco annotations")

    regs: list[BGCRegion] = []
    for _, cols in iter_gff_rows(path):
        contig, source, _ftype, start_s, end_s, *_rest, attr_s = cols
        attrs = parse_attr_str(attr_s)
        if "Type" in attrs:
            regs.append(
                BGCRegion(
                    contig=contig,
                    start=int(start_s),
                    end=int(end_s),
                    tool="gecco",
                    source=source,  # preserve original column 2
                    attrs={"gecco_bgc_type": attrs["Type"]},
                )
            )
    log.info("Parsed GECCO regions: %d", len(regs))
    return regs


def parse_sanntis_regions(path: Path) -> list[BGCRegion]:
    print("Parsing sanntis annotations")

    regs: list[BGCRegion] = []
    for _, cols in iter_gff_rows(path):
        contig, source, _ftype, start_s, end_s, *_rest, attr_s = cols
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
                    source=source,  # preserve original column 2
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
         - IMPORTANT: region.source preserves original antiSMASH column 2 (often includes version)
      2) gene_ann_by_id: gene-level annotations keyed by gene ID (for CDS entries only)
    """

    print("Parsing antismash annotations")

    regions_by_id: dict[str, tuple[str, int, int, str | None, str]] = {}
    gene_parent_by_id: dict[str, str] = {}
    gene_ann_by_id: dict[str, dict[str, str]] = {}

    for _, cols in iter_gff_rows(path):
        contig, source, ftype, start_s, end_s, *_rest, attr_s = cols
        attrs = parse_attr_str(attr_s)

        if ftype in ("region", "biosynthetic-gene-cluster"):
            rid = attrs.get("ID")
            if not rid:
                continue
            product = attrs.get("product")
            # store source too
            regions_by_id[rid] = (contig, int(start_s), int(end_s), product, source)
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
    for rid, (contig, start, end, product, source) in regions_by_id.items():
        ra: dict[str, str] = {"antismash_region_id": rid}
        if product:
            ra["antismash_product"] = product
        regions.append(
            BGCRegion(
                contig=contig,
                start=start,
                end=end,
                tool="antismash",
                source=source,  # preserve original column 2
                attrs=ra,
            )
        )

    for gid, parent in gene_parent_by_id.items():
        if parent in regions_by_id:
            _c, _s, _e, product, _src = regions_by_id[parent]
            if product:
                gene_ann_by_id.setdefault(gid, {})
                gene_ann_by_id[gid]["antismash_product"] = product

    log.info("Parsed antiSMASH regions: %d", len(regions))
    return regions, gene_ann_by_id


# ──────────────────────────────────────────────────────────────────────────────
# Merge overlaps + support score + metadata propagation to CDS
# ──────────────────────────────────────────────────────────────────────────────


def merge_overlaps(regions: Sequence[BGCRegion]) -> list[MergedRegion]:
    """Merge any overlapping intervals per contig (data-portal rule)."""

    print("Merging overlapping regions")

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
    return sorted(
        {
            m.tool
            for m in members
            if cds_within_region(cds_start, cds_end, m.start, m.end)
        }
    )


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


def support_and_filter_cds(
    contig_to_cds: dict[str, list[CDSRec]],
    merged_regions: Sequence[MergedRegion],
    antismash_gene_ann_by_id: dict[str, dict[str, str]],
    n_tools: int,
) -> dict[str, list[str]]:
    """
    For each CDS fully inside a merged region:
      - compute bgc_support = (# DISTINCT tools covering CDS) / (n_tools provided as input)
      - add bgc_tools = comma-separated list of tools covering the CDS
      - add tool metadata keys (from covering member regions)
      - apply antiSMASH gene-level keys ONLY for the matching CDS ID
      - output only those CDS (filtering out non-BGC CDS)
    """
    if n_tools < 1:
        raise ValueError("n_tools must be >= 1")

    out: dict[str, list[str]] = {}
    for mr in merged_regions:
        cds_list = contig_to_cds.get(mr.contig, [])
        if not cds_list or not mr.members:
            continue

        for cds in cds_list:
            if cds.start > mr.end:
                break
            if not cds_within_region(cds.start, cds.end, mr.start, mr.end):
                continue

            tools = _tools_covering_cds(mr.members, cds.start, cds.end)
            if not tools:
                continue

            bgc_support = len(tools) / n_tools

            cols = cds.line.split("\t")
            attr = cols[8]

            # bgc_tools
            attr = replace_or_append_attr(attr, "bgc_tools", ",".join(tools))

            # metadata from members that cover CDS
            meta = _collect_member_meta_for_cds(mr.members, cds.start, cds.end)

            # antiSMASH gene-level metadata ONLY for this CDS (ID match)
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

            # bgc_support last
            attr = replace_or_append_attr(attr, "bgc_support", f"{bgc_support:.2f}")

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
      - members == 1 => source EXACT original source string from predictor GFF (member.source)
      - type 'bgc_region'
      - ID first: contig|bgc:start-end
      - Add bgc_tools= (merged and non-merged)
    """

    print("building region lines")

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
            # preserve the *exact* source string from the original predictor file (col2)
            source = mr.members[0].source
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
    print("Writing the output in gff format")

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
# antiSMASH sideloader JSON output (subregions only) + jsonschema validation
# ──────────────────────────────────────────────────────────────────────────────


def _derive_output_json_path(output_gff: Path) -> Path:
    """Derive a JSON output path with the same basename as *output_gff*.

    Examples:
      - out.gff   -> out.json
      - out.gff3  -> out.json
      - out.gff.gz -> out.json
    """
    name = output_gff.name
    # Handle compressed gff names (even if we don't write compressed output ourselves).
    if name.endswith(".gff.gz") or name.endswith(".gff3.gz"):
        base = name.rsplit(".", 2)[0]  # drop .gff(.3).gz
        return output_gff.with_name(f"{base}.json")
    if name.endswith(".gff") or name.endswith(".gff3"):
        base = name.rsplit(".", 1)[0]
        return output_gff.with_name(f"{base}.json")
    return output_gff.with_suffix(".json")


def _default_antismash_schema_dir() -> Path:
    """Expected location for vendored antiSMASH sideloader schemas."""
    return Path(__file__).absolute().parent / "antismash_sideload_schemas" / "general"


def _fix_refs(node):
    """Recursively strip 'file:' from $ref values (antiSMASH uses file: URIs)."""
    if isinstance(node, dict):
        out = {}
        for k, v in node.items():
            if k == "$ref" and isinstance(v, str):
                out[k] = v.replace("file:", "")
            else:
                out[k] = _fix_refs(v)
        return out
    if isinstance(node, list):
        return [_fix_refs(x) for x in node]
    return node


def _load_schema_files(schema_dir: Path) -> tuple[dict, dict[str, dict]]:
    """Load root schema.json + subschemas into a resolver store."""
    schema_path = schema_dir / "schema.json"
    subs_dir = schema_dir / "subschemas"
    if not schema_path.exists():
        raise FileNotFoundError(f"schema.json not found at: {schema_path}")
    if not subs_dir.exists():
        raise FileNotFoundError(f"subschemas/ dir not found at: {subs_dir}")

    def _read(p: Path) -> dict:
        return _fix_refs(json.loads(p.read_text(encoding="utf-8")))

    root = _read(schema_path)

    store: dict[str, dict] = {}

    def _register_keys_for_doc(p: Path, doc: dict) -> None:
        # Absolute file URI for *actual* file
        store[p.resolve().as_uri()] = doc

        # Common relative keys
        store[p.name] = doc
        store[str(p)] = doc

        # If this is a subschema, also register it under subschemas/<name>
        if p.parent.name == "subschemas":
            store[f"subschemas/{p.name}"] = doc

            # ALSO: register as if it were located at the root schema_dir
            # (this covers refs like "tool.json" resolving to <schema_dir>/tool.json)
            fake_root_path = (schema_dir / p.name).resolve()
            store[fake_root_path.as_uri()] = doc
            store[str(fake_root_path)] = doc

    # root
    _register_keys_for_doc(schema_path, root)

    # subschemas
    for p in sorted(subs_dir.glob("*.json")):
        doc = _read(p)
        _register_keys_for_doc(p, doc)

    return root, store


def validate_sideload_json(payload: dict, schema_dir: Path | None = None) -> None:
    """Validate *payload* against the official antiSMASH sideloader schema (local copies).

    Uses jsonschema Draft7Validator with a referencing.Registry to resolve $ref
    from in-memory schemas
    """
    print("Validating json format")

    try:
        from jsonschema import Draft7Validator  # type: ignore
        from referencing import Registry, Resource  # type: ignore
    except Exception as e:
        raise RuntimeError(
            "Validation requires 'jsonschema' (and its 'referencing' dependency). "
            "Install it (e.g. pip install jsonschema) or disable validation."
        ) from e

    schema_dir = _default_antismash_schema_dir() if schema_dir is None else schema_dir
    root, store = _load_schema_files(schema_dir)

    # Build a registry that can resolve $ref values seen in antiSMASH schemas.
    # We register *every* key from store as a Resource.
    registry = Registry()
    for key, doc in store.items():
        registry = registry.with_resource(key, Resource.from_contents(doc))

    # Validate (raises jsonschema.ValidationError on failure)
    validator = Draft7Validator(schema=root, registry=registry)
    validator.validate(payload)


def _sanitize_details_value(value: str) -> str:
    """Make sure details values are schema-friendly (string, no leading whitespace)."""
    v = str(value).strip()
    # details.json disallows leading '=' '_' ',' or whitespace; we strip and then replace if needed
    while v and v[0] in ("_", "=", ","):
        v = v[1:]
    return v if v else "NA"


def _choose_subregion_label(mr: MergedRegion) -> str:
    """Choose a conservative <=20 char label for a subregion."""
    if len(mr.members) > 1:
        return "MergedBGC"
    m = mr.members[0]
    if m.tool == "antismash" and m.attrs.get("antismash_product"):
        lbl = str(m.attrs["antismash_product"])
    elif m.tool == "gecco" and m.attrs.get("gecco_bgc_type"):
        lbl = str(m.attrs["gecco_bgc_type"])
    elif m.tool == "sanntis" and m.attrs.get("nearest_MiBIG_class"):
        lbl = str(m.attrs["nearest_MiBIG_class"])
    else:
        lbl = m.tool
    # Keep it short and avoid whitespace-heavy labels
    lbl = re.sub(r"\s+", "_", lbl)
    lbl = re.sub(r"[^A-Za-z0-9_\-\.]", "_", lbl)
    return (lbl[:20] if len(lbl) > 20 else lbl) or "BGC"


def build_sideload_json_payload(
    merged_regions: Sequence[MergedRegion],
    *,
    tool_name: str,
    tool_version: str,
    tool_description: str,
) -> dict:
    """Build an antiSMASH sideloader JSON payload using only records/subregions/details."""
    records_by_name: dict[str, list[dict]] = {}

    for mr in merged_regions:
        if not mr.members:
            continue

        # Convert from GFF (1-based inclusive) to sideloader JSON (0-based start, end-exclusive)
        start0 = max(0, int(mr.start) - 1)
        end_excl = int(mr.end)

        tools = sorted({m.tool for m in mr.members})
        sources = sorted({m.source for m in mr.members if m.source})

        # Union member attrs as list-of-strings (schema allows str or list[str]; lists preserve multiplicity)
        union: dict[str, set[str]] = {}
        for m in mr.members:
            for k, v in m.attrs.items():
                if v is None or str(v).strip() == "":
                    continue
                union.setdefault(k, set()).add(_sanitize_details_value(v))

        details: dict[str, object] = {
            "ID": _sanitize_details_value(f"{mr.contig}|bgc:{mr.start}-{mr.end}"),
            "bgc_tools": [_sanitize_details_value(t) for t in tools],
        }
        if len(mr.members) > 1:
            details["member_bgcs"] = _sanitize_details_value(str(len(mr.members)))
        if sources:
            details["sources"] = [_sanitize_details_value(s) for s in sources]

        for k, vals in union.items():
            # Avoid clobbering our own keys
            if k in details:
                k = f"member_{k}"
            details[k] = sorted(vals)

        subregion = {
            "start": start0,
            "end": end_excl,
            "label": _choose_subregion_label(mr),
            "details": details,
        }
        records_by_name.setdefault(mr.contig, []).append(subregion)

    records: list[dict] = []
    for name in sorted(records_by_name.keys()):
        subregions = sorted(records_by_name[name], key=lambda d: int(d["start"]))
        records.append({"name": name, "subregions": subregions})

    payload = {
        "tool": {
            "name": tool_name,
            "version": tool_version if str(tool_version).strip() else "unknown",
            "description": tool_description,
            "configuration": {},  # keep empty by default; can be filled later if needed
        },
        "records": records,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    return payload


def write_sideload_json(
    out_json: Path,
    merged_regions: Sequence[MergedRegion],
    *,
    schema_dir: Path | None = None,
    validate: bool = False,
    tool_name: str = "MGnify BGC mapper",
    tool_version: str = "mgnify_pipelines_toolkit_v1.4.18",
    tool_description: str = "BGC subregions integrated from GECCO/antiSMASH/SanntiS into a base GFF.",
) -> None:
    """Write sideloader JSON. Optionally validate against the official antiSMASH schema files."""

    print("writing json output")

    payload = build_sideload_json_payload(
        merged_regions,
        tool_name=tool_name,
        tool_version=tool_version,
        tool_description=tool_description,
    )

    # Write first (so you still get an artifact even if validation is slow/fails when enabled)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(
        json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8"
    )
    log.info("Wrote sideloader JSON: %s", out_json)

    # Validate only if requested
    if validate:
        validate_sideload_json(payload, schema_dir=schema_dir)


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.verbose),
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    try:
        optional_inputs = validate_inputs(args)
    except Exception as e:
        log.error(str(e))
        return 2

    # number of provided tools (min=1, max=3) used as denominator for support
    n_tools = len(optional_inputs)

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

    cds_lines_by_contig = support_and_filter_cds(
        contig_to_cds=contig_to_cds,
        merged_regions=merged_regions,
        antismash_gene_ann_by_id=antismash_gene_ann_by_id,
        n_tools=n_tools,
    )

    region_lines = build_region_lines(merged_regions)
    write_gff(args.output_gff, region_lines, cds_lines_by_contig)

    log.info("Wrote integrated GFF: %s", args.output_gff)

    # Write antiSMASH-compatible sideloader JSON (subregions only), same basename as output_gff
    out_json = _derive_output_json_path(args.output_gff)
    try:
        write_sideload_json(out_json, merged_regions, validate=args.validate_json)
    except Exception as e:
        log.error("Failed to write/validate sideloader JSON: %s", e)
        return 3

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
