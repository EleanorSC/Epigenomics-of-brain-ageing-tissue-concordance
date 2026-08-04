#!/usr/bin/env python3
"""
04_extract_sommerer_buccal_brain.py

Extract published CpG-level buccal–brain DNA methylation correlations from:

Sommerer Y, Ohlei O, Dobricic V, et al.
A correlation map of genome-wide DNA methylation patterns between paired
human brain and buccal samples. Clinical Epigenetics. 2022;14:118.
DOI: 10.1186/s13148-022-01357-w

Input
-----
The script is designed for the official supplementary workbook:

    Sommerer_13148_2022_1357_MOESM1_ESM.xlsx

It uses:

- Supplementary Table S6:
  CpG-level Spearman correlations between paired prefrontal cortex and buccal
  samples for probes passing q < 0.05.

- Supplementary Table S2:
  chromosome and genomic coordinate, where available.

Important scope limitation
--------------------------
Supplementary Table S6 is already filtered by the original study to q < 0.05.
This script therefore extracts all published S6 observations, but does not
reconstruct the complete genome-wide correlation distribution.

No additional correlation threshold is applied by default.

Outputs
-------
A harmonised long-format CSV containing one row per CpG, suitable for merging
into the unified brain–peripheral concordance database.

Optional probe annotation can be supplied to add gene, genomic context and
CpG-island relation.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

import pandas as pd


SOURCE = "Sommerer et al. (2022), Clinical Epigenetics 14:118"
SOURCE_DOI = "10.1186/s13148-022-01357-w"
SOURCE_URL = "http://www.liga.uni-luebeck.de/buccal_brain_correlation_results/"
DATABASE_TOOL = "MADRC Buccal–Brain Correlation Map"
PERIPHERAL_TISSUE = "Buccal epithelium"
BRAIN_TISSUE = "Prefrontal cortex"
CORR_TYPE = "Spearman rho"
ANALYSIS_TYPE = "Published q-significant within-subject correlation"
SOURCE_SCOPE = (
    "CpGs reported in Sommerer et al. (2022) Supplementary Table S6 "
    "(q < 0.05); not the complete genome-wide correlation dataset"
)

OUTPUT_COLUMNS = [
    "cpg",
    "chr",
    "pos",
    "gene",
    "context",
    "island_relation",
    "peripheral_tissue",
    "brain_tissue",
    "tissue_pair",
    "corr_type",
    "corr_value",
    "abs_corr",
    "p_value",
    "q_value",
    "n_pairs",
    "analysis_type",
    "buccal_mqtl",
    "blood_mqtl",
    "brain_mqtl",
    "corsiv",
    "braun_et_al",
    "source",
    "source_doi",
    "source_url",
    "database_tool",
    "source_table",
    "source_scope",
]


def canonical(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def numeric(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def normalise_chr(value: Any) -> str:
    return re.sub(r"^chr", "", clean_text(value), flags=re.IGNORECASE)


def normalise_position(value: Any) -> str:
    number = numeric(value)
    if number is not None and number.is_integer():
        return str(int(number))
    return clean_text(value)


def yes_no(value: Any) -> str:
    text = clean_text(value).lower()
    if text in {"yes", "y", "true", "1"}:
        return "yes"
    if text in {"no", "n", "false", "0"}:
        return "no"
    return clean_text(value)


def find_header_row(raw: pd.DataFrame, required_terms: Sequence[str]) -> int:
    best_row = 0
    best_score = -1
    for row_index in range(min(30, len(raw))):
        cells = [canonical(value) for value in raw.iloc[row_index].tolist()]
        score = sum(any(term in cell for cell in cells) for term in required_terms)
        if score > best_score:
            best_score = score
            best_row = row_index
    return best_row


def read_excel_table(
    path: Path,
    sheet_name: str,
    required_terms: Sequence[str],
) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name=sheet_name, header=None)
    header_row = find_header_row(raw, required_terms)
    headers = [
        clean_text(value) or f"unnamed_{index}"
        for index, value in enumerate(raw.iloc[header_row].tolist())
    ]
    frame = raw.iloc[header_row + 1 :].copy()
    frame.columns = headers
    return frame.dropna(how="all")


def map_columns(frame: pd.DataFrame) -> dict[str, str]:
    return {canonical(column): column for column in frame.columns}


def first_match(
    columns: Mapping[str, str],
    aliases: Sequence[str],
    contains: Sequence[str] = (),
) -> str | None:
    for alias in aliases:
        if alias in columns:
            return columns[alias]
    for token in contains:
        for canonical_name, original_name in columns.items():
            if token in canonical_name:
                return original_name
    return None


def read_s2_coordinates(path: Path) -> dict[str, dict[str, str]]:
    """Read CpG chromosome and position from Supplementary Table S2."""
    frame = read_excel_table(
        path,
        "S2",
        required_terms=("cpg", "chr", "mapinfo", "rho"),
    )
    columns = map_columns(frame)

    cpg_col = first_match(columns, ["cpg_id", "cpg", "probe_id"], ["cpg"])
    chr_col = first_match(columns, ["chr", "chromosome"], ["chr", "chromosome"])
    pos_col = first_match(columns, ["mapinfo", "position", "pos"], ["mapinfo", "position"])

    if cpg_col is None:
        raise ValueError("S2: could not identify the CpG identifier column.")

    output: dict[str, dict[str, str]] = {}
    for _, row in frame.iterrows():
        cpg = clean_text(row[cpg_col])
        if not cpg.lower().startswith("cg"):
            continue
        output[cpg] = {
            "chr": normalise_chr(row[chr_col]) if chr_col else "",
            "pos": normalise_position(row[pos_col]) if pos_col else "",
        }
    return output


def read_annotation(path: Path, sheet_name: str | None) -> pd.DataFrame:
    """Read optional Illumina or custom CpG annotation."""
    if path.suffix.lower() in {".xlsx", ".xls"}:
        frame = pd.read_excel(path, sheet_name=sheet_name or 0)
    else:
        separator = "," if path.suffix.lower() == ".csv" else "\t"
        frame = pd.read_csv(path, sep=separator, low_memory=False)

    columns = map_columns(frame)
    aliases = {
        "cpg": ["cpg", "cpg_id", "probe", "probe_id", "ilmnid"],
        "chr": ["chr", "chromosome"],
        "pos": ["pos", "position", "mapinfo"],
        "gene": ["gene", "gene_name", "ucsc_refgene_name"],
        "context": ["context", "genomic_context", "ucsc_refgene_group"],
        "island_relation": [
            "island_relation",
            "relation_to_ucsc_cpg_island",
            "cpg_island_relation",
        ],
    }

    selected: dict[str, str] = {}
    for target, target_aliases in aliases.items():
        selected_col = first_match(columns, target_aliases, target_aliases)
        if selected_col is not None:
            selected[target] = selected_col

    if "cpg" not in selected:
        raise ValueError(
            "Annotation file does not contain a recognisable CpG/probe column."
        )

    output = pd.DataFrame({"cpg": frame[selected["cpg"]].astype(str)})
    for target in ("chr", "pos", "gene", "context", "island_relation"):
        output[target] = frame[selected[target]] if target in selected else ""

    output["chr"] = output["chr"].map(normalise_chr)
    output["pos"] = output["pos"].map(normalise_position)
    return output.drop_duplicates("cpg", keep="first")


def extract_s6(path: Path) -> pd.DataFrame:
    """Extract all q < 0.05 CpGs reported in Supplementary Table S6."""
    frame = read_excel_table(
        path,
        "S6",
        required_terms=("cpg", "spearman", "p_value", "q_value"),
    )
    columns = map_columns(frame)

    cpg_col = first_match(columns, ["cpg", "cpg_id", "probe_id"], ["cpg"])
    rho_col = first_match(
        columns,
        ["spearman_r", "spearman_rho", "rho", "correlation"],
        ["spearman", "rho"],
    )
    p_col = first_match(columns, ["p_value", "pvalue", "p"], ["p_value", "pvalue"])
    q_col = first_match(columns, ["q_value", "qvalue", "fdr"], ["q_value", "qvalue"])
    buccal_mqtl_col = first_match(columns, ["buccal_mqtl"], ["buccal_mqtl"])
    blood_mqtl_col = first_match(columns, ["blood_mqtl"], ["blood_mqtl"])
    brain_mqtl_col = first_match(columns, ["brain_mqtl"], ["brain_mqtl"])
    corsiv_col = first_match(columns, ["corsiv"], ["corsiv"])
    braun_col = first_match(columns, ["braun_et_al"], ["braun"])

    if cpg_col is None or rho_col is None:
        raise ValueError(
            "S6: could not identify the required CpG and Spearman correlation columns. "
            f"Observed columns: {list(frame.columns)}"
        )

    records: list[dict[str, Any]] = []
    for _, row in frame.iterrows():
        cpg = clean_text(row[cpg_col])
        rho = numeric(row[rho_col])
        if not cpg.lower().startswith("cg") or rho is None:
            continue

        p_value = numeric(row[p_col]) if p_col else None
        q_value = numeric(row[q_col]) if q_col else None

        records.append(
            {
                "cpg": cpg,
                "corr_value": rho,
                "abs_corr": abs(rho),
                "p_value": "" if p_value is None else p_value,
                "q_value": "" if q_value is None else q_value,
                "buccal_mqtl": yes_no(row[buccal_mqtl_col]) if buccal_mqtl_col else "",
                "blood_mqtl": yes_no(row[blood_mqtl_col]) if blood_mqtl_col else "",
                "brain_mqtl": yes_no(row[brain_mqtl_col]) if brain_mqtl_col else "",
                "corsiv": yes_no(row[corsiv_col]) if corsiv_col else "",
                "braun_et_al": yes_no(row[braun_col]) if braun_col else "",
            }
        )

    return pd.DataFrame.from_records(records)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract Sommerer et al. (2022) Supplementary Table S6 "
            "buccal–prefrontal cortex correlations."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Official Sommerer supplementary XLSX workbook containing S2 and S6.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("sommerer_buccal_brain_correlations.csv"),
        help="Output harmonised CSV.",
    )
    parser.add_argument(
        "--annotation",
        type=Path,
        default=None,
        help="Optional CpG annotation CSV/TSV/XLSX.",
    )
    parser.add_argument(
        "--annotation-sheet",
        default=None,
        help="Optional worksheet name for an Excel annotation file.",
    )
    parser.add_argument(
        "--min-abs-corr",
        type=float,
        default=None,
        help=(
            "Optional inclusive absolute-correlation threshold. "
            "No additional threshold is applied by default."
        ),
    )
    parser.add_argument(
        "--n-pairs",
        type=int,
        default=None,
        help=(
            "Optional study-level matched-pair count to record in every row. "
            "Leave unset if not being asserted from the source."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.min_abs_corr is not None and not 0 <= args.min_abs_corr <= 1:
        raise ValueError("--min-abs-corr must be between 0 and 1.")

    workbook = pd.ExcelFile(args.input)
    required_sheets = {"S2", "S6"}
    missing = required_sheets - set(workbook.sheet_names)
    if missing:
        raise ValueError(
            f"Workbook is missing required worksheets: {sorted(missing)}. "
            f"Available worksheets: {workbook.sheet_names}"
        )

    coordinates = read_s2_coordinates(args.input)
    result = extract_s6(args.input)

    if result.empty:
        raise RuntimeError("No CpG-level observations were extracted from S6.")

    result["chr"] = result["cpg"].map(
        lambda cpg: coordinates.get(cpg, {}).get("chr", "")
    )
    result["pos"] = result["cpg"].map(
        lambda cpg: coordinates.get(cpg, {}).get("pos", "")
    )
    result["gene"] = ""
    result["context"] = ""
    result["island_relation"] = ""

    if args.annotation:
        annotation = read_annotation(args.annotation, args.annotation_sheet)
        result = result.merge(
            annotation,
            on="cpg",
            how="left",
            suffixes=("", "_annotation"),
        )

        for field in ("chr", "pos", "gene", "context", "island_relation"):
            annotation_field = f"{field}_annotation"
            if annotation_field in result.columns:
                original = result[field].fillna("").astype(str)
                replacement = result[annotation_field].fillna("").astype(str)
                result[field] = original.where(original.str.len() > 0, replacement)
                result = result.drop(columns=[annotation_field])

    if args.min_abs_corr is not None:
        result = result.loc[result["abs_corr"] >= args.min_abs_corr].copy()

    result["peripheral_tissue"] = PERIPHERAL_TISSUE
    result["brain_tissue"] = BRAIN_TISSUE
    result["tissue_pair"] = f"{PERIPHERAL_TISSUE} vs {BRAIN_TISSUE}"
    result["corr_type"] = CORR_TYPE
    result["n_pairs"] = "" if args.n_pairs is None else args.n_pairs
    result["analysis_type"] = ANALYSIS_TYPE
    result["source"] = SOURCE
    result["source_doi"] = SOURCE_DOI
    result["source_url"] = SOURCE_URL
    result["database_tool"] = DATABASE_TOOL
    result["source_table"] = "Supplementary Table S6"
    result["source_scope"] = SOURCE_SCOPE

    result = result[OUTPUT_COLUMNS].sort_values("cpg", kind="stable")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)

    print(
        f"Wrote {len(result):,} Sommerer CpG observations to {args.output}",
        file=sys.stderr,
    )
    print(
        f"Coordinates found for "
        f"{(result['chr'].astype(str).str.len() > 0).sum():,} observations.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
