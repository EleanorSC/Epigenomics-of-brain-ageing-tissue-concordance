#!/usr/bin/env python3
"""
03_extract_braun_image_cpg.py

Extract published CpG-level brain–peripheral DNA methylation correlations from
Braun et al. (2019) Supplementary Table 3 and convert them to a harmonised
long-format CSV for the CpG brain–peripheral concordance database.

Reference
---------
Braun PR et al. Genome-wide DNA methylation comparison between live human
brain and peripheral tissues within individuals. Transl Psychiatry. 2019;9:47.
DOI: 10.1038/s41398-019-0376-y

Scope
-----
Supplementary Table 3 contains CpGs passing the study's Bonferroni threshold
for within-subject Spearman correlations between live brain tissue and:

- whole blood
- saliva
- buccal epithelium

The script does not reconstruct all IMAGE-CpG correlations from GSE111165.
It extracts only the published significant subset present in Supplementary
Table 3.

Because supplementary-workbook column names can vary slightly across downloaded
versions, the parser uses flexible column matching and supports two common layouts:

1. Separate worksheets for blood, saliva and buccal results.
2. One worksheet containing separate correlation/P-value columns for each tissue.

Run `python 03_extract_braun_image_cpg.py --inspect TABLE.xlsx` first to print
worksheet names and headers before extraction.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import pandas as pd


SOURCE = "Braun et al. (2019), Translational Psychiatry 9:47"
SOURCE_DOI = "10.1038/s41398-019-0376-y"
SOURCE_URL = "https://han-lab.org/methylation/default/imageCpG"
DATABASE_TOOL = "IMAGE-CpG"
BRAIN_TISSUE = "Live resected brain tissue"
CORR_TYPE = "Spearman rho"
ANALYSIS_TYPE = "Published Bonferroni-significant within-subject correlation"
SOURCE_SCOPE = (
    "CpGs reported in Braun et al. (2019) Supplementary Table 3 after "
    "Bonferroni correction; not the complete IMAGE-CpG correlation dataset"
)

TISSUES = {
    "blood": "Whole blood",
    "saliva": "Saliva",
    "buccal": "Buccal epithelium",
}

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
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
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
    text = clean_text(value)
    return re.sub(r"^chr", "", text, flags=re.IGNORECASE)


def normalise_position(value: Any) -> str:
    number = numeric(value)
    if number is not None and number.is_integer():
        return str(int(number))
    return clean_text(value)


def first_matching_column(
    columns: Mapping[str, str],
    exact_aliases: Sequence[str],
    contains_aliases: Sequence[str] = (),
) -> str | None:
    for alias in exact_aliases:
        if alias in columns:
            return columns[alias]
    for alias in contains_aliases:
        for canonical_name, original_name in columns.items():
            if alias in canonical_name:
                return original_name
    return None


def find_header_row(raw: pd.DataFrame, max_rows: int = 30) -> int:
    """Locate the likely header row in a worksheet with title/notes above it."""
    terms = (
        "cpg",
        "probe",
        "rho",
        "correlation",
        "p_value",
        "pvalue",
        "chromosome",
        "mapinfo",
        "gene",
    )
    best_index = 0
    best_score = -1

    for index in range(min(max_rows, len(raw))):
        cells = [canonical(value) for value in raw.iloc[index].tolist()]
        score = sum(any(term in cell for cell in cells) for term in terms)
        if score > best_score:
            best_score = score
            best_index = index

    return best_index


def read_sheet(path: Path, sheet_name: str) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name=sheet_name, header=None)
    header_row = find_header_row(raw)
    headers = [
        clean_text(value) or f"unnamed_{index}"
        for index, value in enumerate(raw.iloc[header_row].tolist())
    ]
    frame = raw.iloc[header_row + 1 :].copy()
    frame.columns = headers
    frame = frame.dropna(how="all")
    return frame


def inspect_workbook(path: Path) -> None:
    workbook = pd.ExcelFile(path)
    print(f"Workbook: {path}")
    for sheet in workbook.sheet_names:
        frame = read_sheet(path, sheet)
        print(f"\n[{sheet}]")
        print("Columns:")
        for column in frame.columns:
            print(f"  - {column}")


def infer_tissue_from_text(value: str) -> str | None:
    text = canonical(value)
    if "blood" in text:
        return "blood"
    if "saliva" in text:
        return "saliva"
    if "buccal" in text or "epithelial" in text:
        return "buccal"
    return None


def base_columns(frame: pd.DataFrame) -> dict[str, str | None]:
    columns = {canonical(column): column for column in frame.columns}

    return {
        "cpg": first_matching_column(
            columns,
            ["cpg", "cpg_id", "probe", "probe_id", "ilmnid"],
            ["cpg", "probe"],
        ),
        "chr": first_matching_column(
            columns,
            ["chr", "chromosome"],
            ["chromosome", "chr"],
        ),
        "pos": first_matching_column(
            columns,
            ["pos", "position", "mapinfo", "genomic_coordinate"],
            ["mapinfo", "position", "coordinate"],
        ),
        "gene": first_matching_column(
            columns,
            ["gene", "gene_name", "associated_gene", "ucsc_refgene_name"],
            ["gene"],
        ),
        "context": first_matching_column(
            columns,
            ["context", "genic_region", "gene_region", "ucsc_refgene_group"],
            ["context", "genic", "gene_region", "refgene_group"],
        ),
        "island_relation": first_matching_column(
            columns,
            [
                "island_relation",
                "relation_to_ucsc_cpg_island",
                "cpg_island_relation",
            ],
            ["island"],
        ),
        "n_pairs": first_matching_column(
            columns,
            ["n", "n_pairs", "sample_size"],
            ["sample_size", "n_pairs"],
        ),
    }


def find_tissue_columns(
    frame: pd.DataFrame,
    tissue_key: str,
    *,
    allow_generic: bool,
) -> dict[str, str | None]:
    columns = {canonical(column): column for column in frame.columns}
    tissue_aliases = {
        "blood": ("blood", "bl"),
        "saliva": ("saliva", "sa"),
        "buccal": ("buccal", "bu", "ep"),
    }[tissue_key]

    corr_candidates: list[str] = []
    p_candidates: list[str] = []
    q_candidates: list[str] = []

    for canonical_name, original_name in columns.items():
        tissue_match = any(
            re.search(rf"(^|_){re.escape(alias)}($|_)", canonical_name)
            for alias in tissue_aliases
        )
        if not tissue_match:
            continue

        if any(token in canonical_name for token in ("rho", "corr", "correlation")):
            corr_candidates.append(original_name)
        if any(token in canonical_name for token in ("p_value", "pvalue", "_p", "p_")):
            p_candidates.append(original_name)
        if any(token in canonical_name for token in ("q_value", "qvalue", "fdr", "adjusted_p")):
            q_candidates.append(original_name)

    corr_col = corr_candidates[0] if corr_candidates else None
    p_col = p_candidates[0] if p_candidates else None
    q_col = q_candidates[0] if q_candidates else None

    if allow_generic and corr_col is None:
        corr_col = first_matching_column(
            columns,
            ["rho", "spearman_rho", "correlation", "corr", "r"],
            ["rho", "correlation"],
        )
    if allow_generic and p_col is None:
        p_col = first_matching_column(
            columns,
            ["p", "p_value", "pvalue"],
            ["p_value", "pvalue"],
        )
    if allow_generic and q_col is None:
        q_col = first_matching_column(
            columns,
            ["q", "q_value", "qvalue", "fdr", "adjusted_p"],
            ["q_value", "qvalue", "fdr"],
        )

    return {"corr": corr_col, "p": p_col, "q": q_col}


def make_records(
    frame: pd.DataFrame,
    *,
    sheet_name: str,
    tissue_key: str,
    corr_column: str,
    p_column: str | None,
    q_column: str | None,
) -> list[dict[str, Any]]:
    common = base_columns(frame)
    if common["cpg"] is None:
        raise ValueError(
            f"{sheet_name}: could not identify a CpG/probe identifier column. "
            f"Observed columns: {list(frame.columns)}"
        )

    records: list[dict[str, Any]] = []
    peripheral = TISSUES[tissue_key]

    for _, row in frame.iterrows():
        cpg = clean_text(row[common["cpg"]])
        rho = numeric(row[corr_column])

        if not cpg or not cpg.lower().startswith("cg") or rho is None:
            continue

        p_value = numeric(row[p_column]) if p_column else None
        q_value = numeric(row[q_column]) if q_column else None
        n_pairs_value = numeric(row[common["n_pairs"]]) if common["n_pairs"] else None

        records.append(
            {
                "cpg": cpg,
                "chr": normalise_chr(row[common["chr"]]) if common["chr"] else "",
                "pos": normalise_position(row[common["pos"]]) if common["pos"] else "",
                "gene": clean_text(row[common["gene"]]) if common["gene"] else "",
                "context": clean_text(row[common["context"]]) if common["context"] else "",
                "island_relation": (
                    clean_text(row[common["island_relation"]])
                    if common["island_relation"]
                    else ""
                ),
                "peripheral_tissue": peripheral,
                "brain_tissue": BRAIN_TISSUE,
                "tissue_pair": f"{peripheral} vs {BRAIN_TISSUE}",
                "corr_type": CORR_TYPE,
                "corr_value": rho,
                "abs_corr": abs(rho),
                "p_value": "" if p_value is None else p_value,
                "q_value": "" if q_value is None else q_value,
                "n_pairs": "" if n_pairs_value is None else int(n_pairs_value),
                "analysis_type": ANALYSIS_TYPE,
                "source": SOURCE,
                "source_doi": SOURCE_DOI,
                "source_url": SOURCE_URL,
                "database_tool": DATABASE_TOOL,
                "source_table": f"Supplementary Table 3: {sheet_name}",
                "source_scope": SOURCE_SCOPE,
            }
        )

    return records


def extract_from_separate_sheets(
    path: Path,
    sheet_names: Sequence[str],
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []

    for sheet_name in sheet_names:
        tissue_key = infer_tissue_from_text(sheet_name)
        if tissue_key is None:
            continue

        frame = read_sheet(path, sheet_name)
        tissue_columns = find_tissue_columns(
            frame,
            tissue_key,
            allow_generic=True,
        )
        if tissue_columns["corr"] is None:
            raise ValueError(
                f"{sheet_name}: could not identify a correlation column. "
                f"Observed columns: {list(frame.columns)}"
            )

        records.extend(
            make_records(
                frame,
                sheet_name=sheet_name,
                tissue_key=tissue_key,
                corr_column=tissue_columns["corr"],
                p_column=tissue_columns["p"],
                q_column=tissue_columns["q"],
            )
        )

    return records


def extract_from_combined_sheet(
    path: Path,
    sheet_name: str,
) -> list[dict[str, Any]]:
    frame = read_sheet(path, sheet_name)
    records: list[dict[str, Any]] = []

    for tissue_key in TISSUES:
        tissue_columns = find_tissue_columns(
            frame,
            tissue_key,
            allow_generic=False,
        )
        if tissue_columns["corr"] is None:
            continue

        records.extend(
            make_records(
                frame,
                sheet_name=sheet_name,
                tissue_key=tissue_key,
                corr_column=tissue_columns["corr"],
                p_column=tissue_columns["p"],
                q_column=tissue_columns["q"],
            )
        )

    return records


def deduplicate(records: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    seen: set[tuple[Any, ...]] = set()
    output: list[dict[str, Any]] = []

    for record in records:
        key = (
            record["cpg"],
            record["peripheral_tissue"],
            record["brain_tissue"],
            record["corr_value"],
            record["source_table"],
        )
        if key not in seen:
            seen.add(key)
            output.append(record)

    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract Braun et al. (2019) Supplementary Table 3 "
            "CpG-level brain–peripheral correlations."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        help="Braun et al. Supplementary Table 3 XLSX workbook.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("braun_image_cpg_significant_correlations.csv"),
        help="Output harmonised CSV.",
    )
    parser.add_argument(
        "--sheet",
        action="append",
        default=[],
        help=(
            "Worksheet to process. May be supplied repeatedly. If omitted, all "
            "worksheets are inspected automatically."
        ),
    )
    parser.add_argument(
        "--combined-sheet",
        default=None,
        help=(
            "Name of a single worksheet containing distinct blood, saliva and "
            "buccal correlation columns."
        ),
    )
    parser.add_argument(
        "--inspect",
        type=Path,
        default=None,
        help="Print worksheet names and detected headers, then exit.",
    )
    parser.add_argument(
        "--min-abs-corr",
        type=float,
        default=None,
        help=(
            "Optional inclusive absolute-correlation filter. No additional "
            "threshold is applied by default."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.inspect:
        inspect_workbook(args.inspect)
        return 0

    if args.input is None:
        raise ValueError("--input is required unless --inspect is used.")

    if args.min_abs_corr is not None and not 0 <= args.min_abs_corr <= 1:
        raise ValueError("--min-abs-corr must be between 0 and 1.")

    workbook = pd.ExcelFile(args.input)
    sheets = args.sheet or workbook.sheet_names

    if args.combined_sheet:
        if args.combined_sheet not in workbook.sheet_names:
            raise ValueError(
                f"Worksheet {args.combined_sheet!r} was not found. "
                f"Available sheets: {workbook.sheet_names}"
            )
        records = extract_from_combined_sheet(args.input, args.combined_sheet)
    else:
        # First try worksheets whose names identify a peripheral tissue.
        records = extract_from_separate_sheets(args.input, sheets)

        # If none were found, try every worksheet as a combined table.
        if not records:
            for sheet_name in sheets:
                records.extend(extract_from_combined_sheet(args.input, sheet_name))

    records = deduplicate(records)

    if args.min_abs_corr is not None:
        records = [
            record
            for record in records
            if float(record["abs_corr"]) >= args.min_abs_corr
        ]

    if not records:
        raise RuntimeError(
            "No CpG-level correlations were extracted. Run the script with "
            "--inspect first and check the worksheet/header structure."
        )

    result = pd.DataFrame.from_records(records, columns=OUTPUT_COLUMNS)
    result = result.sort_values(
        ["peripheral_tissue", "cpg"],
        kind="stable",
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)

    counts = result.groupby("peripheral_tissue").size()
    print(f"Wrote {len(result):,} observations to {args.output}", file=sys.stderr)
    for tissue, count in counts.items():
        print(f"  {tissue}: {count:,}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
