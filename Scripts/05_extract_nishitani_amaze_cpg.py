#!/usr/bin/env python3
"""
05_extract_nishitani_amaze_cpg.py

Extract published CpG-level brain–peripheral DNA methylation correlations from:

Nishitani S, Isozaki M, Yao A, et al.
Cross-tissue correlations of genome-wide DNA methylation in Japanese live
human brain and blood, saliva, and buccal epithelial tissues.
Translational Psychiatry. 2023;13:72.
DOI: 10.1038/s41398-023-02370-0

Primary input
-------------
The script is designed for Supplementary Table S4:

    41398_2023_2370_MOESM3_ESM.xlsx

This workbook contains six worksheets:

Raw-correlation significance sets
    BRvsBL_BH
    BRvsSA_BH
    BRvsEP_BH

Cell-proportion-adjusted significance sets
    BRvsBLadj_BH
    BRvsSAadj_BH
    BRvsEPadj_BH

nb. these correspond to:
    Brain ↔ whole blood
    Brain ↔ saliva
    Brain ↔ buccal epithelium

Each worksheet includes both raw and adjusted estimates, but the membership of
the worksheet is determined by whether the raw or adjusted analysis surpassed
the Benjamini–Hochberg significance criterion. To avoid duplicate observations:

- raw estimates are extracted only from BRvsBL_BH, BRvsSA_BH and BRvsEP_BH;
- adjusted estimates are extracted only from BRvsBLadj_BH, BRvsSAadj_BH and
  BRvsEPadj_BH.

By default, both raw and adjusted results are retained as distinct observations.
Use --mode raw or --mode adjusted to restrict extraction.

Important scope limitation
--------------------------
Supplementary Table S4 contains only CpGs surpassing the study's
Benjamini–Hochberg significance threshold. It is not the complete genome-wide
AMAZE-CpG correlation dataset.

The R scripts supplied with the article for GSE59685 and GSE95049 reproduce
pre-processing of external validation datasets. They are not required to extract
the primary Nishitani/AMAZE-CpG results from Supplementary Table S4.

No additional absolute-correlation threshold is applied by default.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Any, Mapping

import pandas as pd


SOURCE = "Nishitani et al. (2023), Translational Psychiatry 13:72"
SOURCE_DOI = "10.1038/s41398-023-02370-0"
SOURCE_URL = "https://snishit-amaze-cpg.web.app/"
DATABASE_TOOL = "AMAZE-CpG"
BRAIN_TISSUE = "Living neurosurgical brain tissue"
CORR_TYPE = "Spearman rho"

SHEET_SPECS = {
    "BRvsBL_BH": {
        "peripheral_tissue": "Whole blood",
        "estimate": "raw",
        "corr_column": "rho_brain_blood",
        "p_column": "p_rho_brain_blood",
    },
    "BRvsSA_BH": {
        "peripheral_tissue": "Saliva",
        "estimate": "raw",
        "corr_column": "rho_brain_saliva",
        "p_column": "p_rho_brain_saliva",
    },
    "BRvsEP_BH": {
        "peripheral_tissue": "Buccal epithelium",
        "estimate": "raw",
        "corr_column": "rho_brain_buccal",
        "p_column": "p_rho_brain_buccal",
    },
    "BRvsBLadj_BH": {
        "peripheral_tissue": "Whole blood",
        "estimate": "adjusted",
        "corr_column": "rho_brain_blood_adj",
        "p_column": "p_rho_brain_blood_adj",
    },
    "BRvsSAadj_BH": {
        "peripheral_tissue": "Saliva",
        "estimate": "adjusted",
        "corr_column": "rho_brain_saliva_adj",
        "p_column": "p_rho_brain_saliva_adj",
    },
    "BRvsEPadj_BH": {
        "peripheral_tissue": "Buccal epithelium",
        "estimate": "adjusted",
        "corr_column": "rho_brain_buccal_adj",
        "p_column": "p_rho_brain_buccal_adj",
    },
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
    "cell_composition_adjusted",
    "snp_flag",
    "mqtl_flag",
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


def boolean_text(value: Any) -> str:
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    text = clean_text(value).lower()
    if text in {"true", "t", "1", "yes", "y"}:
        return "TRUE"
    if text in {"false", "f", "0", "no", "n"}:
        return "FALSE"
    return clean_text(value)


def find_header_row(raw: pd.DataFrame) -> int:
    terms = ("chr", "mapinfo", "gene_name", "rho_brain", "snp_flag", "mqtl_flag")
    best_index = 0
    best_score = -1

    for index in range(min(20, len(raw))):
        cells = [canonical(value) for value in raw.iloc[index].tolist()]
        score = sum(any(term in cell for cell in cells) for term in terms)
        if score > best_score:
            best_index = index
            best_score = score

    return best_index


def read_sheet(path: Path, sheet_name: str) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name=sheet_name, header=None)
    header_row = find_header_row(raw)

    headers = []
    for index, value in enumerate(raw.iloc[header_row].tolist()):
        header = clean_text(value)
        if not header:
            header = "cpg" if index == 0 else f"unnamed_{index}"
        headers.append(header)

    frame = raw.iloc[header_row + 1 :].copy()
    frame.columns = headers
    return frame.dropna(how="all")


def columns_by_canonical_name(frame: pd.DataFrame) -> dict[str, str]:
    return {canonical(column): column for column in frame.columns}


def required_column(
    columns: Mapping[str, str],
    canonical_name: str,
    *,
    sheet_name: str,
) -> str:
    if canonical_name not in columns:
        raise ValueError(
            f"{sheet_name}: required column {canonical_name!r} was not found. "
            f"Observed columns: {list(columns.values())}"
        )
    return columns[canonical_name]


def extract_sheet(
    path: Path,
    sheet_name: str,
    *,
    min_abs_corr: float | None,
    n_pairs: int | None,
) -> list[dict[str, Any]]:
    spec = SHEET_SPECS[sheet_name]
    frame = read_sheet(path, sheet_name)
    columns = columns_by_canonical_name(frame)

    cpg_col = required_column(columns, "cpg", sheet_name=sheet_name)
    chr_col = required_column(columns, "chr", sheet_name=sheet_name)
    pos_col = required_column(columns, "mapinfo", sheet_name=sheet_name)
    gene_col = required_column(columns, "gene_name", sheet_name=sheet_name)
    corr_col = required_column(
        columns,
        canonical(spec["corr_column"]),
        sheet_name=sheet_name,
    )
    p_col = required_column(
        columns,
        canonical(spec["p_column"]),
        sheet_name=sheet_name,
    )
    snp_col = required_column(columns, "snp_flag", sheet_name=sheet_name)
    mqtl_col = required_column(columns, "mqtl_flag", sheet_name=sheet_name)
    bh_col = columns.get("bh")

    adjusted = spec["estimate"] == "adjusted"
    analysis_type = (
        "Published BH-significant cell-proportion-adjusted within-subject correlation"
        if adjusted
        else "Published BH-significant raw within-subject correlation"
    )
    source_scope = (
        f"CpGs included in Nishitani et al. (2023) Supplementary Table S4 "
        f"worksheet {sheet_name}, representing the "
        f"{'cell-proportion-adjusted' if adjusted else 'raw'} "
        f"brain–{spec['peripheral_tissue'].lower()} analysis surpassing the "
        f"Benjamini–Hochberg significance threshold; not the complete "
        f"genome-wide AMAZE-CpG dataset"
    )

    records: list[dict[str, Any]] = []

    for _, row in frame.iterrows():
        cpg = clean_text(row[cpg_col])
        rho = numeric(row[corr_col])

        if not cpg.lower().startswith("cg") or rho is None:
            continue
        if min_abs_corr is not None and abs(rho) < min_abs_corr:
            continue

        p_value = numeric(row[p_col])
        bh_value = numeric(row[bh_col]) if bh_col else None
        peripheral = spec["peripheral_tissue"]

        records.append(
            {
                "cpg": cpg,
                "chr": normalise_chr(row[chr_col]),
                "pos": normalise_position(row[pos_col]),
                "gene": clean_text(row[gene_col]),
                "context": "",
                "island_relation": "",
                "peripheral_tissue": peripheral,
                "brain_tissue": BRAIN_TISSUE,
                "tissue_pair": f"{peripheral} vs {BRAIN_TISSUE}",
                "corr_type": CORR_TYPE,
                "corr_value": rho,
                "abs_corr": abs(rho),
                "p_value": "" if p_value is None else p_value,
                "q_value": "" if bh_value is None else bh_value,
                "n_pairs": "" if n_pairs is None else n_pairs,
                "analysis_type": analysis_type,
                "cell_composition_adjusted": "TRUE" if adjusted else "FALSE",
                "snp_flag": boolean_text(row[snp_col]),
                "mqtl_flag": boolean_text(row[mqtl_col]),
                "source": SOURCE,
                "source_doi": SOURCE_DOI,
                "source_url": SOURCE_URL,
                "database_tool": DATABASE_TOOL,
                "source_table": f"Supplementary Table S4: {sheet_name}",
                "source_scope": source_scope,
            }
        )

    return records


def deduplicate(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    seen: set[tuple[Any, ...]] = set()
    output: list[dict[str, Any]] = []

    for record in records:
        key = (
            record["cpg"],
            record["peripheral_tissue"],
            record["cell_composition_adjusted"],
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
            "Extract Nishitani et al. (2023) AMAZE-CpG correlations from "
            "Supplementary Table S4."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Supplementary Table S4 workbook (MOESM3 XLSX).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("nishitani_amaze_cpg_correlations.csv"),
        help="Output harmonised CSV.",
    )
    parser.add_argument(
        "--mode",
        choices=("raw", "adjusted", "both"),
        default="both",
        help="Extract raw estimates, adjusted estimates, or both.",
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
            "Optional matched-pair count to record in every row. "
            "Leave unset unless being asserted from the source."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.min_abs_corr is not None and not 0 <= args.min_abs_corr <= 1:
        raise ValueError("--min-abs-corr must be between 0 and 1.")
    if args.n_pairs is not None and args.n_pairs < 3:
        raise ValueError("--n-pairs must be at least 3.")

    workbook = pd.ExcelFile(args.input)
    missing = set(SHEET_SPECS) - set(workbook.sheet_names)
    if missing:
        raise ValueError(
            f"Workbook is missing required worksheets: {sorted(missing)}. "
            f"Available worksheets: {workbook.sheet_names}"
        )

    if args.mode == "raw":
        sheets = [name for name, spec in SHEET_SPECS.items() if spec["estimate"] == "raw"]
    elif args.mode == "adjusted":
        sheets = [
            name for name, spec in SHEET_SPECS.items() if spec["estimate"] == "adjusted"
        ]
    else:
        sheets = list(SHEET_SPECS)

    records: list[dict[str, Any]] = []
    for sheet_name in sheets:
        sheet_records = extract_sheet(
            args.input,
            sheet_name,
            min_abs_corr=args.min_abs_corr,
            n_pairs=args.n_pairs,
        )
        records.extend(sheet_records)
        print(
            f"{sheet_name}: extracted {len(sheet_records):,} observations",
            file=sys.stderr,
        )

    records = deduplicate(records)
    if not records:
        raise RuntimeError("No CpG-level observations were extracted.")

    result = pd.DataFrame.from_records(records, columns=OUTPUT_COLUMNS)
    result = result.sort_values(
        ["cell_composition_adjusted", "peripheral_tissue", "cpg"],
        kind="stable",
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)

    print(f"Wrote {len(result):,} observations to {args.output}", file=sys.stderr)
    print("Counts by analysis and tissue:", file=sys.stderr)
    counts = result.groupby(
        ["cell_composition_adjusted", "peripheral_tissue"]
    ).size()
    for (adjusted, tissue), count in counts.items():
        label = "adjusted" if adjusted == "TRUE" else "raw"
        print(f"  {label}, {tissue}: {count:,}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
