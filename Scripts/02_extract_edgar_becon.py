#!/usr/bin/env python3
"""
02_extract_edgar_becon.py

Reconstruct CpG-level blood–brain DNA methylation correlations for the
BECon study (Edgar et al., 2017) from a local GSE95049 processed matrix.

The script:

1. reads a processed CpG × sample matrix from CSV/TSV/TXT, optionally gzipped;
2. identifies BECon blood, BA7, BA10 and BA20 samples from their names;
3. removes known technical replicate columns by default;
4. aligns matched subjects across blood and each brain region;
5. calculates CpG-level Spearman correlations;
6. optionally joins Illumina probe annotation;
7. writes a harmonised long-format CSV for downstream database construction.

This script does not reproduce every historical BECon preprocessing step
(BMIQ, ComBat and cell-composition residualisation). Its default output should
therefore be labelled:

    "Recomputed from public processed GSE95049 data"

Reference:
Edgar RD et al. BECon: a tool for interpreting DNA methylation findings
from blood in the context of brain. Transl Psychiatry. 2017;7:e1187.
DOI: 10.1038/tp.2017.171
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


SOURCE_STUDY = "Edgar et al. (2017), Translational Psychiatry 7:e1187"
SOURCE_DOI = "10.1038/tp.2017.171"
SOURCE_URL = "https://redgar598.shinyapps.io/BECon/"
DATABASE_TOOL = "BECon"
ANALYSIS_TYPE = "Recomputed from public processed GSE95049 data"

DEFAULT_REPLICATE_COLUMNS = {
    "BA7250rep",
    "BLOOD169rep",
    "PBMC169rep",
    "BLOOD250rep",
    "PBMC250rep",
    "BA102482",
}

TISSUE_LABELS = {
    "blood": "Whole blood",
    "ba7": "Brodmann area 7",
    "ba10": "Brodmann area 10",
    "ba20": "Brodmann area 20",
}


@dataclass(frozen=True)
class SampleInfo:
    original_name: str
    tissue_code: str
    subject_id: str
    is_replicate: bool = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reconstruct BECon CpG-level blood–brain correlations from GSE95049."
    )
    parser.add_argument(
        "--matrix",
        type=Path,
        required=True,
        help="Processed CpG × sample matrix (.csv, .tsv, .txt, optionally .gz).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("edgar_becon_cpg_correlations.csv"),
        help="Output harmonised CSV.",
    )
    parser.add_argument(
        "--annotation",
        type=Path,
        default=None,
        help=(
            "Optional probe annotation CSV/TSV/XLSX. Expected columns are matched "
            "flexibly to CpG ID, chromosome, position, gene, genomic context and "
            "CpG-island relation."
        ),
    )
    parser.add_argument(
        "--annotation-sheet",
        default=None,
        help="Optional worksheet name when --annotation is an Excel workbook.",
    )
    parser.add_argument(
        "--probe-column",
        default=None,
        help="Name of probe/CpG column in the matrix. Defaults to the first column.",
    )
    parser.add_argument(
        "--separator",
        default=None,
        help="Optional matrix delimiter. Defaults to automatic detection.",
    )
    parser.add_argument(
        "--matrix-kind",
        choices=("beta", "signal-pairs"),
        default="beta",
        help=(
            "'beta' expects one value per CpG/sample. 'signal-pairs' expects "
            "methylated/unmethylated column pairs and requires --methylated-suffix "
            "and --unmethylated-suffix."
        ),
    )
    parser.add_argument(
        "--methylated-suffix",
        default="_M",
        help="Suffix identifying methylated-intensity columns in signal-pairs mode.",
    )
    parser.add_argument(
        "--unmethylated-suffix",
        default="_U",
        help="Suffix identifying unmethylated-intensity columns in signal-pairs mode.",
    )
    parser.add_argument(
        "--beta-offset",
        type=float,
        default=100.0,
        help="Offset used for beta = M / (M + U + offset).",
    )
    parser.add_argument(
        "--min-pairs",
        type=int,
        default=3,
        help="Minimum number of complete subject pairs required for a correlation.",
    )
    parser.add_argument(
        "--min-abs-corr",
        type=float,
        default=None,
        help=(
            "Optional inclusive absolute-correlation threshold. Default retains all "
            "finite correlations."
        ),
    )
    parser.add_argument(
        "--drop-known-replicates",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Remove known BECon technical replicate columns.",
    )
    parser.add_argument(
        "--replicate-column",
        action="append",
        default=[],
        help="Additional sample column to exclude; may be supplied repeatedly.",
    )
    parser.add_argument(
        "--sample-map",
        type=Path,
        default=None,
        help=(
            "Optional CSV with columns sample,tissue,subject_id,is_replicate. "
            "Use this when sample names cannot be parsed automatically."
        ),
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=20000,
        help="Number of CpGs processed per chunk.",
    )
    return parser.parse_args()


def canonical(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def infer_separator(path: Path, explicit: str | None) -> str:
    if explicit is not None:
        return explicit
    suffixes = "".join(path.suffixes).lower()
    if ".csv" in suffixes:
        return ","
    return "\t"


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8-sig") if path.suffix.lower() == ".gz" else path.open(
        "rt", encoding="utf-8-sig"
    )


def read_header(path: Path, separator: str) -> list[str]:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter=separator)
        try:
            return [str(x).strip() for x in next(reader)]
        except StopIteration as exc:
            raise ValueError(f"{path} is empty.") from exc


def parse_sample_name(name: str) -> SampleInfo | None:
    """
    Parse common GSE95049/BECon sample names.

    Handles examples such as:
      BLOOD101
      PBMC101
      BA7101
      BA10101
      BA20101
      BA7250rep
      BA102481 / BA102482  -> subject 248, replicate indicator
    """
    clean = re.sub(r"[\s._-]+", "", str(name)).upper()

    is_replicate = clean.endswith("REP")
    clean_without_rep = re.sub(r"REP$", "", clean)

    if clean_without_rep.startswith(("BLOOD", "PBMC")):
        subject = re.sub(r"^(BLOOD|PBMC)", "", clean_without_rep)
        if subject.isdigit():
            return SampleInfo(name, "blood", subject, is_replicate)
        return None

    # Explicit historical BA10 duplicate names for subject 248.
    if clean_without_rep in {"BA102481", "BA102482"}:
        return SampleInfo(
            name,
            "ba10",
            "248",
            is_replicate=(clean_without_rep == "BA102482"),
        )

    match = re.match(r"^(BA7|BA10|BA20)(\d+)$", clean_without_rep)
    if not match:
        return None

    tissue_code = match.group(1).lower()
    subject_id = match.group(2)
    return SampleInfo(name, tissue_code, subject_id, is_replicate)


def load_sample_map(path: Path) -> dict[str, SampleInfo]:
    frame = pd.read_csv(path)
    columns = {canonical(c): c for c in frame.columns}
    required = {"sample", "tissue", "subject_id"}
    missing = required - set(columns)
    if missing:
        raise ValueError(
            f"Sample map is missing required columns: {sorted(missing)}. "
            f"Observed columns: {list(frame.columns)}"
        )

    output: dict[str, SampleInfo] = {}
    replicate_col = columns.get("is_replicate")
    for _, row in frame.iterrows():
        tissue = canonical(row[columns["tissue"]])
        aliases = {
            "whole_blood": "blood",
            "blood": "blood",
            "pbmc": "blood",
            "ba7": "ba7",
            "brodmann_area_7": "ba7",
            "ba10": "ba10",
            "brodmann_area_10": "ba10",
            "ba20": "ba20",
            "brodmann_area_20": "ba20",
        }
        if tissue not in aliases:
            raise ValueError(f"Unsupported tissue label in sample map: {tissue}")
        sample = str(row[columns["sample"]]).strip()
        output[sample] = SampleInfo(
            original_name=sample,
            tissue_code=aliases[tissue],
            subject_id=str(row[columns["subject_id"]]).strip(),
            is_replicate=bool(row[replicate_col]) if replicate_col else False,
        )
    return output


def build_sample_index(
    sample_columns: Sequence[str],
    *,
    sample_map: Mapping[str, SampleInfo] | None,
    excluded_columns: set[str],
) -> dict[str, dict[str, str]]:
    tissue_to_subject: dict[str, dict[str, str]] = {
        "blood": {},
        "ba7": {},
        "ba10": {},
        "ba20": {},
    }

    for column in sample_columns:
        if column in excluded_columns:
            continue
        info = sample_map.get(column) if sample_map is not None else parse_sample_name(column)
        if info is None:
            continue
        if info.is_replicate:
            continue

        existing = tissue_to_subject[info.tissue_code].get(info.subject_id)
        if existing is not None:
            raise ValueError(
                f"Multiple non-replicate columns map to {info.tissue_code}, "
                f"subject {info.subject_id}: {existing!r} and {column!r}. "
                "Use --sample-map or --replicate-column to resolve this."
            )
        tissue_to_subject[info.tissue_code][info.subject_id] = column

    return tissue_to_subject


def convert_signal_pairs_to_beta(
    frame: pd.DataFrame,
    *,
    probe_column: str,
    methylated_suffix: str,
    unmethylated_suffix: str,
    offset: float,
) -> pd.DataFrame:
    columns = [c for c in frame.columns if c != probe_column]
    methylated = {
        c[: -len(methylated_suffix)]: c
        for c in columns
        if c.endswith(methylated_suffix)
    }
    unmethylated = {
        c[: -len(unmethylated_suffix)]: c
        for c in columns
        if c.endswith(unmethylated_suffix)
    }
    shared = sorted(set(methylated) & set(unmethylated))
    if not shared:
        raise ValueError(
            "No methylated/unmethylated column pairs were detected. "
            "Check --methylated-suffix and --unmethylated-suffix."
        )

    result = pd.DataFrame({probe_column: frame[probe_column]})
    for sample in shared:
        m = pd.to_numeric(frame[methylated[sample]], errors="coerce")
        u = pd.to_numeric(frame[unmethylated[sample]], errors="coerce")
        denominator = m + u + offset
        result[sample] = np.where(denominator > 0, m / denominator, np.nan)
    return result


def correlation_for_pair(
    x: np.ndarray,
    y: np.ndarray,
    *,
    min_pairs: int,
) -> tuple[float, float, int]:
    complete = np.isfinite(x) & np.isfinite(y)
    n_pairs = int(complete.sum())
    if n_pairs < min_pairs:
        return math.nan, math.nan, n_pairs

    x_complete = x[complete]
    y_complete = y[complete]

    # Spearman correlation is undefined if either vector is constant.
    if np.nanstd(x_complete) == 0 or np.nanstd(y_complete) == 0:
        return math.nan, math.nan, n_pairs

    result = spearmanr(x_complete, y_complete)
    return float(result.statistic), float(result.pvalue), n_pairs


def process_chunk(
    frame: pd.DataFrame,
    *,
    probe_column: str,
    sample_index: Mapping[str, Mapping[str, str]],
    min_pairs: int,
    min_abs_corr: float | None,
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []

    for brain_code in ("ba7", "ba10", "ba20"):
        shared_subjects = sorted(
            set(sample_index["blood"]) & set(sample_index[brain_code]),
            key=lambda value: (len(value), value),
        )
        if len(shared_subjects) < min_pairs:
            continue

        blood_columns = [sample_index["blood"][subject] for subject in shared_subjects]
        brain_columns = [sample_index[brain_code][subject] for subject in shared_subjects]

        blood_values = frame[blood_columns].apply(pd.to_numeric, errors="coerce").to_numpy(
            dtype=float
        )
        brain_values = frame[brain_columns].apply(pd.to_numeric, errors="coerce").to_numpy(
            dtype=float
        )
        probes = frame[probe_column].astype(str).to_numpy()

        for row_index, cpg in enumerate(probes):
            rho, p_value, n_pairs = correlation_for_pair(
                blood_values[row_index],
                brain_values[row_index],
                min_pairs=min_pairs,
            )
            if not math.isfinite(rho):
                continue
            if min_abs_corr is not None and abs(rho) < min_abs_corr:
                continue

            brain_label = TISSUE_LABELS[brain_code]
            output.append(
                {
                    "cpg": cpg,
                    "peripheral_tissue": TISSUE_LABELS["blood"],
                    "brain_tissue": brain_label,
                    "tissue_pair": f"{TISSUE_LABELS['blood']} vs {brain_label}",
                    "corr_type": "Spearman rho",
                    "corr_value": rho,
                    "abs_corr": abs(rho),
                    "p_value": p_value,
                    "q_value": "",
                    "n_pairs": n_pairs,
                    "analysis_type": ANALYSIS_TYPE,
                    "source": SOURCE_STUDY,
                    "source_doi": SOURCE_DOI,
                    "source_url": SOURCE_URL,
                    "database_tool": DATABASE_TOOL,
                    "source_table": "GSE95049 processed matrix",
                    "source_scope": "All finite CpG-level correlations recomputed from supplied matrix",
                }
            )

    return output


def read_annotation(path: Path, sheet: str | None) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        frame = pd.read_excel(path, sheet_name=sheet or 0)
    else:
        separator = "," if suffix == ".csv" else "\t"
        frame = pd.read_csv(path, sep=separator, low_memory=False)

    canonical_columns = {canonical(c): c for c in frame.columns}

    aliases = {
        "cpg": ["cpg", "probeid", "probe_id", "name", "ilmnid"],
        "chr": ["chr", "chromosome", "chromosome_37", "chromosome_38"],
        "pos": ["pos", "position", "mapinfo", "coordinate"],
        "gene": ["gene", "gene_name", "ucsc_refgene_name"],
        "context": ["context", "genomic_context", "ucsc_refgene_group"],
        "island_relation": [
            "island_relation",
            "relation_to_ucsc_cpg_island",
            "relation_to_ucsc_cpg_island_name",
        ],
    }

    selected: dict[str, str] = {}
    for target, candidates in aliases.items():
        for candidate in candidates:
            if candidate in canonical_columns:
                selected[target] = canonical_columns[candidate]
                break

    if "cpg" not in selected:
        raise ValueError(
            f"Could not identify a CpG/probe column in annotation file. "
            f"Observed columns: {list(frame.columns)}"
        )

    output = pd.DataFrame({"cpg": frame[selected["cpg"]].astype(str)})
    for target in ("chr", "pos", "gene", "context", "island_relation"):
        output[target] = frame[selected[target]] if target in selected else ""
    return output.drop_duplicates(subset=["cpg"], keep="first")


def main() -> int:
    args = parse_args()

    if args.min_pairs < 3:
        raise ValueError("--min-pairs must be at least 3.")
    if args.min_abs_corr is not None and not 0 <= args.min_abs_corr <= 1:
        raise ValueError("--min-abs-corr must be between 0 and 1.")

    separator = infer_separator(args.matrix, args.separator)
    header = read_header(args.matrix, separator)
    probe_column = args.probe_column or header[0]
    if probe_column not in header:
        raise ValueError(
            f"Probe column {probe_column!r} was not found. Available columns: {header[:20]}"
        )

    excluded = set(args.replicate_column)
    if args.drop_known_replicates:
        excluded |= DEFAULT_REPLICATE_COLUMNS

    sample_map = load_sample_map(args.sample_map) if args.sample_map else None

    # In signal-pairs mode, sample names are reconstructed after beta conversion.
    if args.matrix_kind == "signal-pairs":
        sample_columns = sorted(
            {
                column[: -len(args.methylated_suffix)]
                for column in header
                if column.endswith(args.methylated_suffix)
            }
            & {
                column[: -len(args.unmethylated_suffix)]
                for column in header
                if column.endswith(args.unmethylated_suffix)
            }
        )
    else:
        sample_columns = [column for column in header if column != probe_column]

    sample_index = build_sample_index(
        sample_columns,
        sample_map=sample_map,
        excluded_columns=excluded,
    )

    print("Detected matched subjects:", file=sys.stderr)
    for brain_code in ("ba7", "ba10", "ba20"):
        shared = sorted(set(sample_index["blood"]) & set(sample_index[brain_code]))
        print(
            f"  Blood vs {TISSUE_LABELS[brain_code]}: {len(shared)} subjects "
            f"({', '.join(shared)})",
            file=sys.stderr,
        )

    if not any(
        len(set(sample_index["blood"]) & set(sample_index[brain_code])) >= args.min_pairs
        for brain_code in ("ba7", "ba10", "ba20")
    ):
        raise ValueError(
            "No brain region had enough matched blood–brain subjects. "
            "Inspect the sample names or provide --sample-map."
        )

    records: list[dict[str, Any]] = []
    read_kwargs = {
        "sep": separator,
        "chunksize": args.chunksize,
        "low_memory": False,
        "compression": "infer",
    }

    for chunk_number, chunk in enumerate(pd.read_csv(args.matrix, **read_kwargs), start=1):
        if args.matrix_kind == "signal-pairs":
            chunk = convert_signal_pairs_to_beta(
                chunk,
                probe_column=probe_column,
                methylated_suffix=args.methylated_suffix,
                unmethylated_suffix=args.unmethylated_suffix,
                offset=args.beta_offset,
            )

        chunk_records = process_chunk(
            chunk,
            probe_column=probe_column,
            sample_index=sample_index,
            min_pairs=args.min_pairs,
            min_abs_corr=args.min_abs_corr,
        )
        records.extend(chunk_records)
        print(
            f"Processed chunk {chunk_number}; cumulative observations: {len(records):,}",
            file=sys.stderr,
        )

    result = pd.DataFrame.from_records(records)
    if result.empty:
        raise RuntimeError("No finite CpG-level correlations were generated.")

    if args.annotation:
        annotation = read_annotation(args.annotation, args.annotation_sheet)
        result = result.merge(annotation, on="cpg", how="left")
    else:
        for column in ("chr", "pos", "gene", "context", "island_relation"):
            result[column] = ""

    ordered_columns = [
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
    result = result[ordered_columns].sort_values(
        ["brain_tissue", "cpg"],
        kind="stable",
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)

    print(
        f"Wrote {len(result):,} observations to {args.output}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
