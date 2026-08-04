#!/usr/bin/env python3
"""
00_run_cpg_concordance_pipeline.py

Run the five source-specific extraction scripts, harmonise their CSV outputs,
and generate the two master files consumed by build_cpg_database.py:

    cpg_brain_peripheral_correlation_database.csv
    gene_level_summary.csv

Source extractors
-----------------
01_extract_hannon.py
02_extract_edgar_becon.py
03_extract_braun_image_cpg.py
04_extract_sommerer_buccal_brain.py
05_extract_nishitani_amaze_cpg.py

The source-specific scripts retain the original study provenance and filtering
scope. This orchestration script does not impose an additional correlation
threshold unless --min-abs-corr is supplied.

Recommended repository structure
--------------------------------
Scripts/
    00_run_cpg_concordance_pipeline.py
    01_extract_hannon.py
    02_extract_edgar_becon.py
    03_extract_braun_image_cpg.py
    04_extract_sommerer_buccal_brain.py
    05_extract_nishitani_amaze_cpg.py
    build_cpg_database.py

data/raw/
    Hannon_SupplementaryTables.xlsx
    GSE95049_matrix_signal_BECon.txt.gz
    Braun_Supplementary_Table_3.xlsx
    Sommerer_13148_2022_1357_MOESM1_ESM.xlsx
    41398_2023_2370_MOESM3_ESM.xlsx

Database/intermediate/
Database/

Important methodological note
-----------------------------
The five inputs do not have identical scopes:

- Hannon: published subsets from Supplementary Tables 6 and 7.
- Edgar: all finite correlations recomputed from the supplied GSE95049 matrix.
- Braun: Bonferroni-significant Supplementary Table 3 results.
- Sommerer: q < 0.05 Supplementary Table S6 results.
- Nishitani: BH-significant AMAZE-CpG results; by default raw analyses only.

The master table therefore preserves `analysis_type`, `source_table` and
`source_scope`. It should not be described as five equivalent genome-wide
distributions.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable, Sequence

import pandas as pd


CORE_COLUMNS = [
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

OPTIONAL_COLUMNS = [
    # Hannon
    "refgene_group",
    "cpg_island",
    "variance_explained_percent",
    "mean_methylation_difference",
    "mean_difference_p_value",
    "published_subset",
    # Sommerer
    "buccal_mqtl",
    "blood_mqtl",
    "brain_mqtl",
    "corsiv",
    "braun_et_al",
    # Nishitani
    "cell_composition_adjusted",
    "snp_flag",
    "mqtl_flag",
    "cohort",
    "population",
]

MASTER_COLUMNS = CORE_COLUMNS + OPTIONAL_COLUMNS + [
    "record_id",
    "source_key",
]


SOURCE_KEYS = {
    "hannon": "Hannon_2015",
    "edgar": "Edgar_2017_BECon",
    "braun": "Braun_2019_IMAGE_CpG",
    "sommerer": "Sommerer_2022",
    "nishitani": "Nishitani_2023_AMAZE_CpG",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run all five CpG concordance extractors and build harmonised "
            "master CpG and gene-level CSV files."
        )
    )

    parser.add_argument(
        "--scripts-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory containing the five numbered extraction scripts.",
    )
    parser.add_argument(
        "--intermediate-dir",
        type=Path,
        default=Path("Database/intermediate"),
        help="Directory for source-specific extracted CSV files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("Database"),
        help="Directory for final master CSV files.",
    )
    parser.add_argument(
        "--python",
        default=sys.executable,
        help="Python interpreter used to run the source-specific scripts.",
    )

    # Hannon
    parser.add_argument("--hannon-workbook", type=Path, required=True)
    parser.add_argument(
        "--hannon-table",
        choices=("6", "7", "both"),
        default="both",
    )

    # Edgar
    parser.add_argument("--edgar-matrix", type=Path, required=True)
    parser.add_argument("--edgar-annotation", type=Path, default=None)
    parser.add_argument("--edgar-annotation-sheet", default=None)
    parser.add_argument(
        "--edgar-matrix-kind",
        choices=("beta", "signal-pairs"),
        default="beta",
    )
    parser.add_argument("--edgar-sample-map", type=Path, default=None)
    parser.add_argument("--edgar-probe-column", default=None)
    parser.add_argument("--edgar-separator", default=None)
    parser.add_argument("--edgar-min-pairs", type=int, default=3)

    # Braun
    parser.add_argument("--braun-workbook", type=Path, required=True)
    parser.add_argument("--braun-combined-sheet", default=None)
    parser.add_argument(
        "--braun-sheet",
        action="append",
        default=[],
        help="Optional Braun sheet; may be repeated.",
    )

    # Sommerer
    parser.add_argument("--sommerer-workbook", type=Path, required=True)
    parser.add_argument("--sommerer-annotation", type=Path, default=None)
    parser.add_argument("--sommerer-annotation-sheet", default=None)
    parser.add_argument("--sommerer-n-pairs", type=int, default=None)

    # Nishitani
    parser.add_argument("--nishitani-workbook", type=Path, required=True)
    parser.add_argument(
        "--nishitani-mode",
        choices=("raw", "adjusted", "both"),
        default="raw",
        help=(
            "Raw is recommended for the primary database. `both` retains "
            "cell-composition-adjusted rows as separate observations."
        ),
    )
    parser.add_argument("--nishitani-n-pairs", type=int, default=None)

    # Master-level behaviour
    parser.add_argument(
        "--min-abs-corr",
        type=float,
        default=None,
        help=(
            "Optional inclusive threshold applied after all source CSVs are "
            "combined. Default retains all extracted observations."
        ),
    )
    parser.add_argument(
        "--skip-extraction",
        action="store_true",
        help=(
            "Do not rerun source extractors; combine existing intermediate CSVs."
        ),
    )
    parser.add_argument(
        "--strict",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Fail on missing scripts, inputs, outputs or malformed source CSVs.",
    )
    return parser.parse_args()


def require_file(path: Path, label: str) -> None:
    if not path.exists() or not path.is_file():
        raise FileNotFoundError(f"{label} was not found: {path}")


def run_command(command: Sequence[str], label: str) -> None:
    print(f"\n[{label}]", file=sys.stderr)
    print(" ".join(str(part) for part in command), file=sys.stderr)
    completed = subprocess.run(command, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            f"{label} failed with exit code {completed.returncode}."
        )


def add_optional_argument(
    command: list[str],
    flag: str,
    value: Any,
) -> None:
    if value is not None:
        command.extend([flag, str(value)])


def extraction_paths(intermediate_dir: Path) -> dict[str, Path]:
    return {
        "hannon": intermediate_dir / "hannon_blood_brain_correlations.csv",
        "edgar": intermediate_dir / "edgar_becon_cpg_correlations.csv",
        "braun": intermediate_dir / "braun_image_cpg_significant_correlations.csv",
        "sommerer": intermediate_dir / "sommerer_buccal_brain_correlations.csv",
        "nishitani": intermediate_dir / "nishitani_amaze_cpg_correlations.csv",
    }


def run_extractors(args: argparse.Namespace, outputs: dict[str, Path]) -> None:
    scripts = {
        "hannon": args.scripts_dir / "01_extract_hannon.py",
        "edgar": args.scripts_dir / "02_extract_edgar_becon.py",
        "braun": args.scripts_dir / "03_extract_braun_image_cpg.py",
        "sommerer": args.scripts_dir / "04_extract_sommerer_buccal_brain.py",
        "nishitani": args.scripts_dir / "05_extract_nishitani_amaze_cpg.py",
    }

    for key, script in scripts.items():
        require_file(script, f"{key} extraction script")

    input_files = {
        "Hannon workbook": args.hannon_workbook,
        "Edgar matrix": args.edgar_matrix,
        "Braun workbook": args.braun_workbook,
        "Sommerer workbook": args.sommerer_workbook,
        "Nishitani workbook": args.nishitani_workbook,
    }
    for label, path in input_files.items():
        require_file(path, label)

    if args.edgar_annotation:
        require_file(args.edgar_annotation, "Edgar annotation")
    if args.edgar_sample_map:
        require_file(args.edgar_sample_map, "Edgar sample map")
    if args.sommerer_annotation:
        require_file(args.sommerer_annotation, "Sommerer annotation")

    # 1. Hannon
    command = [
        args.python,
        str(scripts["hannon"]),
        "--input",
        str(args.hannon_workbook),
        "--output",
        str(outputs["hannon"]),
        "--table",
        args.hannon_table,
    ]
    run_command(command, "Hannon extraction")

    # 2. Edgar
    command = [
        args.python,
        str(scripts["edgar"]),
        "--matrix",
        str(args.edgar_matrix),
        "--output",
        str(outputs["edgar"]),
        "--matrix-kind",
        args.edgar_matrix_kind,
        "--min-pairs",
        str(args.edgar_min_pairs),
    ]
    add_optional_argument(command, "--annotation", args.edgar_annotation)
    add_optional_argument(
        command,
        "--annotation-sheet",
        args.edgar_annotation_sheet,
    )
    add_optional_argument(command, "--sample-map", args.edgar_sample_map)
    add_optional_argument(command, "--probe-column", args.edgar_probe_column)
    add_optional_argument(command, "--separator", args.edgar_separator)
    run_command(command, "Edgar/BECon extraction")

    # 3. Braun
    command = [
        args.python,
        str(scripts["braun"]),
        "--input",
        str(args.braun_workbook),
        "--output",
        str(outputs["braun"]),
    ]
    if args.braun_combined_sheet:
        command.extend(["--combined-sheet", args.braun_combined_sheet])
    for sheet in args.braun_sheet:
        command.extend(["--sheet", sheet])
    run_command(command, "Braun/IMAGE-CpG extraction")

    # 4. Sommerer
    command = [
        args.python,
        str(scripts["sommerer"]),
        "--input",
        str(args.sommerer_workbook),
        "--output",
        str(outputs["sommerer"]),
    ]
    add_optional_argument(command, "--annotation", args.sommerer_annotation)
    add_optional_argument(
        command,
        "--annotation-sheet",
        args.sommerer_annotation_sheet,
    )
    add_optional_argument(command, "--n-pairs", args.sommerer_n_pairs)
    run_command(command, "Sommerer extraction")

    # 5. Nishitani
    command = [
        args.python,
        str(scripts["nishitani"]),
        "--input",
        str(args.nishitani_workbook),
        "--output",
        str(outputs["nishitani"]),
    ]

    # Supports both versions of the Nishitani script:
    # - the broader script with --mode;
    # - the AMAZE raw-only script without --mode.
    script_text = scripts["nishitani"].read_text(encoding="utf-8")
    if "--mode" in script_text:
        command.extend(["--mode", args.nishitani_mode])
    elif args.nishitani_mode != "raw":
        raise ValueError(
            "The installed Nishitani extractor supports raw AMAZE-CpG rows only, "
            f"but --nishitani-mode {args.nishitani_mode!r} was requested."
        )

    add_optional_argument(command, "--n-pairs", args.nishitani_n_pairs)
    run_command(command, "Nishitani/AMAZE-CpG extraction")


def clean_string_series(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .replace({"nan": "", "None": "", "<NA>": ""})
        .str.strip()
    )


def standardise_frame(
    frame: pd.DataFrame,
    *,
    source_key: str,
    strict: bool,
) -> pd.DataFrame:
    frame = frame.copy()
    frame.columns = [str(column).strip() for column in frame.columns]

    required = {
        "cpg",
        "peripheral_tissue",
        "brain_tissue",
        "corr_type",
        "corr_value",
        "source",
    }
    missing = required - set(frame.columns)
    if missing:
        message = (
            f"{source_key}: missing required columns {sorted(missing)}; "
            f"observed columns: {list(frame.columns)}"
        )
        if strict:
            raise ValueError(message)
        print(f"WARNING: {message}", file=sys.stderr)

    for column in MASTER_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""

    string_columns = [
        column
        for column in MASTER_COLUMNS
        if column not in {"corr_value", "abs_corr", "p_value", "q_value", "n_pairs"}
    ]
    for column in string_columns:
        frame[column] = clean_string_series(frame[column])

    frame["cpg"] = frame["cpg"].str.replace(r"\*$", "", regex=True)
    frame = frame.loc[frame["cpg"].str.match(r"^cg\d+$", na=False)].copy()

    frame["corr_value"] = pd.to_numeric(frame["corr_value"], errors="coerce")
    frame = frame.loc[frame["corr_value"].notna()].copy()
    frame = frame.loc[frame["corr_value"].between(-1, 1)].copy()

    calculated_abs = frame["corr_value"].abs()
    supplied_abs = pd.to_numeric(frame["abs_corr"], errors="coerce")
    frame["abs_corr"] = supplied_abs.fillna(calculated_abs)

    frame["p_value"] = pd.to_numeric(frame["p_value"], errors="coerce")
    frame["q_value"] = pd.to_numeric(frame["q_value"], errors="coerce")
    frame["n_pairs"] = pd.to_numeric(frame["n_pairs"], errors="coerce").astype("Int64")

    blank_pair = frame["tissue_pair"].eq("")
    frame.loc[blank_pair, "tissue_pair"] = (
        frame.loc[blank_pair, "peripheral_tissue"]
        + " vs "
        + frame.loc[blank_pair, "brain_tissue"]
    )

    frame["source_key"] = source_key
    frame["record_id"] = (
        source_key
        + "|"
        + frame["analysis_type"]
        + "|"
        + frame["cpg"]
        + "|"
        + frame["peripheral_tissue"]
        + "|"
        + frame["brain_tissue"]
        + "|"
        + frame["corr_type"]
    )

    return frame[MASTER_COLUMNS]


def load_and_combine(
    outputs: dict[str, Path],
    *,
    strict: bool,
    min_abs_corr: float | None,
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []

    for key, path in outputs.items():
        if not path.exists():
            if strict:
                raise FileNotFoundError(f"Expected source CSV was not found: {path}")
            print(f"WARNING: skipping missing source CSV: {path}", file=sys.stderr)
            continue

        frame = pd.read_csv(path, low_memory=False)
        standardised = standardise_frame(
            frame,
            source_key=SOURCE_KEYS[key],
            strict=strict,
        )
        print(
            f"Loaded {len(standardised):,} valid observations from {path.name}",
            file=sys.stderr,
        )
        frames.append(standardised)

    if not frames:
        raise RuntimeError("No source CSVs were available for combination.")

    combined = pd.concat(frames, ignore_index=True, sort=False)

    if min_abs_corr is not None:
        combined = combined.loc[combined["abs_corr"] >= min_abs_corr].copy()

    # Only remove exact duplicate records. Distinct studies, analyses, tissue
    # pairs and source tables remain separate evidence rows.
    duplicate_key = [
        "source_key",
        "analysis_type",
        "source_table",
        "cpg",
        "peripheral_tissue",
        "brain_tissue",
        "corr_type",
        "corr_value",
    ]
    before = len(combined)
    combined = combined.drop_duplicates(subset=duplicate_key, keep="first")
    removed = before - len(combined)
    if removed:
        print(f"Removed {removed:,} exact duplicate observations", file=sys.stderr)

    return combined.sort_values(
        [
            "source_key",
            "analysis_type",
            "peripheral_tissue",
            "brain_tissue",
            "cpg",
        ],
        kind="stable",
    ).reset_index(drop=True)


def split_genes(value: Any) -> list[str]:
    text = "" if pd.isna(value) else str(value).strip()
    if not text:
        return []

    output: list[str] = []
    for token in re_split_genes(text):
        gene = token.strip()
        if gene and gene not in output:
            output.append(gene)
    return output


def re_split_genes(value: str) -> list[str]:
    # Most source annotations use semicolons. Commas are also accepted.
    import re
    return re.split(r"[;,]", value)


def build_gene_summary(master: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, Any]] = []

    for row in master.itertuples(index=False):
        genes = split_genes(row.gene)
        for gene in genes:
            records.append(
                {
                    "gene": gene,
                    "chr": row.chr,
                    "cpg": row.cpg,
                    "abs_corr": row.abs_corr,
                    "corr_value": row.corr_value,
                    "tissue_pair": row.tissue_pair,
                    "peripheral_tissue": row.peripheral_tissue,
                    "brain_tissue": row.brain_tissue,
                    "context": row.context,
                    "source_key": row.source_key,
                    "source": row.source,
                    "analysis_type": row.analysis_type,
                }
            )

    columns = [
        "gene",
        "chr",
        "n_unique_cpgs",
        "n_observations",
        "max_abs_corr",
        "mean_abs_corr",
        "tissue_pairs",
        "peripheral_tissues",
        "brain_tissues",
        "sources",
        "analysis_types",
        "contexts",
    ]

    if not records:
        return pd.DataFrame(columns=columns)

    long = pd.DataFrame.from_records(records)

    def join_unique(series: pd.Series) -> str:
        values = sorted(
            {
                str(value).strip()
                for value in series
                if pd.notna(value) and str(value).strip()
            }
        )
        return "; ".join(values)

    summary = (
        long.groupby("gene", dropna=False)
        .agg(
            chr=("chr", join_unique),
            n_unique_cpgs=("cpg", "nunique"),
            n_observations=("cpg", "size"),
            max_abs_corr=("abs_corr", "max"),
            mean_abs_corr=("abs_corr", "mean"),
            tissue_pairs=("tissue_pair", join_unique),
            peripheral_tissues=("peripheral_tissue", join_unique),
            brain_tissues=("brain_tissue", join_unique),
            sources=("source", join_unique),
            analysis_types=("analysis_type", join_unique),
            contexts=("context", join_unique),
        )
        .reset_index()
    )

    summary["max_abs_corr"] = summary["max_abs_corr"].round(6)
    summary["mean_abs_corr"] = summary["mean_abs_corr"].round(6)

    return summary.sort_values(
        ["n_unique_cpgs", "max_abs_corr", "gene"],
        ascending=[False, False, True],
        kind="stable",
    ).reset_index(drop=True)[columns]


def build_source_summary(master: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        master.groupby(
            [
                "source_key",
                "source",
                "analysis_type",
                "peripheral_tissue",
                "brain_tissue",
                "source_scope",
            ],
            dropna=False,
        )
        .agg(
            n_observations=("cpg", "size"),
            n_unique_cpgs=("cpg", "nunique"),
            min_corr=("corr_value", "min"),
            max_corr=("corr_value", "max"),
            median_abs_corr=("abs_corr", "median"),
        )
        .reset_index()
    )
    return grouped.sort_values(
        ["source_key", "analysis_type", "peripheral_tissue", "brain_tissue"],
        kind="stable",
    )


def write_manifest(
    path: Path,
    *,
    args: argparse.Namespace,
    outputs: dict[str, Path],
    master: pd.DataFrame,
) -> None:
    manifest = {
        "master_rows": int(len(master)),
        "unique_cpgs": int(master["cpg"].nunique()),
        "sources": {
            key: {
                "csv": str(csv_path),
                "rows": int((master["source_key"] == SOURCE_KEYS[key]).sum()),
            }
            for key, csv_path in outputs.items()
        },
        "settings": {
            "hannon_table": args.hannon_table,
            "nishitani_mode": args.nishitani_mode,
            "min_abs_corr": args.min_abs_corr,
            "edgar_matrix_kind": args.edgar_matrix_kind,
            "edgar_min_pairs": args.edgar_min_pairs,
        },
    }
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> int:
    args = parse_args()

    if args.min_abs_corr is not None and not 0 <= args.min_abs_corr <= 1:
        raise ValueError("--min-abs-corr must be between 0 and 1.")

    args.intermediate_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    outputs = extraction_paths(args.intermediate_dir)

    if not args.skip_extraction:
        run_extractors(args, outputs)

    master = load_and_combine(
        outputs,
        strict=args.strict,
        min_abs_corr=args.min_abs_corr,
    )
    gene_summary = build_gene_summary(master)
    source_summary = build_source_summary(master)

    master_path = args.output_dir / "cpg_brain_peripheral_correlation_database.csv"
    gene_path = args.output_dir / "gene_level_summary.csv"
    source_path = args.output_dir / "source_level_summary.csv"
    manifest_path = args.output_dir / "pipeline_manifest.json"

    master.to_csv(master_path, index=False)
    gene_summary.to_csv(gene_path, index=False)
    source_summary.to_csv(source_path, index=False)
    write_manifest(
        manifest_path,
        args=args,
        outputs=outputs,
        master=master,
    )

    print("\nPipeline complete", file=sys.stderr)
    print(f"Master observations: {len(master):,}", file=sys.stderr)
    print(f"Unique CpGs: {master['cpg'].nunique():,}", file=sys.stderr)
    print(f"Genes: {len(gene_summary):,}", file=sys.stderr)
    print(f"Master CSV: {master_path}", file=sys.stderr)
    print(f"Gene summary: {gene_path}", file=sys.stderr)
    print(f"Source summary: {source_path}", file=sys.stderr)
    print(f"Manifest: {manifest_path}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
