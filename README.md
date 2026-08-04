# Epigenomics-of-brain-ageing-tissue-concordance

This repository curates a database of 2,937 CpG-site correlations between peripheral tissues (blood, buccal) and brain regions (prefrontal cortex, superior temporal gyrus, entorhinal cortex, cerebellum), compiled from published cross-tissue methylation studies:

- Hannon E, et al. (2015). *Epigenetics*, **10**(11), 1024–1032. DOI: [10.1080/15592294.2015.1100786](https://doi.org/10.1080/15592294.2015.1100786)

- Edgar RD, et al. (2017). *Translational Psychiatry*, **7**, e1187. DOI: [10.1038/tp.2017.171](https://doi.org/10.1038/tp.2017.171)

- Braun PR, et al. (2019). *Translational Psychiatry*, **9**, 47. DOI: [10.1038/s41398-019-0376-y](https://doi.org/10.1038/s41398-019-0376-y)

- Sommerer Y, et al. (2022). *Clinical Epigenetics*, **14**, 118. DOI: [10.1186/s13148-022-01357-w](https://doi.org/10.1186/s13148-022-01357-w)

- Nishitani S, et al. (2023). *Translational Psychiatry*, **13**, 72. DOI: [10.1038/s41398-023-02370-0](https://doi.org/10.1038/s41398-023-02370-0)


<img width="2440" height="1620" alt="image" src="https://github.com/user-attachments/assets/5fa73741-29bf-4576-9abb-1bedd274e973" />



shiny_app/
├── app.R
├── data/
│   └── cpg_brain_peripheral_correlation_database.csv
└── www/
    ├── brain.css
    └── brain.js

The database integrates all published human studies (through 2026) that provide genome-wide, CpG-level cross-tissue DNA methylation correlations between matched brain and peripheral tissues in a format suitable for systematic extraction. At the time of compilation, this comprised five primary resources: Hannon et al. (2015), Edgar et al. (2017), Braun et al. (2019), Sommerer et al. (2022), and Nishitani et al. (2023). Other blood–brain methylation studies exist but either examined limited candidate loci, did not publish genome-wide CpG-level concordance statistics, or did not provide supplementary data amenable to systematic extraction.

All 970 unique qualifying CpG sites (r/ρ > 0.6) are annotated with gene, chromosomal position (hg19), genomic context, and source study. Includes Python scripts to reproduce the full figure set: correlation distributions, tissue-pair comparisons, chromosomal and genomic-context breakdowns, top-gene rankings, and a genome-wide Manhattan-style plot of individual CpG sites — supporting research into reliable peripheral biomarkers of brain ageing and biological aging clocks

In the database, each row represents a distinct tissue-pair comparison (not just a CpG site). For some CpGs therefore, there will be multiple reported Pearson r-values depending on what correlations are reported; e.g. Hannon et al.'s blood-brain dataset tested correlations between blood methylation and methylation in three separate brain regions — entorhinal cortex (EC), superior temporal gyrus (STG), and prefrontal cortex (PFC) — for every probe in their dataset (Supplementary Table 7). So a single CpG like cg07249765 gets one Pearson r value per brain region it was compared against, since blood-EC correlation, blood-STG correlation, and blood-PFC correlation are three independent statistical results, each computed from the same set of paired samples but against a different brain region's methylation values.

This is why the tissue_pair column exists as a distinct field from CpG — the true unique key in the database is the combination of (CpG ID, tissue pair), not CpG ID alone. 

Using cg07249765 as an example, this site happens to correlate strongly and consistently across all three brain regions (0.994–0.996), which itself is a meaningful finding — some CpGs are region-specific proxies while others, like this one, are robust blood-based proxies for cortical methylation state regardless of which brain region is being modelled.

## Database
### Database of CpG–Brain Peripheral Correlations

This repository contains two scripts with distinct purposes:

- **`database_extract_source_data.py`** — extracts, harmonises, and processes the original published supplementary datasets into standardised CSV files.
- **`build_cpg_database.py`** — builds the final Excel workbook from those processed CSVs, applying formatting and adding manually curated content.

---

## 1. Database Data Extraction

### `database_extract_source_data.py`

This script extracts and harmonises data from the original published supplementary tables and generates the source CSV files used by the database.

### Workflow

### 1. Load probe annotation

Loads the shared probe annotation manifest from:

- **Hannon et al. (2015), Supplementary Table 2**

This contains annotation for **22,458 CpG probes**, including:

- chromosome
- genomic position
- associated gene(s)
- genomic-context group
- CpG-island relation

These annotations are reused for both the blood–brain and buccal–brain datasets.

---

### 2. Extract blood–brain correlations

Reads:

- Hannon et al. (2015)
  - Supplementary Table 6
  - Supplementary Table 7

The script:

- extracts all CpG correlations
- filters to **|r| > 0.60**
- joins chromosome, position, gene and genomic annotation from the shared probe manifest

Output:

- **7,252 qualifying CpG–brain correlation rows**

---

### 3. Extract buccal–brain correlations

Reads:

- Sommerer et al. (2022)
  - Supplementary Table S6

The script:

- extracts CpG correlations
- filters to **|ρ| > 0.60**
- retrieves chromosome and genomic position primarily from:
  - Sommerer Supplementary Table S2 (re-annotation)
- falls back to the Hannon probe manifest for probes not present in S2

Output:

- **1,875 qualifying CpG–brain correlation rows**

---

### 4. Build the master database

The blood–brain and buccal–brain datasets are merged into a single master dataset.

Additional processing includes:

- genomic-context classification
- tissue-pair labels
- citation metadata
- numeric rounding
- consistent variable naming

Output:

- **9,127 CpG × tissue-pair observations**

Saved as:

```
cpg_brain_peripheral_correlation_database.csv
```

---

### 5. Generate gene-level summary

The master database is aggregated to gene level by:

- removing duplicated gene symbols
- splitting multi-gene annotations
- collapsing repeated transcript annotations

Output:

- **815 unique genes**

Saved as:

```
gene_level_summary.csv
```

---

# 2. Database Workbook Construction

## `build_cpg_database.py`

This standalone script reproduces the complete Excel workbook directly from the CSV files generated by `database_extract_source_data.py`.

### Workflow

### 1. Build the CpG database sheet

Reads:

```
cpg_brain_peripheral_correlation_database.csv
```

The script:

- imports all **9,127** CpG × tissue-pair observations
- filters to the **2,937** rows with gene annotations
- removes duplicated gene symbols (e.g. `SDK1;SDK1` → `SDK1`)
- prefixes chromosome identifiers with `chr`

These data populate the **CpG Database** worksheet.

---

### 2. Build the Gene Summary sheet

Reads:

```
gene_level_summary.csv
```

containing:

- **815 unique genes**

This file is imported directly into the **Gene Summary** worksheet without further modification.

---

### 3. Add curated literature examples

The script inserts:

- **12 manually verified example studies**
- **5-resource access guide**

These entries are hard-coded because they were curated manually from the literature rather than generated through automated extraction.

---

### 4. Generate the README worksheet

Automatically creates a README worksheet containing:

- database overview
- data provenance
- construction workflow
- usage guidance
- citation information

---

### 5. Apply workbook formatting

The script applies consistent formatting across all worksheets, including:

- Nexus-style teal header formatting
- frozen header rows
- autofilters
- column widths
- numeric formatting
- workbook styling

---

# Outputs

Running both scripts produces the complete database.

```
database_extract_source_data.py
    │
    ├── cpg_brain_peripheral_correlation_database.csv
    └── gene_level_summary.csv
            │
            ▼
build_cpg_database.py
            │
            ▼
CpG–Brain Peripheral Correlation Database.xlsx
```
