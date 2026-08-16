# CpG Brain-Peripheral Concordance Project

## Licensing

The source code in this repository is released under the MIT License.

Published supplementary datasets included in `data/raw/` remain the copyright of
their original authors and publishers and are distributed only for the purpose of
reproducing the published pipeline. Users should consult the original publications
for the applicable terms of use.

Any downloaded GEO datasets are not redistributed with this repository and must
be obtained from NCBI GEO.

# XTRACE-CpG:
**XTRACE-CpG**: **Cross**-**T**issue **R**eference for **A**nalysing **C**ross-**T**issue **E**pigenetics

There is currently no unified, curated resource that systematically collates published evidence on cross-tissue DNA methylation concordance at individual CpG sites. This project seeks to systemically review this research landscape and develop a harmonised and versioned database of published brain-peripheral DNA methylation measured from matched individuals, designed to be reproducibly updated as new studies become available.

The objective is to collate findings from published work that report correlations between DNAm measured in a peripheral tissue (e.g. blood, saliva) and brain tissue (e.g. as obtained from postmortem or biopsy samples). These genome-wide brain-peripheral DNA methylation concordance findings will then be aggregated into an openly available version-controlled database. 

## Repositry

This repository curates a database of 2,937 CpG-site correlations between peripheral tissues (blood, buccal) and brain regions (prefrontal cortex, superior temporal gyrus, entorhinal cortex, cerebellum), compiled from published cross-tissue methylation studies:

- Hannon E, et al. (2015). *Epigenetics*, **10**(11), 1024–1032. DOI: [10.1080/15592294.2015.1100786](https://doi.org/10.1080/15592294.2015.1100786)

- Edgar RD, et al. (2017). *Translational Psychiatry*, **7**, e1187. DOI: [10.1038/tp.2017.171](https://doi.org/10.1038/tp.2017.171)

- Braun PR, et al. (2019). *Translational Psychiatry*, **9**, 47. DOI: [10.1038/s41398-019-0376-y](https://doi.org/10.1038/s41398-019-0376-y)

- Sommerer Y, et al. (2022). *Clinical Epigenetics*, **14**, 118. DOI: [10.1186/s13148-022-01357-w](https://doi.org/10.1186/s13148-022-01357-w)

- Nishitani S, et al. (2023). *Translational Psychiatry*, **13**, 72. DOI: [10.1038/s41398-023-02370-0](https://doi.org/10.1038/s41398-023-02370-0)



<p align="left">
  <img src="Figures/glass_brain.png" width="800">
</p>

<p align="left">
<b>Figure 1.</b> Overview of the CpG Concordance Database. Published brain–peripheral DNA methylation concordance resources were harmonised into a unified searchable database and interactive web application.
</p>


## Overview


The database integrates all published human studies (through 2026) that provide genome-wide, CpG-level cross-tissue DNA methylation correlations between matched brain and peripheral tissues in a format suitable for systematic extraction. At the time of compilation, this comprised five primary resources: Hannon et al. (2015), Edgar et al. (2017), Braun et al. (2019), Sommerer et al. (2022), and Nishitani et al. (2023). Other blood–brain methylation studies exist but either examined limited candidate loci, did not publish genome-wide CpG-level concordance statistics, or did not provide supplementary data amenable to systematic extraction.

All 970 unique qualifying CpG sites (r/ρ > 0.6) are annotated with gene, chromosomal position (hg19), genomic context, and source study. Includes Python scripts to reproduce the full figure set: correlation distributions, tissue-pair comparisons, chromosomal and genomic-context breakdowns, top-gene rankings, and a genome-wide Manhattan-style plot of individual CpG sites — supporting research into reliable peripheral biomarkers of brain ageing and biological aging clocks

In the database, each row represents a distinct tissue-pair comparison (not just a CpG site). For some CpGs therefore, there will be multiple reported Pearson r-values depending on what correlations are reported; e.g. Hannon et al.'s blood-brain dataset tested correlations between blood methylation and methylation in three separate brain regions — entorhinal cortex (EC), superior temporal gyrus (STG), and prefrontal cortex (PFC) — for every probe in their dataset (Supplementary Table 7). So a single CpG like cg07249765 gets one Pearson r value per brain region it was compared against, since blood-EC correlation, blood-STG correlation, and blood-PFC correlation are three independent statistical results, each computed from the same set of paired samples but against a different brain region's methylation values.

<p align="left">
  <img src="Figures/database_overview.png" width="800">
</p>

<p align="left">
<b>Figure 2.</b> Overview of the CpG Concordance Database. Published brain–peripheral DNA methylation concordance resources were harmonised into a unified searchable database and interactive web application.
</p>

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

<p align="left">
  <img src="Figures/access_guide.png" width="800">
</p>

<p align="left">
<b>Figure 3.</b> User Guide for the searchable database and interactive web application.
</p>

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

# CpG Brain–Peripheral Concordance Database

This repository contains two Python scripts that extract published brain-peripheral DNA methylation concordance data and assemble the results into a searchable Excel workbook.

## Scripts

### `extract_source_data.py`

This is the **data-generation script**. It reads supplementary Excel files from the source studies, extracts CpG-level cross-tissue correlation estimates, annotates the CpGs, merges the datasets, and generates:

- `cpg_brain_peripheral_correlation_database.csv`
- `gene_level_summary.csv`
- `essex_high_corr.csv`
- `buccal_brain_high_corr.csv`

### `build_cpg_database.py`

This is the **presentation/export script**. It reads the two master CSV files and creates:

- `CpG_Brain_Peripheral_Correlation_Database.xlsx`

The workbook contains:

1. `README`
2. `CpG_Database`
3. `Gene_Summary`
4. `Additional_Curated_Examples`
5. `Database_Access_Guide`

## Current source data

The automated extraction currently incorporates two published supplementary datasets.

### Hannon et al. (2015)

**Study:** Hannon E, Lunnon K, Schalkwyk L, Mill J. *Interindividual methylomic variation across blood, cortex, and cerebellum: implications for epigenetic studies of neurological and neuropsychiatric phenotypes.* Epigenetics. 2015;10(11):1024–1032.

**Input workbook expected:** `essex_supp.xlsx`

**Sheets used:**

- `Supplementary Table 2`
  - CpG/probe annotation
  - chromosome
  - genomic position
  - UCSC RefGene symbol
  - UCSC RefGene group
  - relation to CpG island

- `Supplementary Table 6`
  - blood–brain Pearson correlations for the low-variance probe set

- `Supplementary Table 7`
  - blood–brain Pearson correlations for the highly variable probe set

**Brain regions extracted:**

- prefrontal cortex (PFC)
- entorhinal cortex (EC)
- superior temporal gyrus (STG)
- cerebellum (CER)

Each retained row represents one CpG × brain-region comparison.

### Sommerer et al. (2022)

**Study:** Sommerer Y, Ohlei O, Dobricic V, et al. *A correlation map of genome-wide DNA methylation patterns between paired human brain and buccal samples.* Clinical Epigenetics. 2022;14:118.

**Input workbook expected:** `sommerer_supp1.xlsx`

**Sheets used:**

- `S2`
  - CpG identifier
  - chromosome
  - genomic position
  - publication/reanalysis metadata

- `S6`
  - CpG-level Spearman correlations between buccal epithelium and prefrontal cortex
  - p values and q values
  - mQTL, CoRSIV, and IMAGE-CpG flags

Important: `S6` is already a study-defined subset of CpGs (the top q < 0.05 correlated probes). Removing the repository’s absolute-correlation threshold would retain all rows available in `S6`, but it would **not** recover the full genome-wide Sommerer correlation distribution.

## Current filtering rule

The current extraction script applies an absolute-correlation threshold:

```python
CORR_THRESHOLD = 0.6
```

Hannon rows are retained only when:

```python
abs(r) > 0.6
```

This is implemented by excluding rows meeting:

```python
if r is None or abs(r) <= CORR_THRESHOLD:
    continue
```

Sommerer rows are retained only when:

```python
abs(rho) > 0.6
```

This is implemented by excluding rows meeting:

```python
if not probe or rho is None or abs(rho) <= CORR_THRESHOLD:
    continue
```

Accordingly, the present database contains only observations with:

```text
|r| > 0.60 or |rho| > 0.60
```

Values exactly equal to 0.60 are excluded.

## Where threshold-dependent wording also appears

The threshold is not only implemented in `extract_source_data.py`; it is also hard-coded into labels and descriptive text in `build_cpg_database.py`, including:

- the module docstring
- the `Gene_Summary` column label `# Qualifying CpGs (r>0.6)`
- the workbook title
- the workbook README text
- references to “validated” or “high-correlation” CpGs

If the threshold is removed or changed, these descriptions must also be updated.

## Recommendation

The recommended design is to separate **data retention** from **default presentation**.

### Preferred approach

1. Retain every correlation value available in each source table.
2. Preserve source-specific significance statistics such as p values and q values where available.
3. Store the unfiltered observations in the master CSV.
4. Apply a user-selectable threshold only in the Shiny application or downstream analysis.
5. Use 0.60 as a default display filter, not as an irreversible extraction criterion.

This approach preserves the published evidence, supports sensitivity analyses, and allows the empirical distribution of correlation coefficients to be examined without rerunning extraction.

### Important limitation

For Sommerer et al., the currently used `S6` sheet is itself filtered to q < 0.05. A truly unfiltered Sommerer dataset requires a supplementary table or downloadable resource containing genome-wide CpG-level correlations. The script cannot reconstruct values that are absent from the input files.

## Suggested implementation

Replace:

```python
CORR_THRESHOLD = 0.6
```

with:

```python
CORR_THRESHOLD = None
```

Then replace the Hannon filter with:

```python
if r is None:
    continue

if CORR_THRESHOLD is not None and abs(r) <= CORR_THRESHOLD:
    continue
```

Replace the Sommerer filter with:

```python
if not probe or rho is None:
    continue

if CORR_THRESHOLD is not None and abs(rho) <= CORR_THRESHOLD:
    continue
```

This allows the same script to support either mode:

```python
CORR_THRESHOLD = None   # retain all available correlations
```

or:

```python
CORR_THRESHOLD = 0.6    # retain only |correlation| > 0.6
```

A more explicit alternative is:

```python
APPLY_CORRELATION_FILTER = False
CORR_THRESHOLD = 0.6
```

with:

```python
if APPLY_CORRELATION_FILTER and abs(r) <= CORR_THRESHOLD:
    continue
```

## Current annotation workflow

The extraction script builds a shared annotation lookup from Hannon Supplementary Table 2.

For each CpG, it retrieves or derives:

- CpG identifier
- chromosome
- genomic coordinate
- associated gene symbol(s)
- UCSC RefGene group
- CpG-island relation
- tissue-pair label
- correlation metric
- correlation coefficient
- source study
- source URL
- source database/tool

Genomic context is recoded as follows:

- `Promoter-associated`
  - TSS200
  - TSS1500
  - 5'UTR
  - 1stExon

- `Gene body/intragenic`
  - Body
  - 3'UTR

- `Unannotated/Intergenic`
  - no matching RefGene-group annotation

For Sommerer rows, chromosome and position are preferentially obtained from Sommerer `S2`; if unavailable, the script falls back to the Hannon/450K annotation lookup.

## Outputs

### Master CpG table

`cpg_brain_peripheral_correlation_database.csv`

One row per CpG × tissue-pair comparison, containing:

- `cpg`
- `chr`
- `pos`
- `gene`
- `context`
- `island_relation`
- `tissue_pair`
- `corr_type`
- `corr_value`
- `source`
- `source_url`
- `database_tool`

Correlation coefficients are rounded to three decimal places in the master CSV.

### Gene-level summary

`gene_level_summary.csv`

One row per gene, containing:

- chromosome
- number of unique retained CpGs
- maximum absolute correlation
- tissue pairs represented
- genomic contexts represented

### Excel workbook

`CpG_Brain_Peripheral_Correlation_Database.xlsx`

The primary `CpG_Database` worksheet excludes rows without a gene annotation. Those rows remain in the master CSV but are not displayed in the main gene-searchable Excel table.

The workbook also contains manually entered `Additional_Curated_Examples`. These are separate from the automated bulk extraction and should not be described as script-extracted observations.
