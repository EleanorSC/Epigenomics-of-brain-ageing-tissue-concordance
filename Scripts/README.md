# CpG Brain–Peripheral Concordance Database

This repository contains two Python scripts that extract published matched brain DNAm with peripheral DNAm (tissue concordance data) and assemble the results into a searchable Excel workbook.

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

## Reproducibility notes

- The scripts currently use absolute paths under `/home/user/workspace`.
- Input filenames and sheet names must match those expected by the code.
- The source workbooks must be downloaded separately.
- No raw methylation data are scraped from websites.
- The scripts read local supplementary Excel files.
- The scripts do not query the live BECon, IMAGE-CpG, AMAZE-CpG, Essex, or MADRC web interfaces programmatically.
- Correlation values are not manually changed during bulk extraction, but they are filtered and rounded.
- The manually curated examples in `build_cpg_database.py` are hard-coded and should be maintained separately from the reproducible bulk-extraction pipeline.

## Recommended next revisions

1. Remove irreversible correlation filtering from extraction.
2. Retain p and q values in the master schema where available.
3. Add a `source_table` column for provenance.
4. Preserve the unrounded correlation coefficient in the master CSV.
5. Add a derived `abs_corr` column for filtering and plotting.
6. Move manually curated examples into a separate CSV rather than hard-coding them.
7. Replace absolute workspace paths with command-line arguments or paths relative to the repository.
8. Update the Shiny app so that 0.60 is the default interactive filter rather than an extraction rule.
9. Generate threshold sensitivity summaries at 0.50, 0.60, 0.70, and 0.80.
10. Clearly distinguish:
   - all available source-table observations
   - study-significant observations
   - high-concordance observations selected by the application

## Running the pipeline

After placing the required supplementary workbooks in the configured workspace:

```bash
python3 extract_source_data.py
python3 build_cpg_database.py
```

The first script generates the CSV data products. The second script assembles the formatted Excel workbook.

## NOTES
Given the substantial tissue specificity of DNA methylation, strong cross-tissue correlations are relatively uncommon. A threshold of |r|≥0.60 was therefore chosen to enrich the database for CpGs demonstrating reproducible concordance across tissues, rather than ubiquitous tissue-specific methylation patterns. 

The objective of the database is to identify peripheral methylation markers that are likely to provide informative surrogates of brain methylation, rather than to catalogue all statistically significant correlations. Researchers can remove the threshold if they wish. The threshold is currently set in `extract_source_data.py`:

```python
CORR_THRESHOLD = 0.6
```

This is enforced twice within the script for specific papers:
For Hannon:

```python
if r is None or abs(r) <= CORR_THRESHOLD:
    continue
```

For Sommerer:
```python
if not probe or rho is None or abs(rho) <= CORR_THRESHOLD:
    continue
```
Such that the current script retains only observations with ∣r∣>0.60or ∣ρ∣>0.60
This will be clearly outlined in the ‘access guide’ on the interactive web tool itself.

To change this, set 
```python
CORR_THRESHOLD = None
```
and then:
```python
if r is None:
    continue

if CORR_THRESHOLD is not None and abs(r) <= CORR_THRESHOLD:
    continue
```

Nb. Sommerer S6 is already a study-filtered table of q-significant CpGs, so retaining all rows from S6 still means retaining all rows from a preselected subset, not all assayed CpGs. The extraction script reads local supplementary workbooks; it does not scrape the live web tools or raw methylation data

