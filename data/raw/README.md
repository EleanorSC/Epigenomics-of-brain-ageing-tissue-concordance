# data/raw

This directory contains the original source files used to construct the **CpG Brain–Peripheral Correlation Database**.

All files in this directory are **read-only source data**. They should never be modified manually.

The extraction scripts in `Scripts/` parse these files and produce harmonised intermediate CSV files in `Database/intermediate/`, which are subsequently merged into the master database.

---

## Directory contents

https://epigenetics.essex.ac.uk/bloodbrain/

| Original File Name | Source | Citation | doi | Used by |
| `SupplementaryTables.xlsx` | https://epigenetics.essex.ac.uk/bloodbrain/ | Hannon et al. (2015) *Epigenetics* |10.1080/15592294.2015.1100786 |`01_extract_hannon.py` |


| Original File Name                                                       | Source                                                                                                                     | Citation                                                 | DOI                           | Used by                               |
| ------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------- | ----------------------------- | ------------------------------------- |
| `SupplementaryTables.xlsx`                                               | [https://epigenetics.essex.ac.uk/bloodbrain/](https://epigenetics.essex.ac.uk/bloodbrain/)                                 | Hannon E, *et al.* (2015). *Epigenetics*                 | 10.1080/15592294.2015.1100786 | `01_extract_hannon.py`                |
| `GSE95049_series_matrix.txt.gz` *(or whichever GEO file you downloaded)* | [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95049](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95049) | Edgar RD, *et al.* (2017). *Translational Psychiatry*    | 10.1038/tp.2017.171           | `02_extract_edgar_becon.py`           |
| `Supplementary_Table_3.xlsx`                                             | Supplementary material                                                                                                     | Braun PR, *et al.* (2019). *Translational Psychiatry*    | 10.1038/s41398-019-0376-y     | `03_extract_braun_image_cpg.py`       |
| `13148_2022_1357_MOESM1_ESM.xlsx`                                        | Supplementary material                                                                                                     | Sommerer Y, *et al.* (2022). *Clinical Epigenetics*      | 10.1186/s13148-022-01357-w    | `04_extract_sommerer_buccal_brain.py` |
| `41398_2023_2370_MOESM3_ESM.xlsx`                                        | Supplementary material                                                                                                     | Nishitani S, *et al.* (2023). *Translational Psychiatry* | 10.1038/s41398-023-02370-0    | `05_extract_nishitani_amaze_cpg.py`   |



|------|--------|---------|
| `Hannon_SupplementaryTables.xlsx` | Hannon et al. (2015) *Epigenetics* | `01_extract_hannon.py` |
| `GSE95049_processed_matrix.*` | Edgar et al. (2017) BECon / GSE95049 | `02_extract_edgar_becon.py` |
| `Braun_Supplementary_Table_3.xlsx` | Braun et al. (2019) *Translational Psychiatry* | `03_extract_braun_image_cpg.py` |
| `Sommerer_13148_2022_1357_MOESM1_ESM.xlsx` | Sommerer et al. (2022) *Clinical Epigenetics* | `04_extract_sommerer_buccal_brain.py` |
| `41398_2023_2370_MOESM3_ESM.xlsx` | Nishitani et al. (2023) *Translational Psychiatry* | `05_extract_nishitani_amaze_cpg.py` |

---

## Data sources

### Hannon et al. (2015)

**Reference**

> Hannon E, Lunnon K, Schalkwyk L, Mill J. *Interindividual methylomic variation across blood, cortex, and cerebellum: implications for epigenetic studies of neurological and neuropsychiatric phenotypes.* Epigenetics. 2015;10(11):1024–1032.

**Required file**

```
Hannon_SupplementaryTables.xlsx
```

**Used tables**

- Supplementary Table 2 (probe annotation)
- Supplementary Table 6 (published blood–brain correlations)
- Supplementary Table 7 (additional published concordance subset)

---

### Edgar et al. (2017) (BECon)

**Reference**

> Edgar RD et al. *BECon: a tool for interpreting DNA methylation findings from blood in the context of brain.* Translational Psychiatry. 2017;7:e1187.

**Required file**

```
GSE95049_processed_matrix.*
```

This is the processed methylation matrix downloaded from GEO (accession **GSE95049**).

Unlike the other studies, CpG correlations are **recomputed** from the processed matrix rather than extracted directly from supplementary tables.

Optional Illumina probe annotation files may also be placed in this directory.

---

### Braun et al. (2019)

**Reference**

> Braun PR et al. *Genome-wide DNA methylation comparison between live human brain and peripheral tissues within individuals.* Translational Psychiatry. 2019;9:47.

**Required file**

```
Braun_Supplementary_Table_3.xlsx
```

**Used table**

- Supplementary Table 3

The extraction script imports the published Bonferroni-significant CpGs.

---

### Sommerer et al. (2022)

**Reference**

> Sommerer Y et al. *A correlation map of genome-wide DNA methylation patterns between paired human brain and buccal samples.* Clinical Epigenetics. 2022;14:118.

**Required file**

```
Sommerer_13148_2022_1357_MOESM1_ESM.xlsx
```

**Used tables**

- Supplementary Table S2
- Supplementary Table S6

The extraction script imports all published q < 0.05 CpGs.

---

### Nishitani et al. (2023)

**Reference**

> Nishitani S et al. *Cross-tissue correlations of genome-wide DNA methylation in Japanese live human brain and blood, saliva, and buccal epithelial tissues.* Translational Psychiatry. 2023;13:72.

**Required file**

```
41398_2023_2370_MOESM3_ESM.xlsx
```

The script extracts only the **independent AMAZE-CpG cohort**.

By default, the database uses the **raw** brain–blood, brain–saliva and brain–buccal analyses.

The supplied R scripts for GSE59685 and GSE95049 are **not** used, as they reproduce analyses from previously published cohorts rather than generating independent observations.

---

## Workflow

```
data/raw/
        │
        ▼
01_extract_hannon.py
02_extract_edgar_becon.py
03_extract_braun_image_cpg.py
04_extract_sommerer_buccal_brain.py
05_extract_nishitani_amaze_cpg.py
        │
        ▼
Database/intermediate/
        │
        ▼
00_run_cpg_concordance_pipeline.py
        │
        ▼
Database/
    ├── cpg_brain_peripheral_correlation_database.csv
    ├── gene_level_summary.csv
    ├── source_level_summary.csv
    └── pipeline_manifest.json
        │
        ▼
build_cpg_database.py
        │
        ▼
CpG_Brain_Peripheral_Correlation_Database.xlsx
```

---

## Notes

- Original source files should **never be edited**.
- Any updates or corrections should be implemented in the extraction scripts rather than by modifying the downloaded supplementary material.
- Each extraction script records the original publication, supplementary table, and analysis type to preserve full provenance throughout the pipeline.
- The five studies differ in scope. Some provide published significant subsets, whereas the BECon dataset is reconstructed from public processed data. These differences are retained in the database metadata and should be considered when interpreting coverage.
