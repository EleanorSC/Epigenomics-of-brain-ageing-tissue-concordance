# TRACE CpG: cross-tissue DNA methylation concordance in human brain and peripheral tissue samples

## Licensing

The source code in this repository is released under the MIT License.

Published supplementary datasets included in `data/raw/` remain the copyright of
their original authors and publishers and are distributed only for the purpose of
reproducing the published pipeline. Users should consult the original publications
for the applicable terms of use.

Any downloaded GEO datasets are not redistributed with this repository and must
be obtained from NCBI GEO.

# TRACE-CpG:
**TRACE-CpG**: **T**issue **R**eference-database for **A**nalysing **C**oncordant-**T**issue **E**pigenetics

There is currently no unified, curated resource that systematically collates published evidence on cross-tissue DNA methylation concordance at individual CpG sites. This project seeks to systemically review this research landscape and develop a harmonised and versioned database of published brain-peripheral DNA methylation measured from matched individuals, designed to be reproducibly updated as new studies become available.

The objective is to collate findings from published work that report correlations between DNAm measured in a peripheral tissue (e.g. blood, saliva) and brain tissue (e.g. as obtained from postmortem or biopsy samples). These genome-wide brain-peripheral DNA methylation concordance findings will then be aggregated into an openly available version-controlled database. 

## Repositry

This repository aims to curate a database of all published CpG-site correlations between peripheral tissues (blood, buccal) and brain regions (prefrontal cortex, superior temporal gyrus, entorhinal cortex, cerebellum), compiled from published cross-tissue methylation studies to date. At present the prototype database is built from published correlations from the following studies:

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

## Data Acquisition Scripts 
### Automated literature identification and screening

### `Extract_pubmed_screening.R`

Runs the preregistered PubMed search, retrieves citation metadata
and abstracts using the Entrez API, automatically prioritises
potentially eligible studies, and generates the screening CSV and
Excel workbooks used during study selection.

### `Extract_scopus_screening.R`

Runs the corresponding Scopus search using the Elsevier Scopus API and retrieves available citation metadata and abstracts. The script applies the same gate-based automated screening and prioritisation framework used for PubMed, while retaining Scopus-specific identifiers (EID and Scopus ID) and using DOI, EID, and title/year information for within-database deduplication.

Both pipelines:

- implement equivalent searches across PubMed and Scopus;
- retrieve and standardise bibliographic metadata and abstracts;
- remove duplicate records within each database;
- apply conservative automated exclusion rules according to the prespecified eligibility gates;
- retain records that cannot be confidently excluded for manual review;
- prioritise records containing evidence of brain tissue, peripheral tissue, matched/cross-tissue sampling, and genome-wide DNA methylation;
- use evidence of downloadable CpG-level results to further prioritise manual screening;
- generate screening CSV files and formatted Excel workbooks for independent review and adjudication; and
- generate summary statistics and interactive Sankey diagrams describing the screening process.

Automated screening is used for reproducible triage and prioritisation rather than as a substitute for investigator-led eligibility assessment. Records for which eligibility cannot be determined reliably from bibliographic metadata and abstracts are retained for manual screening.

### `Search_ewas_catalog.R`

Finally, to ensure relevant published findings are not missed from the literature searches, `Search_ewas_catalog.R`, Uses the EWAS Catalog as a supplementary source for identifying potentially eligible studies that may not be captured through the primary PubMed and Scopus searches.

The script:

- imports the downloaded EWAS Catalog study metadata;
- audits the tissue terminology used across Catalog records;
- preserves the original EWAS Catalog tissue labels while creating cleaned and harmonised tissue classifications;
- harmonises relevant tissue labels into broad `brain` and `peripheral` categories;
- identifies unique publications containing brain-tissue EWAS analyses for manual eligibility screening;
- identifies unique peripheral-tissue EWAS publications and their associated traits for additional checking of potentially relevant brain-related phenotypes; and
- generates tissue-frequency summaries that can be used to audit the coverage and harmonisation of EWAS Catalog tissue labels.

The EWAS Catalog search is used as a supplementary identification strategy rather than as a substitute for the preregistered PubMed and Scopus searches. Candidate publications identified through the Catalog are subsequently assessed against the same eligibility criteria and deduplicated against records identified through the bibliographic database searches.


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
