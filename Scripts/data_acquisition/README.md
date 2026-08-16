# Data Acquisition Scripts 
## Automated literature identification and screening

### `Search_ewas_catalog.R`

Uses the EWAS Catalog as a supplementary source for identifying potentially eligible studies that may not be captured through the primary PubMed and Scopus searches.

The script:

- imports the downloaded EWAS Catalog study metadata;
- audits the tissue terminology used across Catalog records;
- preserves the original EWAS Catalog tissue labels while creating cleaned and harmonised tissue classifications;
- harmonises relevant tissue labels into broad `brain` and `peripheral` categories;
- identifies unique publications containing brain-tissue EWAS analyses for manual eligibility screening;
- identifies unique peripheral-tissue EWAS publications and their associated traits for additional checking of potentially relevant brain-related phenotypes; and
- generates tissue-frequency summaries that can be used to audit the coverage and harmonisation of EWAS Catalog tissue labels.

The EWAS Catalog search is used as a supplementary identification strategy rather than as a substitute for the preregistered PubMed and Scopus searches. Candidate publications identified through the Catalog are subsequently assessed against the same eligibility criteria and deduplicated against records identified through the bibliographic database searches.

### Extract_pubmed_screening.R

Runs the preregistered PubMed search, retrieves citation metadata
and abstracts using the Entrez API, automatically prioritises
potentially eligible studies, and generates the screening CSV and
Excel workbooks used during study selection.

### Extract_scopus_screening.R

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

