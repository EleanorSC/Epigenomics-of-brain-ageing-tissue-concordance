# Data Acquisition Scripts 
## Automated literature identification and screening

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

To ensure relevant published findings are not missed from the literature searches, `Search_ewas_catalog.R`, Uses the EWAS Catalog as a supplementary source for identifying potentially eligible studies that may not be captured through the primary PubMed and Scopus searches.

The script:

- imports the downloaded EWAS Catalog study metadata;
- audits the tissue terminology used across Catalog records;
- preserves the original EWAS Catalog tissue labels while creating cleaned and harmonised tissue classifications;
- harmonises relevant tissue labels into broad `brain` and `peripheral` categories;
- identifies unique publications containing brain-tissue EWAS analyses for manual eligibility screening;
- identifies unique peripheral-tissue EWAS publications and their associated traits for additional checking of potentially relevant brain-related phenotypes; and
- generates tissue-frequency summaries that can be used to audit the coverage and harmonisation of EWAS Catalog tissue labels.

### `Search_GEO_resource.R`

Finally, the `Search_GEO_resource.R` script uses the NCBI Gene Expression Omnibus (GEO) as a final supplementary study-identification resource to identify potentially relevant human DNA methylation datasets containing both brain and peripheral tissue samples, with particular attention to GEO Series not linked to an indexed PubMed publication.

The script:

- searches GEO DataSets programmatically for human DNA methylation Series containing evidence of both brain and peripheral tissue samples;
- retrieves GEO Series-level metadata using the NCBI Entrez API;
- cross-references identified GEO Series against PubMed to determine whether an indexed publication is linked to the dataset;
- inspects sample-level GEO metadata to confirm the presence of brain and peripheral samples and identify terminology indicative of matched, paired, or common-donor sampling;
- flags potentially relevant GEO Series without a linked PubMed/Scopus publication record;
- prioritises candidate Series containing stronger evidence of matched brain–peripheral sampling; and
- exports candidate datasets and associated metadata for manual eligibility assessment by Independent reviewers.


#### Of Note
Both The EWAS Catalog search and GEO search are used as a supplementary identifications strategies. If new studies are found here, the absence of a linked publication record indexed in either PubMed or Scopus is not interpreted as evidence that a dataset is unpublished; such Series are retained for manual investigation to determine whether an associated publication exists and whether the dataset satisfies the predefined study eligibility criteria.
