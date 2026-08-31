# TRACE CpG Amendment Log (OSF pre-registration)

**Project:** TRACE CpG: cross-tissue DNA methylation concordance between human brain and peripheral tissues  
**OSF registration:** [insert DOI/URL]  
**Original registration date:** [DD Month YYYY]

## Purpose

This document provides a prospective record of substantive amendments to the
preregistered TRACE CpG protocol following initial registration.

Amendments are documented to distinguish analyses and methodological decisions
specified in the original preregistration from those introduced subsequently.
For each amendment, we record the date, section(s) affected, nature of the
change, rationale for the change, and whether the amendment was made before or
after inspection of the relevant study results/data.

Minor editorial changes that do not alter study eligibility, data extraction,
variable definitions, statistical analyses, or interpretation are not recorded.

---

## Amendment 001

**Date:** [31 August 2026]

**Protocol section(s):** Eligibility criteria; data acquisition; data extraction

**Type of amendment:** Expansion/clarification of eligible data sources

**Original approach:**  
The initial draft protocol primarily anticipated extraction of published CpG-level
brain–peripheral correlation or concordance statistics from manuscripts and
supplementary materials ONLY.

**Amendment:**  
Where suitable CpG-level concordance statistics are not reported directly,
eligible studies may also contribute raw or processed genome-wide DNAm data
from an associated public repository (e.g., GEO), provided that sufficient
metadata are available to identify matched brain and peripheral samples from
the same individuals and CpG-level concordance can be reproducibly derived.

The least computationally intensive suitable data representation will be used
preferentially. In order of preference:

1. Complete CpG-level correlation/concordance results;
2. Processed CpG × sample methylation matrices (e.g., β-values, M-values, or
   site-level methylation fractions);
3. Raw or minimally processed platform-specific data (e.g., IDAT files for
   methylation arrays or raw signal/modified-base alignment data for long-read
   sequencing), only where necessary and computationally feasible.

**Rationale:**  
Prototype database analyses conducted in 08/26 demonstrated differences in the
reporting scope of available published datasets. Some studies provide
genome-wide/non-result-selected CpG-level correlations, whereas others report
only significance-filtered or deliberately curated CpG subsets. Restricting
TRACE CpG to published supplementary correlation tables could therefore
introduce substantial reporting/selection bias and prevent recovery of
genome-wide concordance estimates where the underlying matched-tissue
methylation data are publicly available.

**Timing relative to data/results:**  
This amendment was motivated by inspection of the reporting structure and
distribution of CpG-level results during database prototyping. It was made
before the final planned analyses of the completed TRACE CpG database.

**Potential impact:**  
The amendment may increase the number of CpGs and studies for which
non-result-selected concordance estimates can be obtained and improve
comparability across datasets. Results derived de novo from underlying
methylation data will be distinguishable from correlation estimates extracted
directly from published sources.

---

## Amendment 002

**Date:** [31 August 2026]

**Protocol section(s):** Data harmonisation; statistical analysis

**Type of amendment:** Classification of source-data reporting scope

**Original approach:**  
The original draft OSF protocol did not formally distinguish CpG-level datasets according
to whether the reported loci represented the complete measured CpG set or a
selected subset.

**Amendment:**  
Each contributing CpG-level dataset will be classified as:

1. genome-wide/non-result-selected;
2. significance-filtered; or
3. otherwise curated/selected.

This classification will be retained as a database variable and incorporated
into descriptive analyses.

**Rationale:**  
Initial analyses demonstrated that distributions of published CpG correlations
can be strongly determined by the criteria used to select CpGs for reporting.
For example, a dataset containing only statistically significant CpGs cannot
be interpreted as representing the genome-wide distribution of cross-tissue
DNAm concordance.

**Timing relative to data/results:**  
Introduced following prototype visualisation of the available CpG-level
datasets and before final analysis.

**Potential impact:**  
Between-study descriptive comparisons will explicitly account for differences
in reporting scope. Significance-filtered and curated datasets will not be
interpreted as unbiased estimates of genome-wide concordance.

---

## Amendment 003

**Date:** [DD Month YYYY]

**Protocol section(s):** Analysis plan; reproducibility analyses

**Type of amendment:** [Clarification / addition]

**Original approach:**  
[Describe original preregistered approach.]

**Amendment:**  
[Describe exactly what has changed.]

**Rationale:**  
[Explain why the change was necessary.]

**Timing relative to data/results:**  
[State what data/results had been inspected when this decision was made.]

**Potential impact:**  
[State how the amendment could affect eligibility, analyses, results or
interpretation.]

---

# Amendment summary

| ID | Date | Section(s) | Amendment | Result-informed? |
|----|------|------------|-----------|------------------|
| 001 | [date] | Eligibility; data acquisition | Permitted derivation of concordance from repository DNAm data | Yes — prompted by prototype assessment of reporting scope |
| 002 | [date] | Harmonisation; analysis | Added source-data reporting-scope classification | Yes — prompted by prototype distributions |
| 003 | [date] | — | — | — |

---

## Versioning policy

Substantive amendments will be added sequentially and will not overwrite
previous entries. Where an amendment affects an analysis, the final manuscript
will distinguish the preregistered analysis from the amended or additional
analysis where relevant.

Exploratory analyses not specified in either the original preregistration or a
prospectively documented amendment will be identified as exploratory in the
resulting manuscript.
