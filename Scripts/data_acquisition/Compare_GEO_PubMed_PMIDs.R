###############################################################
# Compare_GEO_PubMed_PMIDs.R
#
# Purpose:
# Compare PubMed IDs linked to candidate GEO Series against
# records already identified in the PubMed screening pipeline.
#
# Inputs:
#   GEO_manual_eligibility_review.csv
#   PubMed_screening_manual_review.xlsx
#
# Outputs:
#   GEO_PubMed_PMID_crossreference.csv
#   GEO_PMIDs_already_in_PubMed_search.csv
#   GEO_PMIDs_not_in_PubMed_search.csv
###############################################################

library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)


###############################################################
# 1. Load files
###############################################################

geo <- read_csv(
  "GEO_manual_eligibility_review.csv",
  show_col_types = FALSE
)

pubmed <- read_excel(
  "PubMed_screening_manual_review.xlsx",
  sheet = "Manual screening"
)


###############################################################
# 2. Standardise PubMed PMIDs
###############################################################

pubmed_pmids <- pubmed %>%
  transmute(
    PMID = as.character(PMID)
  ) %>%
  mutate(
    PMID = str_trim(PMID)
  ) %>%
  filter(
    !is.na(PMID),
    PMID != ""
  ) %>%
  distinct(PMID)


###############################################################
# 3. Extract linked PMIDs from GEO
#
# A GEO Series may contain more than one linked PMID.
###############################################################

geo_pmids <- geo %>%
  select(
    GSE,
    GEO_Title = Title,
    Linked_PMIDs,
    Has_linked_PubMed,
    GEO_triage
  ) %>%
  filter(
    !is.na(Linked_PMIDs),
    Linked_PMIDs != ""
  ) %>%
  separate_rows(
    Linked_PMIDs,
    sep = ";\\s*"
  ) %>%
  mutate(
    PMID = str_trim(
      as.character(Linked_PMIDs)
    )
  ) %>%
  filter(
    PMID != ""
  ) %>%
  select(
    -Linked_PMIDs
  )


###############################################################
# 4. Compare GEO PMIDs against PubMed search
###############################################################

comparison <- geo_pmids %>%
  mutate(
    Found_in_PubMed_search =
      PMID %in% pubmed_pmids$PMID
  ) %>%
  left_join(
    pubmed %>%
      mutate(
        PMID = str_trim(
          as.character(PMID)
        )
      ) %>%
      select(
        PMID,
        PubMed_Title = Title,
        DOI,
        Reference,
        Year,
        Final_verdict,
        Final_gate,
        Final_reason
      ),
    by = "PMID"
  ) %>%
  arrange(
    desc(Found_in_PubMed_search),
    GSE,
    PMID
  )


###############################################################
# 5. Separate matched and unmatched GEO publications
###############################################################

geo_found_in_pubmed <- comparison %>%
  filter(
    Found_in_PubMed_search
  )

geo_not_found_in_pubmed <- comparison %>%
  filter(
    !Found_in_PubMed_search
  )


###############################################################
# 6. Export results
###############################################################

write_csv(
  comparison,
  "GEO_PubMed_PMID_crossreference.csv"
)

write_csv(
  geo_found_in_pubmed,
  "GEO_PMIDs_already_in_PubMed_search.csv"
)

write_csv(
  geo_not_found_in_pubmed,
  "GEO_PMIDs_not_in_PubMed_search.csv"
)


###############################################################
# 7. Report results
###############################################################

cat(
  "\n=========================================================\n",
  "GEO–PubMed PMID Cross-reference\n",
  "=========================================================\n\n",
  
  "GEO Series screened:                 ",
  n_distinct(geo$GSE), "\n",
  
  "GEO Series with linked PMID(s):      ",
  n_distinct(geo_pmids$GSE), "\n",
  
  "Unique GEO-linked PMIDs:             ",
  n_distinct(geo_pmids$PMID), "\n\n",
  
  "GEO-linked PMIDs found in PubMed:    ",
  n_distinct(geo_found_in_pubmed$PMID), "\n",
  
  "GEO-linked PMIDs NOT found:          ",
  n_distinct(geo_not_found_in_pubmed$PMID), "\n\n",
  
  "GEO Series represented in PubMed:    ",
  n_distinct(geo_found_in_pubmed$GSE), "\n",
  
  "=========================================================\n",
  sep = ""
)


###############################################################
# 8. Print publications missed by PubMed search
###############################################################

if (nrow(geo_not_found_in_pubmed) > 0) {
  
  cat(
    "\nGEO-linked publications NOT identified by the ",
    "PubMed search:\n\n"
  )
  
  print(
    geo_not_found_in_pubmed %>%
      select(
        GSE,
        PMID,
        GEO_Title,
        GEO_triage
      )
  )
  
} else {
  
  cat(
    "\nAll GEO-linked PMIDs were already identified ",
    "by the PubMed search.\n"
  )
}


###############################################################
# End
###############################################################
