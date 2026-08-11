###############################################################
# Extract_pubmed_screening.R
#
# Purpose:
# Run the preregistered PubMed search, retrieve citation
# metadata and abstracts, and create an Excel screening file
# for manual review of eligibility by independent screeners
###############################################################

#install.packages("rentrez")
#install.packages("openxlsx")

###############################################################
# Start timer
###############################################################

start_time <- Sys.time()

library(rentrez)
library(xml2)
library(dplyr)
library(stringr)
library(purrr)
library(openxlsx)

###############################################################
# 1. Define PubMed search
###############################################################

pubmed_query <- '

(
  "DNA methylation"[Title/Abstract]
  OR DNAm[Title/Abstract]
  OR methylome[Title/Abstract]
  OR EWAS[Title/Abstract]
  OR epigenome-wide[Title/Abstract]
)

AND
(
  blood[Title/Abstract]
  OR saliva[Title/Abstract]
  OR buccal[Title/Abstract]
)

AND 
( 
  brain[Title/Abstract] 
  OR "brain tissue"[Title/Abstract] 
  OR cortex[Title/Abstract] 
  OR "cerebral cortex"[Title/Abstract] 
  OR "prefrontal cortex"[Title/Abstract] 
  OR "dorsolateral prefrontal cortex"[Title/Abstract] 
  OR DLPFC[Title/Abstract] 
  OR "entorhinal cortex"[Title/Abstract] 
  OR "superior temporal gyrus"[Title/Abstract] 
  OR STG[Title/Abstract] 
  OR "inferior frontal gyrus"[Title/Abstract] 
  OR cerebellum[Title/Abstract] 
  OR cerebellar[Title/Abstract] 
  OR hippocampus[Title/Abstract] 
  OR hippocampal[Title/Abstract] 
  OR amygdala[Title/Abstract] 
  OR "nucleus accumbens"[Title/Abstract]
  OR striatum[Title/Abstract]
  OR putamen[Title/Abstract]
  OR caudate[Title/Abstract]
  OR "substantia nigra"[Title/Abstract]
)



'

###############################################################
# 2. Search PubMed
###############################################################

search_results <- entrez_search(
  db = "pubmed",
  term = pubmed_query,
  retmax = 10000,
  use_history = TRUE
)

cat("PubMed records identified:", search_results$count, "\n")

###############################################################
# 3. Retrieve full PubMed XML
###############################################################

pubmed_xml <- entrez_fetch(
  db = "pubmed",
  web_history = search_results$web_history,
  rettype = "xml",
  parsed = FALSE
)

doc <- read_xml(pubmed_xml)

articles <- xml_find_all(doc, ".//PubmedArticle")

###############################################################
# 4. Helper functions
###############################################################

get_text <- function(node, xpath) {
  x <- xml_find_first(node, xpath)
  
  if (inherits(x, "xml_missing")) {
    return(NA_character_)
  }
  
  text <- xml_text(x, trim = TRUE)
  
  ifelse(text == "", NA_character_, text)
}

get_abstract <- function(node) {
  
  abstract_nodes <- xml_find_all(
    node,
    ".//Article/Abstract/AbstractText"
  )
  
  if (length(abstract_nodes) == 0) {
    return(NA_character_)
  }
  
  paste(
    xml_text(abstract_nodes, trim = TRUE),
    collapse = " "
  )
}

get_authors <- function(node) {
  
  authors <- xml_find_all(
    node,
    ".//Article/AuthorList/Author"
  )
  
  if (length(authors) == 0) {
    return(NA_character_)
  }
  
  author_names <- map_chr(
    authors,
    function(a) {
      
      surname <- get_text(a, "./LastName")
      initials <- get_text(a, "./Initials")
      
      paste(
        na.omit(c(surname, initials)),
        collapse = " "
      )
    }
  )
  
  paste(author_names, collapse = "; ")
}

get_first_author <- function(node) {
  
  author <- xml_find_first(
    node,
    ".//Article/AuthorList/Author[1]"
  )
  
  if (inherits(author, "xml_missing")) {
    return(NA_character_)
  }
  
  surname <- get_text(author, "./LastName")
  initials <- get_text(author, "./Initials")
  
  paste(
    na.omit(c(surname, initials)),
    collapse = " "
  )
}

get_doi <- function(node) {
  
  doi_nodes <- xml_find_all(
    node,
    './/ArticleIdList/ArticleId[@IdType="doi"]'
  )
  
  if (length(doi_nodes) == 0) {
    return(NA_character_)
  }
  
  xml_text(doi_nodes[1], trim = TRUE)
}

###############################################################
# 5. Extract citation information
###############################################################

pubmed_records <- map_dfr(
  articles,
  function(article) {
    
    tibble(
      PMID = get_text(
        article,
        ".//MedlineCitation/PMID"
      ),
      
      DOI = get_doi(article),
      
      First_author = get_first_author(article),
      
      Authors = get_authors(article),
      
      Year = coalesce(
        get_text(
          article,
          ".//Article/Journal/JournalIssue/PubDate/Year"
        ),
        get_text(
          article,
          ".//Article/ArticleDate/Year"
        )
      ),
      
      Journal = get_text(
        article,
        ".//Article/Journal/Title"
      ),
      
      Title = get_text(
        article,
        ".//Article/ArticleTitle"
      ),
      
      Abstract = get_abstract(article)
    )
  }
)

###############################################################
# 6. Construct compact study citation
###############################################################

pubmed_records <- pubmed_records %>%
  mutate(
    Study = paste0(
      First_author,
      " et al. ",
      Year,
      ", ",
      Journal
    )
  )

###############################################################
# 7. Add useful automatic screening flags
#
# IMPORTANT:
# These are prioritisation aids only.
# They do NOT determine eligibility.
###############################################################

pubmed_records <- pubmed_records %>%
  mutate(
    
    search_text = str_to_lower(
      paste(Title, Abstract)
    ),
    
    mentions_brain = str_detect(
      search_text,
      paste(
        c(
          "brain",
          "cortex",
          "cortical",
          "prefrontal",
          "entorhinal",
          "temporal gyrus",
          "cerebell",
          "hippocamp"
        ),
        collapse = "|"
      )
    ),
    
    mentions_peripheral = str_detect(
      search_text,
      paste(
        c(
          "whole blood",
          "peripheral blood",
          "blood",
          "saliva",
          "buccal",
          "leukocyte"
        ),
        collapse = "|"
      )
    ),
    
    mentions_matching = str_detect(
      search_text,
      paste(
        c(
          "matched",
          "paired",
          "same individuals",
          "same subjects",
          "same donors",
          "cross-tissue",
          "cross tissue",
          "concordance",
          "correlation",
          "surrogate",
          "blood-brain",
          "brain-blood"
        ),
        collapse = "|"
      )
    ),
    
    mentions_genomewide = str_detect(
      search_text,
      paste(
        c(
          "epigenome-wide",
          "epigenome wide",
          "genome wide",
          "genome-wide",
         # "genome-wide DNA methylation",
         # "genome wide DNA methylation",
         # "genome-wide DNAm",
         # "genome wide DNAm",
         # "genome-wide methylation",
         # "genome wide methylation",
          "450k",
          "epic array",
          "methylation array"
        ),
        collapse = "|"
      )
    )
  )

###############################################################
# 8. Create screening columns
###############################################################

screening <- pubmed_records %>%
  select(
    DOI,
    PMID,
    Study,
    Title,
    Year,
    Journal,
    
    mentions_brain,
    mentions_peripheral,
    mentions_matching,
    mentions_genomewide,
    
    Abstract
  ) %>%
  mutate(
    
    Verdict = "",
    
    Reason = "",
    
    Human_samples = "",
    
    Brain_tissue = "",
    
    Peripheral_tissue = "",
    
    Matched_samples = "",
    
    Genome_wide_DNAm = "",
    
    CpG_level_concordance = "",
    
    Downloadable_results = "",
    
    Full_text_checked = "",
    
    Notes = ""
  )

###############################################################
# 9. Remove any duplicate PubMed records
###############################################################

screening <- screening %>%
  distinct(PMID, .keep_all = TRUE)


###############################################################
# 9b. Identify no. eligible studies from search
###############################################################


screening_possibily_eligible <- screening %>%
  filter(mentions_brain == "TRUE" &
           mentions_peripheral == "TRUE" &
              mentions_matching == "TRUE" &
                mentions_genomewide == "TRUE")


# n = 129 are possibly eligible studies

library(dplyr)
library(stringr)


#Sanity check

nishitani_study <- screening %>%
  filter(
    str_detect(Study, regex("Nishitani", ignore_case = TRUE)) |
      str_detect(DOI, regex("10\\.1038/s41398-023-02370-0", ignore_case = TRUE))
  )

nishitani_study



###############################################################
# Combine title and abstract into one searchable field
###############################################################

screening_possibily_eligible_exc <- screening_possibily_eligible %>%
  mutate(
    search_text = str_to_lower(
      paste(
        coalesce(Title, ""),
        coalesce(Abstract, "")
      )
    )
  )

###############################################################
# Define exclusion / screening flags
###############################################################

screening_possibily_eligible_exc <- screening_possibily_eligible_exc %>%
  mutate(
    
    flag_review = str_detect(
      search_text,
      regex(
        "\\breview\\b|\\bsystematic review\\b|\\bscoping review\\b",
        ignore_case = TRUE
      )
    ),
    
    flag_meta_analysis = str_detect(
      search_text,
      regex(
        "\\bmeta[- ]?analysis\\b",
        ignore_case = TRUE
      )
    ),
    
    flag_non_human = str_detect(
      search_text,
      regex(
        "\\bmouse\\b",
        ignore_case = TRUE
      )
    ),
    
    flag_candidate_gene = str_detect(
      search_text,
      regex(
        "\\bcandidate[- ]?gene\\b|\\bgene[- ]specific\\b",
        ignore_case = TRUE
      )
    ),
    
    flag_differential_methylation = str_detect(
      search_text,
      regex(
        "differential methylation|differentially methylated|\\bDMPs?\\b|\\bDMRs?\\b",
        ignore_case = TRUE
      )
    )
  )


screening_possibily_eligible_exc  <- screening_possibily_eligible_exc  %>%
  mutate(
    exclusion_flag = case_when(
      flag_review ~ "Review",
      flag_non_human ~ "Non human sample",
      flag_meta_analysis ~ "Meta-analysis",
      flag_candidate_gene ~ "Candidate gene",
      flag_differential_methylation ~ "Differential methylation",
      TRUE ~ ""
    )
  )


screening_possibily_eligible_exc <- screening_possibily_eligible_exc %>%
  rowwise() %>%
  mutate(
    exclusion_flags = paste(
      c(
        if (flag_review) "Review" else NULL,
        if (flag_meta_analysis) "Meta-analysis" else NULL,
        if (flag_candidate_gene) "Candidate gene" else NULL,
        if (flag_differential_methylation) "Differential methylation" else NULL,
        if (flag_non_human) "EWAS" else NULL
      ),
      collapse = "; "
    )
  ) %>%
  ungroup()

###

screening_possibily_eligible_EXCLUDED <- screening_possibily_eligible_exc %>%
  filter(exclusion_flags != "") %>%
  select(
    DOI,	
    PMID,
    Study,
    Title,
    Year,
    Journal,
    Abstract,
    Verdict,
    Reason,
    exclusion_flags
  )

# leaves us with n = 49 to exclude, leaving 35 to manually screen

screening_possibily_eligible_MANUAL <- screening_possibily_eligible_exc %>%
  mutate(
    Verdict = if_else(
      exclusion_flags != "",
      "Exclude",
      ""
    ),
    Reason = if_else(
      exclusion_flags != "",
      exclusion_flags,
      ""
    )
  ) %>%
  select(
    DOI,	
    PMID,
    Study,
    Title,
    Year,
    Journal,
    Abstract,
    Verdict,
    Reason
  )


write.csv(
  screening_possibily_eligible_MANUAL,
  "PubMed_screening_manual.csv",
  row.names = FALSE
)

###############################################################
# 10. Save raw PubMed retrieval
###############################################################

write.csv(
  screening,
  "PubMed_screening_records.csv",
  row.names = FALSE
)

###############################################################
# 11. Create formatted Excel workbook
###############################################################

wb <- createWorkbook()

addWorksheet(
  wb,
  "Screening"
)

writeData(
  wb,
  sheet = "Screening",
  x = screening,
  withFilter = TRUE
)

freezePane(
  wb,
  sheet = "Screening",
  firstRow = TRUE
)

setColWidths(
  wb,
  sheet = "Screening",
  cols = 1:ncol(screening),
  widths = "auto"
)

# Keep long fields manageable
setColWidths(
  wb,
  sheet = "Screening",
  cols = which(names(screening) == "Title"),
  widths = 50
)

setColWidths(
  wb,
  sheet = "Screening",
  cols = which(names(screening) == "Abstract"),
  widths = 80
)

setColWidths(
  wb,
  sheet = "Screening",
  cols = which(names(screening) == "Reason"),
  widths = 45
)

setColWidths(
  wb,
  sheet = "Screening",
  cols = which(names(screening) == "Notes"),
  widths = 45
)

wrapStyle <- createStyle(
  wrapText = TRUE,
  valign = "top"
)

addStyle(
  wb,
  sheet = "Screening",
  style = wrapStyle,
  rows = 2:(nrow(screening) + 1),
  cols = 1:ncol(screening),
  gridExpand = TRUE
)

###############################################################
# 12. Add dropdowns for screening decisions
###############################################################

dataValidation(
  wb,
  sheet = "Screening",
  cols = which(names(screening) == "Verdict"),
  rows = 2:(nrow(screening) + 1),
  type = "list",
  value = '"Include,Exclude,Maybe,Full-text screen"'
)

yes_no_cols <- which(
  names(screening) %in% c(
    "Human_samples",
    "Brain_tissue",
    "Peripheral_tissue",
    "Matched_samples",
    "Genome_wide_DNAm",
    "CpG_level_concordance",
    "Downloadable_results",
    "Full_text_checked"
  )
)

for (col in yes_no_cols) {
  
  dataValidation(
    wb,
    sheet = "Screening",
    cols = col,
    rows = 2:(nrow(screening) + 1),
    type = "list",
    value = '"Yes,No,Unclear"'
  )
}

###############################################################
# 13. Add search information sheet
###############################################################

addWorksheet(
  wb,
  "Search_information"
)

search_info <- data.frame(
  Field = c(
    "Database",
    "Search date",
    "Number retrieved",
    "Search string"
  ),
  
  Value = c(
    "PubMed",
    as.character(Sys.Date()),
    search_results$count,
    pubmed_query
  )
)

writeData(
  wb,
  "Search_information",
  search_info
)

setColWidths(
  wb,
  "Search_information",
  cols = 1:2,
  widths = c(20, 100)
)

###############################################################
# 14. Save workbook
###############################################################

saveWorkbook(
  wb,
  "PubMed_screening.xlsx",
  overwrite = TRUE
)

cat(
  "Screening workbook created with",
  nrow(screening),
  "unique PubMed records.\n"
)

###############################################################
# Report runtime
###############################################################

end_time <- Sys.time()

runtime <- as.numeric(
  difftime(end_time, start_time, units = "secs")
)

hours <- runtime %/% 3600
minutes <- (runtime %% 3600) %/% 60
seconds <- round(runtime %% 60)

cat(
  "\n=================================================\n",
  "PubMed screening completed successfully\n",
  "Unique PubMed records:",
  nrow(screening), "\n",
  "Runtime:",
  sprintf("%02d:%02d:%02d", hours, minutes, seconds),
  "\n=================================================\n"
)


#####
