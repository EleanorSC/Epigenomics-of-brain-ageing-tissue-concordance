###############################################################
# Extract_pubmed_screening.R
#
# Purpose:
# Run the preregistered lit. search, retrieve citation
# metadata and abstracts, apply conservative automated
# screening flags, organise exclusions by prespecified
# screening gate, and create CSV files for manual review.
#
# Screening gates
# Gate 1: Appropriate study design
# Gate 2: Appropriate tissue sampling
# Gate 3: Appropriate analysis
# Gate 4: Appropriate data availability
#
# IMPORTANT:
# Automated exclusions are assigned only where the PubMed
# title/abstract/publication type provides positive evidence.
# Absence of a keyword is NOT treated as evidence of exclusion.
# Gate 2 can include conservative automatic exclusions where the
# title/abstract gives evidence of a study ill-fitted to main aims.
# Gate 4 will usually require manual full-text review.
###############################################################

# install.packages(c(
#   "rentrez", "xml2", "dplyr", "stringr",
#   "purrr"
# ))

###############################################################
# Start timer
###############################################################

start_time <- Sys.time()

cat(
  "\n=========================================================\n",
  "Brain–Peripheral DNA Methylation Concordance Project\n",
  "Gate-based PubMed Screening Pipeline\n",
  "Started: ", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n",
  "=========================================================\n",
  sep = ""
)

###############################################################
# Progress reporting 
###############################################################

stage_start <- function(stage_number, stage_name) {
  cat(
    "\n---------------------------------------------------------\n",
    sprintf("[%d/9] %s\n", stage_number, stage_name),
    "Started: ", format(Sys.time(), "%H:%M:%S"), "\n",
    "---------------------------------------------------------\n",
    sep = ""
  )
  Sys.time()
}

stage_end <- function(stage_time, message = "Completed") {
  elapsed <- as.numeric(
    difftime(Sys.time(), stage_time, units = "secs")
  )

  cat(
    "✓ ", message,
    " (", round(elapsed, 1), " s)\n",
    sep = ""
  )
}

###############################################################
# Load packages
###############################################################

library(rentrez)
library(xml2)
library(dplyr)
library(stringr)
library(purrr)

###############################################################
# 1. Define search terms
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
  OR serum[Title/Abstract]
  OR placenta[Title/Abstract]
  OR plasma[Title/Abstract]
  OR saliva[Title/Abstract]
  OR buccal[Title/Abstract]
  OR cheek[Title/Abstract]
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
# 2. Search PubMed - as first test-case
###############################################################

t <- stage_start(1, "Searching PubMed")

search_results <- entrez_search(
  db = "pubmed",
  term = pubmed_query,
  retmax = 10000,
  use_history = TRUE
)

stage_end(
  t,
  paste0(
    format(search_results$count, big.mark = ","),
    " records identified"
  )
)

###############################################################
# 3. Retrieve PubMed XML in batches
#
# PubMed records are downloaded in batches rather than as one
# large transfer. This reduces the chance that a transient NCBI
# or network interruption terminates my retrieval attempt.
###############################################################

t <- stage_start(2, "Downloading PubMed records")

batch_size <- 500
max_retries <- 3
n_records <- search_results$count

batch_starts <- seq(
  from = 0,
  to = n_records - 1,
  by = batch_size
)

pubmed_xml_batches <- vector(
  "list",
  length(batch_starts)
)

for (i in seq_along(batch_starts)) {

  retstart <- batch_starts[i]

  retmax <- min(
    batch_size,
    n_records - retstart
  )

  cat(
    "  Downloading batch ",
    i,
    "/",
    length(batch_starts),
    " (records ",
    retstart + 1,
    "-",
    retstart + retmax,
    ")...\n",
    sep = ""
  )

  success <- FALSE

  for (attempt in seq_len(max_retries)) {

    result <- tryCatch(
      {
        entrez_fetch(
          db = "pubmed",
          web_history = search_results$web_history,
          rettype = "xml",
          parsed = FALSE,
          retstart = retstart,
          retmax = retmax
        )
      },
      error = function(e) {

        cat(
          "    Attempt ",
          attempt,
          " failed: ",
          conditionMessage(e),
          "\n",
          sep = ""
        )

        NULL
      }
    )

    if (!is.null(result)) {

      pubmed_xml_batches[[i]] <- result
      success <- TRUE

      break
    }

    if (attempt < max_retries) {

      retry_delay <- 2 * attempt

      cat(
        "    Retrying in ",
        retry_delay,
        " seconds...\n",
        sep = ""
      )

      Sys.sleep(retry_delay)
    }
  }

  if (!success) {
    stop(
      paste0(
        "Failed to download PubMed batch ",
        i,
        " after ",
        max_retries,
        " attempts."
      )
    )
  }

  # Brief delay between successful requests
  Sys.sleep(0.4)
}

stage_end(
  t,
  paste0(
    length(pubmed_xml_batches),
    " PubMed batches downloaded"
  )
)

###############################################################
# 4. Parse PubMed XML
###############################################################

t <- stage_start(3, "Parsing PubMed XML")

articles <- map(
  pubmed_xml_batches,
  function(x) {

    doc <- read_xml(x)

    xml_find_all(
      doc,
      ".//PubmedArticle"
    )
  }
) %>%
  flatten()

###############################################################
# Retrieval integrity check
###############################################################

if (length(articles) != search_results$count) {

  warning(
    paste0(
      "PubMed retrieval count mismatch: search returned ",
      search_results$count,
      " records but ",
      length(articles),
      " records were downloaded and parsed."
    )
  )

} else {

  cat(
    "✓ Retrieval integrity check passed: ",
    format(length(articles), big.mark = ","),
    " / ",
    format(search_results$count, big.mark = ","),
    " records retrieved\n",
    sep = ""
  )
}

stage_end(
  t,
  paste0(
    format(length(articles), big.mark = ","),
    " articles parsed"
  )
)

###############################################################
# 5. Helper functions
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

get_publication_types <- function(node) {
  type_nodes <- xml_find_all(
    node,
    ".//Article/PublicationTypeList/PublicationType"
  )

  if (length(type_nodes) == 0) {
    return(NA_character_)
  }

  paste(
    xml_text(type_nodes, trim = TRUE),
    collapse = "; "
  )
}

###############################################################
# 6. Extract citation information
###############################################################

t <- stage_start(4, "Extracting PubMed citation metadata")

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

      Abstract = get_abstract(article),

      Publication_types = get_publication_types(article)
    )
  }
)

stage_end(
  t,
  paste0(
    format(nrow(pubmed_records), big.mark = ","),
    " records extracted"
  )
)

###############################################################
# 7. search/prioritisation fields insofar as I can automate this
###############################################################

t <- stage_start(5, "Creating screening and prioritisation fields")

pubmed_records <- pubmed_records %>%
  mutate(
    Study = paste0(
      First_author,
      " et al. ",
      Year,
      ", ",
      Journal
    ),

    search_text = str_to_lower(
      paste(
        coalesce(Title, ""),
        coalesce(Abstract, ""),
        coalesce(Publication_types, "")
      )
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
          "hippocamp",
          "amygdala",
          "striat",
          "putamen",
          "caudate",
          "substantia nigra"
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
          "serum",
          "plasma",
          "placenta",
          "saliva",
          "buccal",
          "cheek",
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
          "450k",
          "450 k",
          "epic array",
          "methylation array",
          "infinium"
        ),
        collapse = "|"
      )
    )
  )

stage_end(
  t,
  paste0(
    "screening fields created for ",
    format(nrow(pubmed_records), big.mark = ","),
    " records"
  )
)

###############################################################
# 8. Deduplicate to stop bug
###############################################################

t <- stage_start(6, "Deduplicating PubMed records")

n_before_dedup <- nrow(pubmed_records)

screening <- pubmed_records %>%
  distinct(PMID, .keep_all = TRUE)

n_duplicates_removed <- n_before_dedup - nrow(screening)

stage_end(
  t,
  paste0(
    format(nrow(screening), big.mark = ","),
    " unique records retained; ",
    format(n_duplicates_removed, big.mark = ","),
    " duplicates removed"
  )
)

###############################################################
# 9. Apply conservative automated gate exclusions
###############################################################

t <- stage_start(7, "Applying gate-based automated screening")

screening <- screening %>%
  mutate(

    ###########################################################
    # Gate 1: Appropriate study design?
    ###########################################################

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

    flag_candidate_gene = str_detect(
      search_text,
      regex(
        "\\bcandidate[- ]?gene\\b|\\bgene[- ]specific\\b",
        ignore_case = TRUE
      )
    ),

    flag_non_human = str_detect(
      search_text,
      regex(
        "\\bmouse\\b|\\bmice\\b|\\brat\\b|\\brats\\b|\\bmonkey\\b|\\bmonkeys\\b|\\bmacaque\\b|\\bnon[- ]human primate",
        ignore_case = TRUE
      )
    ),

    flag_abstract_only = str_detect(
      str_to_lower(
        paste(
          coalesce(Publication_types, ""),
          coalesce(Title, "")
        )
      ),
      regex(
        "conference abstract|meeting abstract|congress abstract|conference proceedings|meeting proceedings",
        ignore_case = TRUE
      )
    ),

    flag_secondary_analysis = str_detect(
      search_text,
      regex(
        "secondary analysis|secondary analyses|previously published dataset|previously published data|prediction model|predictive model|risk score",
        ignore_case = TRUE
      )
    ),
    
    ###########################################################
    # Gate 2: Appropriate tissue sampling?
    #
    # These rules use positive evidence from the title/abstract.
    # Absence of matching terminology alone is NOT treated as
    # sufficient evidence for exclusion.
    ###########################################################

    mentions_brain_tissue_dnam = str_detect(
      search_text,
      regex(
        paste(
          c(
            "brain tissue.{0,80}methyl",
            "brain.{0,80}DNA methylation",
            "brain.{0,80}DNAm",
            "cortex.{0,80}methyl",
            "cortical.{0,80}methyl",
            "prefrontal cortex.{0,80}methyl",
            "entorhinal cortex.{0,80}methyl",
            "cerebell.{0,80}methyl",
            "hippocamp.{0,80}methyl",
            "amygdala.{0,80}methyl",
            "striat.{0,80}methyl",
            "postmortem brain.{0,80}methyl",
            "post-mortem brain.{0,80}methyl"
          ),
          collapse = "|"
        ),
        ignore_case = TRUE
      )
    ),

    mentions_neuroimaging = str_detect(
      search_text,
      regex(
        paste(
          c(
            "\\bMRI\\b",
            "\\bfMRI\\b",
            "\\bDTI\\b",
            "\\bEEG\\b",
            "magnetic resonance imaging",
            "neuroimaging",
            "brain imaging",
            "cortical thickness",
            "cortical volume",
            "grey matter volume",
            "gray matter volume",
            "white matter integrity",
            "fractional anisotropy",
            "mean diffusivity",
            "\\bPET\\b",
            "positron emission tomography"
          ),
          collapse = "|"
        ),
        ignore_case = TRUE
      )
    ),

    flag_neuroimaging_no_brain_dnam =
      mentions_neuroimaging &
      mentions_peripheral &
      !mentions_brain_tissue_dnam,

    flag_explicit_unmatched_tissues = str_detect(
      search_text,
      regex(
        paste(
          c(
            "independent cohorts",
            "separate cohorts",
            "different cohorts",
            "independent samples",
            "separate samples",
            "different individuals",
            "different participants",
            "unmatched samples",
            "unpaired samples"
          ),
          collapse = "|"
        ),
        ignore_case = TRUE
      )
    ),

    ###########################################################
    # Gate 3: Appropriate analysis?
    #
    # Differential methylation is SOMETIMES identifiable from
    # title/abstract. BUT Lack of concordance terminology is NOT
    # treated as automatic exclusion.
    ###########################################################

    flag_differential_methylation = str_detect(
      search_text,
      regex(
        "differential methylation|differentially methylated|\\bDMPs?\\b|\\bDMRs?\\b",
        ignore_case = TRUE
      )
    ),

    ###########################################################
    # Automatic gate assignment
    #
    # Gate 1 takes precedence over Gate 3.
    # Gate 2 and Gate 4 normally require manual assessment.
    ###########################################################

    Auto_gate = case_when(
      flag_review ~ "Gate 1",
      flag_meta_analysis ~ "Gate 1",
      flag_candidate_gene ~ "Gate 1",
      flag_secondary_analysis ~ "Gate 1",
      flag_abstract_only ~ "Gate 1",
      flag_non_human ~ "Gate 1",
      flag_neuroimaging_no_brain_dnam ~ "Gate 2",
      flag_explicit_unmatched_tissues ~ "Gate 2",
      flag_differential_methylation ~ "Gate 3",
      TRUE ~ ""
    ),

    Auto_reason = case_when(
      flag_review ~ "Review article",
      flag_meta_analysis ~ "Meta-analysis",
      flag_candidate_gene ~ "Candidate-gene study",
      flag_secondary_analysis ~ "Secondary analysis of previously published data",
      flag_abstract_only ~ "Abstract only",
      flag_non_human ~ "Non-human tissue/sample",
      flag_neuroimaging_no_brain_dnam ~ "Neuroimaging phenotype; no brain-tissue DNAm",
      flag_explicit_unmatched_tissues ~ "Brain and peripheral samples explicitly unmatched",
      flag_differential_methylation ~ "Differential methylation only / DMP-DMR analysis",
      TRUE ~ ""
    ),

    Auto_verdict = if_else(
      Auto_gate != "",
      "Exclude",
      ""
    )
  )



## Sanity check, review all studies that have been flagged as non-human
## n = 247 studies 
non_human_studies <- screening %>%
  filter(flag_non_human == TRUE) %>%
  select(
    PMID,
    DOI,
    Study,
    Title,
    Abstract,
    Auto_gate,
    Auto_reason
  )

non_human_studies

## Sanity check, review all studies that have been flagged as unmatched tissues
## n = 28 studies 
unmatched_samples <- screening %>%
  filter(flag_explicit_unmatched_tissues == TRUE) %>%
  select(
    PMID,
    DOI,
    Study,
    Title,
    Abstract,
    Auto_gate,
    Auto_reason
  )

unmatched_samples


## Sanity check, review all studies that have been flagged as secondary analyses
## n = 35 studies

secondary_analysis <- screening %>%
  filter(flag_secondary_analysis == TRUE) %>%
  select(
    PMID,
    DOI,
    Study,
    Title,
    Abstract,
    Auto_gate,
    Auto_reason
  )

## Sanity check, review all studies that have been flagged as candidate gene studies
## n = 37 studies

candidate_gene <- screening %>%
  filter(flag_candidate_gene == TRUE) %>%
  select(
    PMID,
    DOI,
    Study,
    Title,
    Abstract,
    Auto_gate,
    Auto_reason
  )


###############################################################
# Add manual screening fields
###############################################################

screening <- screening %>%
  mutate(
    Gate1_design = "",
    Gate2_tissue_sampling = "",
    Gate3_analysis = "",
    Gate4_data_availability = "",

    Manual_gate = "",
    Manual_reason = "",

    Final_verdict = Auto_verdict,
    Final_gate = Auto_gate,
    Final_reason = Auto_reason,

    Notes = ""
  )

###############################################################
# Automatically prioritised records
###############################################################

screening_priority <- screening %>%
  filter(
    Auto_verdict == "",
    mentions_brain == TRUE,
    mentions_peripheral == TRUE,
    mentions_matching == TRUE,
    mentions_genomewide == TRUE
  )

###############################################################
# Automated exclusions
###############################################################

screening_auto_excluded <- screening %>%
  filter(Auto_verdict == "Exclude")

###############################################################
# Records still (potentially) requiring manual screening
###############################################################

screening_manual_review <- screening %>%
  filter(Auto_verdict == "")

###############################################################
# Gate summary: AUTOMATED exclusions only
###############################################################

auto_gate_summary <- screening %>%
  filter(Auto_verdict == "Exclude") %>%
  count(
    Auto_gate,
    Auto_reason,
    name = "n_excluded"
  ) %>%
  arrange(
    Auto_gate,
    desc(n_excluded),
    Auto_reason
  )

auto_gate_totals <- screening %>%
  filter(Auto_verdict == "Exclude") %>%
  count(
    Auto_gate,
    name = "n_failed"
  ) %>%
  arrange(Auto_gate)

cat(
  "  Automatically excluded: ",
  format(nrow(screening_auto_excluded), big.mark = ","),
  "\n",
  "  Potential for review according to Gate: ",
  format(nrow(screening_manual_review), big.mark = ","),
  "\n",
  "  High-priority manual records: ",
  format(nrow(screening_priority), big.mark = ","),
  "\n",
  sep = ""
)

stage_end(t, "Gate-based automated screening completed")

###############################################################
# 10. Create screening summary tables
###############################################################

t <- stage_start(8, "Creating gate summary tables")

###############################################################
# Gate definitions
###############################################################

gate_definitions <- tibble(
  Gate = c(
    "Gate 1",
    "Gate 2",
    "Gate 3",
    "Gate 4"
  ),
  Question = c(
    "Is this a primary human genome-wide DNA methylation study rather than a candidate-gene study, review, meta-analysis, abstract-only report, secondary analysis, or non-human study?",
    "Does the study measure matched brain and peripheral DNA methylation in the same individuals?",
    "Does the study report CpG-level brain-peripheral correlation or concordance statistics rather than only differential methylation or conventional EWAS analyses?",
    "Are downloadable genome-wide CpG-level concordance statistics available for extraction without re-analysing raw data?"
  )
)

###############################################################
# Prespecified exclusion reasons
###############################################################

exclusion_dictionary <- tibble(
  Gate = c(
    "Gate 1",
    "Gate 1",
    "Gate 1",
    "Gate 1",
    "Gate 1",
    "Gate 1",
    "Gate 2",
    "Gate 2",
    "Gate 2",
    "Gate 3",
    "Gate 3",
    "Gate 4"
  ),

  Exclusion_reason = c(
    "Candidate-gene study",
    "Meta-analysis",
    "Secondary analysis of previously published matched brain-peripheral data",
    "Abstract only",
    "Review article",
    "Non-human tissue/sample",
    "Neuroimaging phenotype; no brain-tissue DNAm",
    "Brain and peripheral samples explicitly unmatched",
    "No matched brain and peripheral tissue DNAm in the same individuals",
    "Differential methylation only / DMP-DMR analysis",
    "No CpG-level brain-peripheral concordance/correlation analysis",
    "No downloadable CpG-level concordance results; raw-data reanalysis required"
  )
)

###############################################################
# Overall automated screening status
###############################################################

screening_status_summary <- tibble(
  Screening_status = c(
    "Records retrieved",
    "Unique PubMed records",
    "Duplicates removed",
    "Automatically excluded",
    "Remaining for manual screening",
    "High-priority manual records"
  ),

  n = c(
    search_results$count,
    nrow(screening),
    n_duplicates_removed,
    nrow(screening_auto_excluded),
    nrow(screening_manual_review),
    nrow(screening_priority)
  )
)

stage_end(t, "Gate summary tables created")

###############################################################
# 11. Write CSV outputs
###############################################################

t <- stage_start(9, "Writing screening CSV files")

write.csv(
  screening,
  "PubMed_screening_records.csv",
  row.names = FALSE
)

write.csv(
  screening_manual_review,
  "PubMed_screening_manual.csv",
  row.names = FALSE
)

write.csv(
  screening_auto_excluded,
  "PubMed_screening_automatic_exclusions.csv",
  row.names = FALSE
)

write.csv(
  auto_gate_summary,
  "PubMed_screening_gate_summary.csv",
  row.names = FALSE
)

write.csv(
  exclusion_dictionary,
  "PubMed_screening_exclusion_dictionary.csv",
  row.names = FALSE
)

write.csv(
  screening_priority,
  "PubMed_screening_priority.csv",
  row.names = FALSE
)

stage_end(
  t,
  "Screening CSV files written"
)


#library(dplyr)
#library(readr)
#library(stringr)

# ============================================================
# CREATE MANUAL SCREENING WORKBOOK
# ============================================================

#library(dplyr)
library(openxlsx)

# ------------------------------------------------------------
# 1. CREATE REVIEWER-FACING DATAFRAME
# ------------------------------------------------------------

manual_screening <- screening %>%
  mutate(
    
    # Short citation used for screening
    Reference = case_when(
      !is.na(First_author) & !is.na(Year) ~
        paste0(First_author, " et al. ", Year),
      !is.na(First_author) ~ First_author,
      TRUE ~ NA_character_
    ),
    
    # --------------------------------------------------------
    # ORDER RECORDS FOR MANUAL REVIEW
    #
    # 1 = High-priority records requiring manual screening
    # 2 = Other records requiring manual screening
    # 3 = Automatically excluded records
    #
    # This variable is used only for sorting and is removed
    # from the final workbook.
    # --------------------------------------------------------
    
    manual_review_order = case_when(
      
      PMID %in% screening_priority$PMID ~ 1L,
      
      Auto_verdict == "" ~ 2L,
      
      Auto_verdict == "Exclude" ~ 3L,
      
      TRUE ~ 4L
    ),
    
    # Blank fields for independent manual screening
    `Reviewer 1` = NA_character_,
    `Reviewer 2` = NA_character_,
    
    # Blank fields for final adjudication
    Final_verdict = NA_character_,
    Final_gate = NA_character_,
    Final_reason = NA_character_,
    Notes = NA_character_
  ) %>%
  
# --------------------------------------------------------
# Put records requiring human attention at the top
# --------------------------------------------------------

arrange(
  manual_review_order,
  desc(suppressWarnings(as.integer(Year))),
  Reference
)%>%
  
  # Select ONLY columns required for human review
  # manual_review_order is deliberately omitted here
  select(
    DOI,
    PMID,
    Reference,
    Year,
    Title,
    Journal,
    Abstract,
    
    Auto_verdict,
    Auto_gate,
    Auto_reason,
    
    `Reviewer 1`,
    `Reviewer 2`,
    
    Final_verdict,
    Final_gate,
    Final_reason,
    Notes
  )


# ------------------------------------------------------------
# MY SANITY CHECK: MANUAL REVIEW ORDER
# ------------------------------------------------------------

n_priority_manual <- nrow(screening_priority)

message(
  "✓ High-priority records placed first in workbook: ",
  n_priority_manual
)

if (n_priority_manual != 58) {
  warning(
    "Expected 58 high-priority manual records, but found ",
    n_priority_manual,
    ". Check screening_priority criteria."
  )
}


# ------------------------------------------------------------
# 2. CREATE WORKBOOK
# ------------------------------------------------------------

wb <- createWorkbook()

addWorksheet(
  wb,
  sheetName = "Manual screening",
  gridLines = FALSE
)


# ------------------------------------------------------------
# 3. WRITE DATA
# ------------------------------------------------------------

writeData(
  wb,
  sheet = "Manual screening",
  x = manual_screening,
  startRow = 1,
  headerStyle = NULL,
  withFilter = TRUE
)

# Set all data rows to approximately 25 px
setRowHeights(
  wb,
  sheet = "Manual screening",
  rows = 2:(nrow(manual_screening) + 1),
  heights = 18.75
)


# ------------------------------------------------------------
# 4. DEFINE STYLES
# ------------------------------------------------------------

study_header <- createStyle(
  fgFill = "#D9EAF7",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

auto_header <- createStyle(
  fgFill = "#FFF2CC",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

reviewer_header <- createStyle(
  fgFill = "#E2F0D9",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

final_header <- createStyle(
  fgFill = "#F4CCCC",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

body_style <- createStyle(
  valign = "top"
)

wrap_style <- createStyle(
  wrapText = TRUE,
  valign = "top"
)


# ------------------------------------------------------------
# 5. APPLY HEADER STYLES BY SECTION
# ------------------------------------------------------------

# Study information
addStyle(
  wb,
  "Manual screening",
  style = study_header,
  rows = 1,
  cols = 1:7,
  gridExpand = TRUE
)

# Automated screening
addStyle(
  wb,
  "Manual screening",
  style = auto_header,
  rows = 1,
  cols = 8:10,
  gridExpand = TRUE
)

# Independent reviewers
addStyle(
  wb,
  "Manual screening",
  style = reviewer_header,
  rows = 1,
  cols = 11:12,
  gridExpand = TRUE
)

# Final adjudication
addStyle(
  wb,
  "Manual screening",
  style = final_header,
  rows = 1,
  cols = 13:16,
  gridExpand = TRUE
)


# ------------------------------------------------------------
# 6. FORMAT DATA ROWS
# ------------------------------------------------------------

n_rows <- nrow(manual_screening) + 1

if (n_rows > 1) {
  
  addStyle(
    wb,
    "Manual screening",
    style = body_style,
    rows = 2:n_rows,
    cols = 1:16,
    gridExpand = TRUE
  )
  
  # Wrap long text columns
  addStyle(
    wb,
    "Manual screening",
    style = wrap_style,
    rows = 2:n_rows,
    cols = c(
      5,   # Title
      7,   # Abstract
      10,  # Auto_reason
      15,  # Final_reason
      16   # Notes
    ),
    gridExpand = TRUE,
    stack = TRUE
  )
}


# ------------------------------------------------------------
# 7. SET COLUMN WIDTHS
# ------------------------------------------------------------

setColWidths(wb, "Manual screening", cols = 1, widths = 24)  # DOI
setColWidths(wb, "Manual screening", cols = 2, widths = 12)  # PMID
setColWidths(wb, "Manual screening", cols = 3, widths = 22)  # Reference
setColWidths(wb, "Manual screening", cols = 4, widths = 9)   # Year
setColWidths(wb, "Manual screening", cols = 5, widths = 55)  # Title
setColWidths(wb, "Manual screening", cols = 6, widths = 28)  # Journal
setColWidths(wb, "Manual screening", cols = 7, widths = 80)  # Abstract

setColWidths(wb, "Manual screening", cols = 8, widths = 16)  # Auto verdict
setColWidths(wb, "Manual screening", cols = 9, widths = 14)  # Auto gate
setColWidths(wb, "Manual screening", cols = 10, widths = 45) # Auto reason

setColWidths(wb, "Manual screening", cols = 11:12, widths = 18)

setColWidths(wb, "Manual screening", cols = 13, widths = 16)
setColWidths(wb, "Manual screening", cols = 14, widths = 14)
setColWidths(wb, "Manual screening", cols = 15, widths = 45)
setColWidths(wb, "Manual screening", cols = 16, widths = 40)


# ------------------------------------------------------------
# 8. FREEZE HEADER + IDENTIFYING COLUMNS
# ------------------------------------------------------------

freezePane(
  wb,
  sheet = "Manual screening",
  firstActiveRow = 2,
  firstActiveCol = 4
)


# ------------------------------------------------------------
# 9. ADD DATA VALIDATION DROPDOWNS
# ------------------------------------------------------------

if (n_rows > 1) {
  
  # Reviewer decisions
  dataValidation(
    wb,
    "Manual screening",
    cols = 11:12,
    rows = 2:n_rows,
    type = "list",
    value = '"Include,Exclude,Unsure"'
  )
  
  # Final verdict
  dataValidation(
    wb,
    "Manual screening",
    cols = 13,
    rows = 2:n_rows,
    type = "list",
    value = '"Include,Exclude"'
  )
  
  # Final gate
  dataValidation(
    wb,
    "Manual screening",
    cols = 14,
    rows = 2:n_rows,
    type = "list",
    value = '"Gate 1,Gate 2,Gate 3,Gate 4"'
  )
}


# ------------------------------------------------------------
# 10. OPTIONAL INSTRUCTIONS SHEET
# ------------------------------------------------------------

addWorksheet(
  wb,
  sheetName = "Instructions",
  gridLines = FALSE
)

instructions <- data.frame(
  Field = c(
    "Reviewer 1",
    "Reviewer 2",
    "Final_verdict",
    "Final_gate",
    "Final_reason",
    "Notes"
  ),
  
  Instructions = c(
    "Independent manual decision: Include, Exclude, or Unsure.",
    "Independent manual decision: Include, Exclude, or Unsure.",
    "Consensus decision after manual review.",
    "Eligibility gate responsible for exclusion.",
    "Concise standardised reason for the final decision.",
    "Optional free-text comments, ambiguities, or adjudication notes."
  )
)

writeData(
  wb,
  "Instructions",
  instructions,
  startRow = 1,
  withFilter = FALSE
)

setColWidths(
  wb,
  "Instructions",
  cols = 1,
  widths = 22
)

setColWidths(
  wb,
  "Instructions",
  cols = 2,
  widths = 80
)

addStyle(
  wb,
  "Instructions",
  style = study_header,
  rows = 1,
  cols = 1:2,
  gridExpand = TRUE
)


# ------------------------------------------------------------
# 11. SAVE WORKBOOK
# ------------------------------------------------------------

output_file <- "PubMed_screening_manual_review.xlsx"

saveWorkbook(
  wb,
  file = output_file,
  overwrite = TRUE
)

message(
  "✓ Manual screening workbook created: ",
  output_file
)


###############################################################
# NOTE TO SELF SCREENING HIERARCHY
#
# +---------------------------+----------------+---------------------------------------------------------------+
# | Group                     | Approx. size   | Definition                                                    |
# +---------------------------+----------------+---------------------------------------------------------------+
# | Automatically excluded    | ~500           | Positive evidence that the study fails a                     |
# |                           |                | prespecified exclusion criterion (e.g. review,               |
# |                           |                | meta-analysis, candidate-gene, non-human, etc.).             |
# +---------------------------+----------------+---------------------------------------------------------------+
# | Manual screening          | ~500–700       | No positive evidence for automatic exclusion.                |
# |                           |                | Eligibility cannot be determined from the                    |
# |                           |                | PubMed title/abstract alone.                                 |
# +---------------------------+----------------+---------------------------------------------------------------+
# | High-priority manual      | 58             | Subset of manual-review records with positive                |
# |                           |                | evidence for:                                                |
# |                           |                |   • brain tissue                                             |
# |                           |                |   • peripheral tissue                                        |
# |                           |                |   • matched/cross-tissue sampling                            |
# |                           |                |   • genome-wide DNA methylation                              |
# |                           |                | These are reviewed first.                                    |
# +---------------------------+----------------+---------------------------------------------------------------+
###############################################################



#.   1238 papers retrieved
#.   │
#.   ├── ~500 automatically excluded
#.   │
#.   └── ~700 not automatically excludable
#.   │
#.   ├── 58
#.   │    Candidate papers to screen
#.   │
#.   └── ~640 Insufficient information in title/abstract to automatically exclude


###############################################################
# CREATE INTERACTIVE SANKEY DIAGRAM OF SCREENING PIPELINE
#
# The Sankey is generated directly from the 'screening'
# dataframe. Automated exclusions use the mutually exclusive
# Auto_gate and Auto_reason variables. Records not automatically
# excluded proceed through the gates and ultimately enter manual
# review.
#
# Manual-review records are divided into:
#   1. High-priority manual records
#   2. Other records requiring manual review
#
# Because Auto_gate is assigned using case_when(), each
# automatically excluded study contributes to ONE exclusion
# reason only, preventing double-counting in the Sankey.
###############################################################

# Install once if required:
install.packages(c("networkD3", "htmlwidgets"))

library(networkD3)
library(htmlwidgets)


###############################################################
# 1. CHECK INPUT DATA
###############################################################

required_sankey_columns <- c(
  "PMID",
  "Auto_verdict",
  "Auto_gate",
  "Auto_reason"
)

missing_sankey_columns <- setdiff(
  required_sankey_columns,
  names(screening)
)

if (length(missing_sankey_columns) > 0) {
  stop(
    "Cannot create Sankey diagram. Missing columns: ",
    paste(missing_sankey_columns, collapse = ", ")
  )
}


###############################################################
# 2. COUNT AUTOMATED EXCLUSIONS BY GATE AND REASON
###############################################################

sankey_exclusions <- screening %>%
  filter(
    Auto_verdict == "Exclude",
    Auto_gate != "",
    Auto_reason != ""
  ) %>%
  count(
    Auto_gate,
    Auto_reason,
    name = "n"
  ) %>%
  arrange(
    Auto_gate,
    desc(n),
    Auto_reason
  )


###############################################################
# 3. CALCULATE RECORD FLOW THROUGH EACH GATE
###############################################################

n_unique <- nrow(screening)

n_gate1_excluded <- sum(
  screening$Auto_gate == "Gate 1" &
    screening$Auto_verdict == "Exclude",
  na.rm = TRUE
)

n_after_gate1 <- n_unique - n_gate1_excluded


n_gate2_excluded <- sum(
  screening$Auto_gate == "Gate 2" &
    screening$Auto_verdict == "Exclude",
  na.rm = TRUE
)

n_after_gate2 <- n_after_gate1 - n_gate2_excluded


n_gate3_excluded <- sum(
  screening$Auto_gate == "Gate 3" &
    screening$Auto_verdict == "Exclude",
  na.rm = TRUE
)

n_after_gate3 <- n_after_gate2 - n_gate3_excluded


n_gate4_excluded <- sum(
  screening$Auto_gate == "Gate 4" &
    screening$Auto_verdict == "Exclude",
  na.rm = TRUE
)

n_after_gate4 <- n_after_gate3 - n_gate4_excluded


###############################################################
# 4. MANUAL-REVIEW COUNTS
###############################################################

n_manual <- sum(
  screening$Auto_verdict == "",
  na.rm = TRUE
)

n_priority <- nrow(screening_priority)

n_other_manual <- n_manual - n_priority


###############################################################
# 5. INTEGRITY CHECKS
###############################################################

if (n_other_manual < 0) {
  stop(
    "Sankey count error: high-priority manual records exceed ",
    "the total number of manual-review records."
  )
}

if (n_after_gate4 != n_manual) {
  warning(
    "Sankey flow check: ",
    n_after_gate4,
    " records remain after automated gates, but ",
    n_manual,
    " records are classified for manual review."
  )
}

n_total_auto_excluded <- sum(
  screening$Auto_verdict == "Exclude",
  na.rm = TRUE
)

if ((n_total_auto_excluded + n_manual) != n_unique) {
  warning(
    "Screening counts do not reconcile: ",
    n_total_auto_excluded,
    " automatically excluded + ",
    n_manual,
    " manual-review records != ",
    n_unique,
    " unique records."
  )
}


###############################################################
# 6. CREATE SANKEY NODES
###############################################################

# Core flow nodes

core_nodes <- tibble(
  name = c(
    "Unique PubMed records",
    "Gate 1 passed",
    "Gate 2 passed",
    "Gate 3 passed",
    "Manual review",
    "High-priority manual",
    "Other manual review"
  )
)


# Exclusion-reason nodes
#
# Prefix the reason with the gate so that identical reason names
# from different gates could never generate duplicate node names.

reason_nodes <- sankey_exclusions %>%
  mutate(
    name = paste0(
      Auto_gate,
      ": ",
      Auto_reason
    )
  ) %>%
  distinct(name) %>%
  select(name)


# Combine all nodes

sankey_nodes <- bind_rows(
  core_nodes,
  reason_nodes
)


###############################################################
# 7. HELPER FUNCTION TO FIND NODE INDEX
#
# networkD3 uses ZERO-BASED indexing, so subtract 1 from the
# ordinary R row number.
###############################################################

node_id <- function(node_name) {
  
  idx <- match(
    node_name,
    sankey_nodes$name
  )
  
  if (is.na(idx)) {
    stop(
      "Sankey node not found: ",
      node_name
    )
  }
  
  idx - 1
}


###############################################################
# 8. CREATE CORE LINKS THROUGH THE SCREENING GATES
###############################################################

sankey_links <- tibble(
  source_name = character(),
  target_name = character(),
  value = numeric()
)


###############################################################
# Gate 1:
# Unique records -> Gate 1 exclusion reasons / Gate 1 passed
###############################################################

gate1_reasons <- sankey_exclusions %>%
  filter(Auto_gate == "Gate 1")

if (nrow(gate1_reasons) > 0) {
  
  sankey_links <- bind_rows(
    sankey_links,
    
    gate1_reasons %>%
      transmute(
        source_name = "Unique PubMed records",
        target_name = paste0(
          "Gate 1: ",
          Auto_reason
        ),
        value = n
      )
  )
}

if (n_after_gate1 > 0) {
  
  sankey_links <- bind_rows(
    sankey_links,
    tibble(
      source_name = "Unique PubMed records",
      target_name = "Gate 1 passed",
      value = n_after_gate1
    )
  )
}


###############################################################
# Gate 2:
# Gate 1 passed -> Gate 2 exclusion reasons / Gate 2 passed
###############################################################

gate2_reasons <- sankey_exclusions %>%
  filter(Auto_gate == "Gate 2")

if (nrow(gate2_reasons) > 0) {
  
  sankey_links <- bind_rows(
    sankey_links,
    
    gate2_reasons %>%
      transmute(
        source_name = "Gate 1 passed",
        target_name = paste0(
          "Gate 2: ",
          Auto_reason
        ),
        value = n
      )
  )
}

if (n_after_gate2 > 0) {
  
  sankey_links <- bind_rows(
    sankey_links,
    tibble(
      source_name = "Gate 1 passed",
      target_name = "Gate 2 passed",
      value = n_after_gate2
    )
  )
}


###############################################################
# Gate 3:
# Gate 2 passed -> Gate 3 exclusion reasons / Gate 3 passed
###############################################################

gate3_reasons <- sankey_exclusions %>%
  filter(Auto_gate == "Gate 3")

if (nrow(gate3_reasons) > 0) {
  
  sankey_links <- bind_rows(
    sankey_links,
    
    gate3_reasons %>%
      transmute(
        source_name = "Gate 2 passed",
        target_name = paste0(
          "Gate 3: ",
          Auto_reason
        ),
        value = n
      )
  )
}

if (n_after_gate3 > 0) {
  
  sankey_links <- bind_rows(
    sankey_links,
    tibble(
      source_name = "Gate 2 passed",
      target_name = "Gate 3 passed",
      value = n_after_gate3
    )
  )
}


###############################################################
# Gate 4
#
# At present your automated pipeline will normally make no
# Gate 4 exclusions because data availability requires manual
# full-text assessment.
#
# The code below nevertheless supports automated Gate 4
# exclusions if these are added to the pipeline in future.
###############################################################

gate4_reasons <- sankey_exclusions %>%
  filter(Auto_gate == "Gate 4")

if (nrow(gate4_reasons) > 0) {
  
  sankey_links <- bind_rows(
    sankey_links,
    
    gate4_reasons %>%
      transmute(
        source_name = "Gate 3 passed",
        target_name = paste0(
          "Gate 4: ",
          Auto_reason
        ),
        value = n
      )
  )
  
  if (n_after_gate4 > 0) {
    
    sankey_links <- bind_rows(
      sankey_links,
      tibble(
        source_name = "Gate 3 passed",
        target_name = "Manual review",
        value = n_after_gate4
      )
    )
  }
  
} else {
  
  # No automated Gate 4 screening at present:
  # all Gate 3 survivors proceed to manual review.
  
  if (n_after_gate3 > 0) {
    
    sankey_links <- bind_rows(
      sankey_links,
      tibble(
        source_name = "Gate 3 passed",
        target_name = "Manual review",
        value = n_after_gate3
      )
    )
  }
}


###############################################################
# Manual review:
# split into high-priority and other manual records
###############################################################

if (n_priority > 0) {
  
  sankey_links <- bind_rows(
    sankey_links,
    tibble(
      source_name = "Manual review",
      target_name = "High-priority manual",
      value = n_priority
    )
  )
}


if (n_other_manual > 0) {
  
  sankey_links <- bind_rows(
    sankey_links,
    tibble(
      source_name = "Manual review",
      target_name = "Other manual review",
      value = n_other_manual
    )
  )
}


###############################################################
# 9. CONVERT NODE NAMES TO networkD3 INDICES
###############################################################

sankey_links <- sankey_links %>%
  mutate(
    source = vapply(
      source_name,
      node_id,
      numeric(1)
    ),
    
    target = vapply(
      target_name,
      node_id,
      numeric(1)
    )
  ) %>%
  select(
    source,
    target,
    value,
    source_name,
    target_name
  )


###############################################################
# 10. REMOVE UNUSED NODES
#
# This prevents nodes for gates/reasons that currently have
# zero observations from appearing in the diagram.
###############################################################

used_node_names <- unique(
  c(
    sankey_links$source_name,
    sankey_links$target_name
  )
)

sankey_nodes <- sankey_nodes %>%
  filter(name %in% used_node_names)


###############################################################
# Recalculate indices because unused nodes were removed
###############################################################

sankey_links <- sankey_links %>%
  mutate(
    source = match(
      source_name,
      sankey_nodes$name
    ) - 1,
    
    target = match(
      target_name,
      sankey_nodes$name
    ) - 1
  )


###############################################################
# 11. ADD NODE GROUPS FOR VISUAL ORGANISATION
###############################################################

sankey_nodes <- sankey_nodes %>%
  mutate(
    group = case_when(
      
      name == "Unique PubMed records" ~
        "Search",
      
      str_detect(name, "^Gate 1") ~
        "Gate 1",
      
      str_detect(name, "^Gate 2") ~
        "Gate 2",
      
      str_detect(name, "^Gate 3") ~
        "Gate 3",
      
      str_detect(name, "^Gate 4") ~
        "Gate 4",
      
      name %in% c(
        "Manual review",
        "High-priority manual",
        "Other manual review"
      ) ~
        "Manual review",
      
      TRUE ~
        "Other"
    )
  )


###############################################################
# 12. ADD LINK GROUPS
###############################################################

sankey_links <- sankey_links %>%
  mutate(
    group = case_when(
      
      str_detect(target_name, "^Gate 1:") ~
        "Excluded Gate 1",
      
      str_detect(target_name, "^Gate 2:") ~
        "Excluded Gate 2",
      
      str_detect(target_name, "^Gate 3:") ~
        "Excluded Gate 3",
      
      str_detect(target_name, "^Gate 4:") ~
        "Excluded Gate 4",
      
      target_name == "High-priority manual" ~
        "High priority",
      
      target_name == "Other manual review" ~
        "Manual review",
      
      TRUE ~
        "Continued"
    )
  )


###############################################################
# 13. CREATE INTERACTIVE SANKEY
###############################################################

screening_sankey <- sankeyNetwork(
  Links = sankey_links,
  Nodes = sankey_nodes,
  
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  
  NodeGroup = "group",
  LinkGroup = "group",
  
  fontSize = 12,
  nodeWidth = 30,
  nodePadding = 15,
  
  sinksRight = FALSE
)


###############################################################
# 14. ADD TITLE ABOVE THE SANKEY
###############################################################

screening_sankey <- htmlwidgets::prependContent(
  screening_sankey,
  
  htmltools::tags$div(
    
    style = "
      font-family: Arial, sans-serif;
      margin-bottom: 15px;
    ",
    
    htmltools::tags$h2(
      "PubMed screening flow"
    ),
    
    htmltools::tags$p(
      paste0(
        format(n_unique, big.mark = ","),
        " unique records identified; ",
        format(n_total_auto_excluded, big.mark = ","),
        " automatically excluded; ",
        format(n_manual, big.mark = ","),
        " retained for manual screening."
      )
    )
  )
)


###############################################################
# 15. SAVE INTERACTIVE HTML FILE
###############################################################

sankey_output_file <- "PubMed_screening_sankey.html"

saveWidget(
  screening_sankey,
  file = sankey_output_file,
  selfcontained = TRUE
)

message(
  "✓ Screening Sankey diagram created: ",
  sankey_output_file
)


###############################################################
# 16. PRINT SANKEY SUMMARY TO CONSOLE
###############################################################

cat(
  "\nSANKEY SCREENING FLOW\n",
  "---------------------\n",
  
  "Unique records:                ",
  format(n_unique, big.mark = ","), "\n",
  
  "Gate 1 exclusions:             ",
  format(n_gate1_excluded, big.mark = ","), "\n",
  
  "Remaining after Gate 1:        ",
  format(n_after_gate1, big.mark = ","), "\n",
  
  "Gate 2 exclusions:             ",
  format(n_gate2_excluded, big.mark = ","), "\n",
  
  "Remaining after Gate 2:        ",
  format(n_after_gate2, big.mark = ","), "\n",
  
  "Gate 3 exclusions:             ",
  format(n_gate3_excluded, big.mark = ","), "\n",
  
  "Remaining after Gate 3:        ",
  format(n_after_gate3, big.mark = ","), "\n",
  
  "Gate 4 automated exclusions:   ",
  format(n_gate4_excluded, big.mark = ","), "\n",
  
  "Manual review:                 ",
  format(n_manual, big.mark = ","), "\n",
  
  "  High-priority manual:        ",
  format(n_priority, big.mark = ","), "\n",
  
  "  Other manual review:         ",
  format(n_other_manual, big.mark = ","), "\n\n",
  
  sep = ""
)





###############################################################
# Final pipeline summary
###############################################################

end_time <- Sys.time()

runtime <- as.numeric(
  difftime(end_time, start_time, units = "secs")
)

hours <- runtime %/% 3600
minutes <- (runtime %% 3600) %/% 60
seconds <- round(runtime %% 60)

###############################################################
# Pull automated failure counts by gate
###############################################################

gate1_auto <- sum(screening$Auto_gate == "Gate 1", na.rm = TRUE)
gate2_auto <- sum(screening$Auto_gate == "Gate 2", na.rm = TRUE)
gate3_auto <- sum(screening$Auto_gate == "Gate 3", na.rm = TRUE)
gate4_auto <- sum(screening$Auto_gate == "Gate 4", na.rm = TRUE)

###############################################################
# Detailed automated exclusion counts
###############################################################

count_auto_reason <- function(reason_text) {
  sum(
    screening$Auto_reason == reason_text &
      screening$Auto_verdict == "Exclude",
    na.rm = TRUE
  )
}

n_review <- count_auto_reason("Review article")
n_meta <- count_auto_reason("Meta-analysis")
n_candidate <- count_auto_reason("Candidate-gene study")
n_secondary <- count_auto_reason("Secondary analysis of previously published data")
n_abstract <- count_auto_reason("Abstract only")
n_nonhuman <- count_auto_reason("Non-human tissue/sample")
n_differential <- count_auto_reason(
  "Differential methylation only / DMP-DMR analysis"
)

n_neuroimaging_no_brain_dnam <- count_auto_reason(
  "Neuroimaging phenotype; no brain-tissue DNAm"
)

n_unmatched_tissues <- count_auto_reason(
  "Brain and peripheral samples explicitly unmatched"
)

cat(
  "\n=========================================================\n",
  "PubMed screening pipeline completed successfully\n",
  "=========================================================\n\n",

  "SEARCH RESULTS\n",
  "--------------\n",
  "Records retrieved:              ",
  format(search_results$count, big.mark = ","), "\n",
  "Unique PubMed records:          ",
  format(nrow(screening), big.mark = ","), "\n",
  "Duplicates removed:             ",
  format(n_duplicates_removed, big.mark = ","), "\n\n",

  "AUTOMATED GATE FAILURES\n",
  "-----------------------\n",
  "Gate 1 - Study design:          ",
  format(gate1_auto, big.mark = ","), "\n",
  "Gate 2 - Tissue sampling:       ",
  format(gate2_auto, big.mark = ","),
  "  [additional manual assessment may be required]\n",
  "Gate 3 - Analysis:              ",
  format(gate3_auto, big.mark = ","), "\n",
  "Gate 4 - Data availability:     ",
  format(gate4_auto, big.mark = ","),
  "  [manual assessment normally required]\n\n",

  "AUTOMATED EXCLUSIONS BY REASON\n",
  "------------------------------\n",
  "Gate 1 - Study design\n",
  "  Review articles:              ",
  format(n_review, big.mark = ","), "\n",
  "  Meta-analyses:                ",
  format(n_meta, big.mark = ","), "\n",
  "  Candidate-gene studies:       ",
  format(n_candidate, big.mark = ","), "\n",
  "  Secondary analyses:           ",
  format(n_secondary, big.mark = ","), "\n",
  "  Abstract-only reports:        ",
  format(n_abstract, big.mark = ","), "\n",
  "  Non-human samples:            ",
  format(n_nonhuman, big.mark = ","), "\n",
  "  Gate 1 total:                 ",
  format(gate1_auto, big.mark = ","), "\n\n",

  "Gate 2 - Tissue sampling\n",
  "  Neuroimaging-epigenomics study; no brain tissue DNAm sampled:  ",
  format(n_neuroimaging_no_brain_dnam, big.mark = ","), "\n",
  "  Explicitly unmatched samples: ",
  format(n_unmatched_tissues, big.mark = ","), "\n",
  "  Gate 2 total:                 ",
  format(gate2_auto, big.mark = ","), "\n\n",

  "Gate 3 - Analysis\n",
  "  Differential methylation:     ",
  format(n_differential, big.mark = ","), "\n",
  "  Gate 3 total:                 ",
  format(gate3_auto, big.mark = ","), "\n\n",

  "SCREENING STATUS\n",
  "----------------\n",
  "Automatically excluded:         ",
  format(nrow(screening_auto_excluded), big.mark = ","), "\n",
  
  "Potential for review according to Gate: ",
  format(nrow(screening_manual_review), big.mark = ","), "\n",
  
  "Exclusions for EC or EL to check: ",
  format(n_unmatched_tissues, big.mark = ","), "\n",

  "High-priority manual records:   ",
  format(nrow(screening_priority), big.mark = ","), "\n\n",


  
  
  
  
  "OUTPUT FILES\n",
  "------------\n",
  "✓ PubMed_screening_records.csv\n",
  "✓ PubMed_screening_manual.csv\n",
  "✓ PubMed_screening_automatic_exclusions.csv\n",
  "✓ PubMed_screening_gate_summary.csv\n",
  "✓ PubMed_screening_exclusion_dictionary.csv\n",
  "✓ PubMed_screening_priority.csv\n\n",

  "Total runtime: ",
  sprintf("%02d:%02d:%02d", hours, minutes, seconds), "\n",
  "Finished: ",
  format(end_time, "%Y-%m-%d %H:%M:%S"), "\n",
  "=========================================================\n",
  sep = ""
)

###############################################################
# NOTE FOR FINAL MANUAL SCREENING
#
# Once manual/full-text screening is complete, use Final_verdict,
# Final_gate and Final_reason as the authoritative fields.
#
# Final gate counts can then be calculated with:
#
# final_gate_summary <- completed_screening %>%
#   filter(Final_verdict == "Exclude") %>%
#   count(Final_gate, Final_reason, name = "n_excluded") %>%
#   arrange(Final_gate, desc(n_excluded))
#
# This ensures that each excluded study contributes to one
# prespecified gate and one final exclusion reason.
###############################################################
