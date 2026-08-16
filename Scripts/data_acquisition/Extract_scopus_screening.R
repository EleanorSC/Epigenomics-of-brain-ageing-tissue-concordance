###############################################################
# Extract_scopus_screening.R
#
# Purpose:
# Run the preregistered Scopus search directly through the
# official Elsevier Scopus Search API, retrieve citation
# metadata and abstracts, apply conservative automated
# screening flags, organise exclusions by prespecified
# screening gate, and create CSV/workbook/Sankey outputs.
#
# REQUIREMENT:
# To do this, you will require an Elsevier API key with Scopus access.
# You can get this via institutional login via the Elsevier Developer Portal:
# https://dev.elsevier.com/apikey/
# Before running this script, set your key once per R session as follows:
#
# Sys.setenv(ELSEVIER_API_KEY = "YOUR INDIVIDUAL API KEY")
# Sys.getenv("ELSEVIER_API_KEY")
# nzchar(Sys.getenv("ELSEVIER_API_KEY"))
#
# Do NOT hard-code the API key into your script if it will later be
# committed to GitHub as I have done with this script
###############################################################

# install.packages(c(
#   "httr2", "jsonlite", "dplyr", "stringr", "purrr",
#   "tibble", "openxlsx", "networkD3", "htmlwidgets",
#   "htmltools", "wesanderson"
# ))
#
###############################################################
# Start timer
###############################################################

start_time <- Sys.time()

cat(
  "\n=========================================================\n",
  "Brain–Peripheral DNA Methylation Concordance Project\n",
  "Gate-based Scopus Screening Pipeline\n",
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

library(httr2)
library(jsonlite)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

###############################################################
# 1. Define Scopus search
#
# TITLE-ABS() searches title and abstract in Scopus.
# This is the direct Scopus translation of the preregistered
# PubMed [Title/Abstract] search.
###############################################################

scopus_query <- paste0(
  "TITLE-ABS(",
  "\"DNA methylation\" OR DNAm OR methylome OR EWAS OR \"epigenome-wide\"",
  ") AND TITLE-ABS(",
  "blood OR serum OR placenta OR plasma OR saliva OR buccal OR cheek",
  ") AND TITLE-ABS(",
  "brain OR \"brain tissue\" OR cortex OR \"cerebral cortex\" OR ",
  "\"prefrontal cortex\" OR \"dorsolateral prefrontal cortex\" OR DLPFC OR ",
  "\"entorhinal cortex\" OR \"superior temporal gyrus\" OR STG OR ",
  "\"inferior frontal gyrus\" OR cerebellum OR cerebellar OR ",
  "hippocampus OR hippocampal OR amygdala OR \"nucleus accumbens\" OR ",
  "striatum OR putamen OR caudate OR \"substantia nigra\"",
  ")"
)

cat(
  "\nSCOPUS QUERY\n",
  "============\n",
  scopus_query,
  "\n\n",
  sep = ""
)

###############################################################
# 2. Check Elsevier API key
###############################################################

api_key <- ELSEVIER_API_KEY

if (api_key == "") {
  stop(
    "\nNo Elsevier API key found.\n\n",
    "Set it in this R session with:\n",
    "Sys.setenv(ELSEVIER_API_KEY = \"ELSEVIER_API_KEY\")\n\n",
    "Then rerun the script."
  )
}

###############################################################
# 3. Search Scopus and retrieve STANDARD records in batches
#
# STANDARD view contains the abstract (dc:description) and
# richer author/affiliation metadata.
###############################################################

t <- stage_start(1, "Searching and downloading Scopus records")

scopus_endpoint <- "https://api.elsevier.com/content/search/scopus"

batch_size <- 25L
max_retries <- 3L

make_scopus_request <- function(start = 0L, count = batch_size) {

  request(scopus_endpoint) %>%
    req_headers(
      `X-ELS-APIKey` = api_key,
      Accept = "application/json"
    ) %>%
    req_url_query(
      query = scopus_query,
      view = "STANDARD",
      start = start,
      count = count,
      suppressNavLinks = "true"
    ) %>%
    req_retry(
      max_tries = max_retries,
      backoff = function(tries) 2^tries
    ) %>%
    req_error(
      body = function(resp) {
        paste0(
          "Scopus API request failed (HTTP ",
          resp_status(resp),
          "): ",
          resp_body_string(resp)
        )
      }
    ) %>%
    req_perform()
}

# First request determines total result count
first_response <- make_scopus_request(start = 0L, count = batch_size)
first_json <- resp_body_json(first_response, simplifyVector = FALSE)

n_scopus_retrieved <- as.integer(
  first_json[["search-results"]][["opensearch:totalResults"]]
)

if (is.na(n_scopus_retrieved)) {
  stop("Could not determine the total number of Scopus search results.")
}

cat(
  "  Scopus search returned ",
  format(n_scopus_retrieved, big.mark = ","),
  " records.\n",
  sep = ""
)

extract_entries <- function(x) {
  entries <- x[["search-results"]][["entry"]]

  if (is.null(entries)) {
    return(list())
  }

  entries
}

scopus_entries <- extract_entries(first_json)

if (n_scopus_retrieved > batch_size) {

  batch_starts <- seq(
    from = batch_size,
    to = n_scopus_retrieved - 1L,
    by = batch_size
  )

  for (i in seq_along(batch_starts)) {

    start_i <- batch_starts[i]

    cat(
      "  Downloading batch ",
      i + 1L,
      "/",
      length(batch_starts) + 1L,
      " (records ",
      start_i + 1L,
      "-",
      min(start_i + batch_size, n_scopus_retrieved),
      ")...\n",
      sep = ""
    )

    response_i <- make_scopus_request(
      start = start_i,
      count = min(batch_size, n_scopus_retrieved - start_i)
    )

    json_i <- resp_body_json(
      response_i,
      simplifyVector = FALSE
    )

    scopus_entries <- c(
      scopus_entries,
      extract_entries(json_i)
    )

    Sys.sleep(0.25)
  }
}

stage_end(
  t,
  paste0(
    format(length(scopus_entries), big.mark = ","),
    " Scopus records downloaded"
  )
)

###############################################################
# 4. Retrieval integrity check
###############################################################

t <- stage_start(2, "Checking Scopus retrieval integrity")

if (length(scopus_entries) != n_scopus_retrieved) {

  warning(
    "Scopus retrieval count mismatch: search returned ",
    n_scopus_retrieved,
    " records but ",
    length(scopus_entries),
    " records were downloaded."
  )

} else {

  cat(
    "✓ Retrieval integrity check passed: ",
    format(length(scopus_entries), big.mark = ","),
    " / ",
    format(n_scopus_retrieved, big.mark = ","),
    " records retrieved\n",
    sep = ""
  )
}

stage_end(t, "Scopus retrieval checked")

###############################################################
# 5. Helper functions for Scopus JSON
###############################################################

scalar_chr <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_character_)
  }

  if (is.list(x) && !is.null(x[["$"]])) {
    x <- x[["$"]]
  }

  x <- unlist(x, use.names = FALSE)

  if (length(x) == 0) {
    return(NA_character_)
  }

  x <- as.character(x[1])

  if (is.na(x) || x == "") {
    return(NA_character_)
  }

  x
}

get_scopus_authors <- function(entry) {

  authors <- entry[["author"]]

  if (is.null(authors) || length(authors) == 0) {
    creator <- scalar_chr(entry[["dc:creator"]])
    return(creator)
  }

  # A single author may arrive as a named list rather than
  # a list of author records.
  if (!is.null(authors[["authname"]])) {
    authors <- list(authors)
  }

  author_names <- vapply(
    authors,
    function(a) {
      name <- scalar_chr(a[["authname"]])

      if (!is.na(name)) {
        return(name)
      }

      surname <- scalar_chr(a[["surname"]])
      initials <- scalar_chr(a[["initials"]])

      paste(
        na.omit(c(surname, initials)),
        collapse = " "
      )
    },
    character(1)
  )

  author_names <- author_names[
    !is.na(author_names) &
      author_names != ""
  ]

  if (length(author_names) == 0) {
    return(NA_character_)
  }

  paste(author_names, collapse = "; ")
}

get_first_scopus_author <- function(entry) {

  authors <- entry[["author"]]

  if (!is.null(authors) && length(authors) > 0) {

    if (!is.null(authors[["authname"]])) {
      first <- authors
    } else {
      first <- authors[[1]]
    }

    name <- scalar_chr(first[["authname"]])

    if (!is.na(name)) {
      return(name)
    }

    surname <- scalar_chr(first[["surname"]])
    initials <- scalar_chr(first[["initials"]])

    out <- paste(
      na.omit(c(surname, initials)),
      collapse = " "
    )

    if (out != "") {
      return(out)
    }
  }

  scalar_chr(entry[["dc:creator"]])
}

###############################################################
# 6. Extract citation metadata and abstracts
###############################################################

t <- stage_start(3, "Extracting Scopus citation metadata")

scopus_records <- map_dfr(
  scopus_entries,
  function(entry) {

    cover_date <- scalar_chr(
      entry[["prism:coverDate"]]
    )

    tibble(
      EID = scalar_chr(
        entry[["eid"]]
      ),

      Scopus_ID = str_remove(
        scalar_chr(entry[["dc:identifier"]]),
        "^SCOPUS_ID:"
      ),

      PMID = scalar_chr(
        entry[["pubmed-id"]]
      ),

      DOI = scalar_chr(
        entry[["prism:doi"]]
      ),

      First_author = get_first_scopus_author(entry),

      Authors = get_scopus_authors(entry),

      Year = if_else(
        !is.na(cover_date),
        str_sub(cover_date, 1, 4),
        NA_character_
      ),

      Journal = scalar_chr(
        entry[["prism:publicationName"]]
      ),

      Title = scalar_chr(
        entry[["dc:title"]]
      ),

      Abstract = scalar_chr(
        entry[["dc:description"]]
      ),

      Publication_types = coalesce(
        scalar_chr(entry[["subtypeDescription"]]),
        scalar_chr(entry[["subtype"]]),
        scalar_chr(entry[["prism:aggregationType"]])
      )
    )
  }
)

# Normalise DOI strings
scopus_records <- scopus_records %>%
  mutate(
    DOI = str_remove(
      str_trim(DOI),
      regex("^https?://(dx\\.)?doi\\.org/", ignore_case = TRUE)
    )
  )

###############################################################
# Metadata completeness report
###############################################################

cat(
  "\nMetadata completeness:\n",
  "  Missing abstracts: ",
  sum(is.na(scopus_records$Abstract) | scopus_records$Abstract == ""),
  "\n",
  "  Missing DOI:       ",
  sum(is.na(scopus_records$DOI) | scopus_records$DOI == ""),
  "\n",
  "  Missing PMID:      ",
  sum(is.na(scopus_records$PMID) | scopus_records$PMID == ""),
  "\n",
  sep = ""
)

stage_end(
  t,
  paste0(
    format(nrow(scopus_records), big.mark = ","),
    " records extracted"
  )
)

###############################################################
# 7. search/prioritisation fields insofar as I can automate this
###############################################################

t <- stage_start(5, "Creating screening and prioritisation fields")

scopus_records <- scopus_records %>%
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
    format(nrow(scopus_records), big.mark = ","),
    " records"
  )
)

###############################################################
# 8. Deduplicate Scopus records
#
# Scopus may contain records with DOI and/or PMID as well as its
# own Electronic Identifier (EID). Deduplication therefore uses:
# DOI -> PMID -> EID -> normalised title/year.
###############################################################

t <- stage_start(6, "Deduplicating Scopus records")

n_before_dedup <- nrow(scopus_records)

screening <- scopus_records %>%
  mutate(
    dedup_key = case_when(
      !is.na(DOI) & DOI != "" ~
        paste0("DOI:", str_to_lower(str_trim(DOI))),
      !is.na(PMID) & PMID != "" ~
        paste0("PMID:", str_trim(PMID)),
      !is.na(EID) & EID != "" ~
        paste0("EID:", EID),
      TRUE ~ paste0(
        "TITLEYEAR:",
        str_squish(str_to_lower(coalesce(Title, ""))),
        ":",
        coalesce(Year, "")
      )
    )
  ) %>%
  distinct(dedup_key, .keep_all = TRUE) %>%
  select(-dedup_key)

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
    # Gate 4: Appropriate data availability?
    #
    # Gate 4 cannot normally be determined reliably from a
    # Scopus title/abstract alone. These variables therefore
    # identify POSITIVE EVIDENCE that extractable CpG-level
    # results may be available and are used for prioritisation
    # during manual/full-text screening.
    #
    # IMPORTANT:
    # These variables do NOT automatically exclude studies.
    ###########################################################
    
    mentions_supplementary_data = str_detect(
      search_text,
      regex(
        paste(
          c(
            "supplementary data",
            "supplementary dataset",
            "supplementary datasets",
            "supplementary table",
            "supplementary tables",
            "supplemental data",
            "supplemental dataset",
            "supplemental table",
            "supplemental tables",
            "supporting information",
            "supporting data"
          ),
          collapse = "|"
        ),
        ignore_case = TRUE
      )
    ),
    
    mentions_data_repository = str_detect(
      search_text,
      regex(
        paste(
          c(
            "gene expression omnibus",
            "\\bGEO\\b",
            "\\bGSE[0-9]+\\b",
            "arrayexpress",
            "figshare",
            "zenodo",
            "dryad",
            "data repository",
            "public repository",
            "publicly available dataset",
            "publicly available data"
          ),
          collapse = "|"
        ),
        ignore_case = TRUE
      )
    ),
    
    mentions_data_availability = str_detect(
      search_text,
      regex(
        paste(
          c(
            "data are available",
            "data is available",
            "data available",
            "publicly available",
            "available online",
            "available in the supplementary",
            "available in supplementary",
            "downloadable"
          ),
          collapse = "|"
        ),
        ignore_case = TRUE
      )
    ),
    
    mentions_cpg_level_results = str_detect(
      search_text,
      regex(
        paste(
          c(
            "CpG[- ]level",
            "individual CpG",
            "individual CpGs",
            "site[- ]specific correlation",
            "site[- ]specific correlations",
            "probe[- ]level",
            "genome[- ]wide correlation",
            "genome[- ]wide correlations",
            "correlation coefficients",
            "correlation coefficient"
          ),
          collapse = "|"
        ),
        ignore_case = TRUE
      )
    ),
    
    mentions_raw_data = str_detect(
      search_text,
      regex(
        paste(
          c(
            "raw data",
            "raw methylation data",
            "raw microarray data",
            "raw intensity data",
            "raw intensity files",
            "raw IDAT",
            "raw IDATs",
            "IDAT files"
          ),
          collapse = "|"
        ),
        ignore_case = TRUE
      )
    ),
    
    ###########################################################
    # Combined Gate 4 evidence indicators
    ###########################################################
    
    gate4_data_evidence =
      mentions_supplementary_data |
      mentions_data_repository |
      mentions_data_availability |
      mentions_cpg_level_results,
    
    # Strongest abstract-level evidence that the paper may
    # provide directly extractable CpG-level results.
    gate4_high_priority =
      mentions_cpg_level_results &
      (
        mentions_supplementary_data |
          mentions_data_availability
      ),
    
    ###########################################################
    # Gate 4 triage classification
    #
    # This classification is descriptive only. Definitive
    # Gate 4 eligibility requires manual inspection of the
    # publication and its supplementary files.
    ###########################################################
    
    Gate4_evidence = case_when(
      
      gate4_high_priority ~
        "Strong evidence of downloadable CpG-level results",
      
      mentions_cpg_level_results &
        mentions_data_repository ~
        "CpG-level results and data repository mentioned",
      
      mentions_raw_data &
        !mentions_cpg_level_results &
        !mentions_supplementary_data ~
        "Raw data mentioned; check for processed concordance results",
      
      mentions_supplementary_data ~
        "Supplementary results mentioned; inspect manually",
      
      mentions_data_repository |
        mentions_data_availability ~
        "Data availability mentioned; inspect manually",
      
      TRUE ~
        "No Gate 4 information in title/abstract"
    ),
    
    ###########################################################
    # Automatic gate assignment

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
    EID,
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
    EID,
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
    EID,
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
    EID,
    PMID,
    DOI,
    Study,
    Title,
    Abstract,
    Auto_gate,
    Auto_reason
  )


###############################################################
# Sanity check: Gate 4 evidence from title/abstract
###############################################################

gate4_evidence_summary <- screening %>%
  count(
    Gate4_evidence,
    name = "n"
  ) %>%
  arrange(desc(n))

gate4_evidence_summary


###############################################################
# Review strongest Gate 4 candidates
###############################################################

gate4_strong_candidates <- screening %>%
  filter(
    gate4_high_priority == TRUE
  ) %>%
  arrange(
    desc(suppressWarnings(as.integer(Year)))
  ) %>%
  select(
    EID,
    PMID,
    DOI,
    Study,
    Year,
    Title,
    Abstract,
    Gate4_evidence
  )

gate4_strong_candidates


###############################################################
# Review studies mentioning raw data
###############################################################

gate4_raw_data_candidates <- screening %>%
  filter(
    mentions_raw_data == TRUE
  ) %>%
  arrange(
    desc(suppressWarnings(as.integer(Year)))
  ) %>%
  select(
    EID,
    PMID,
    DOI,
    Study,
    Year,
    Title,
    Abstract,
    Gate4_evidence
  )

gate4_raw_data_candidates

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
#
# These records contain positive evidence for the main
# eligibility concepts. Within this group, papers with stronger
# evidence of downloadable CpG-level data are reviewed first.
###############################################################

screening_priority <- screening %>%
  filter(
    Auto_verdict == "",
    mentions_brain == TRUE,
    mentions_peripheral == TRUE,
    mentions_matching == TRUE,
    mentions_genomewide == TRUE
  ) %>%
  mutate(
    
    Gate4_priority_order = case_when(
      
      gate4_high_priority ~ 1L,
      
      mentions_cpg_level_results &
        mentions_data_repository ~ 2L,
      
      gate4_data_evidence ~ 3L,
      
      TRUE ~ 4L
    )
  ) %>%
  arrange(
    Gate4_priority_order,
    desc(suppressWarnings(as.integer(Year)))
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
    "Unique Scopus records",
    "Duplicates removed",
    "Automatically excluded",
    "Remaining for manual screening",
    "High-priority manual records"
  ),

  n = c(
    n_scopus_retrieved,
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
  "Scopus_screening_records.csv",
  row.names = FALSE
)

write.csv(
  screening_manual_review,
  "Scopus_screening_manual.csv",
  row.names = FALSE
)

write.csv(
  screening_auto_excluded,
  "Scopus_screening_automatic_exclusions.csv",
  row.names = FALSE
)

write.csv(
  auto_gate_summary,
  "Scopus_screening_gate_summary.csv",
  row.names = FALSE
)

write.csv(
  exclusion_dictionary,
  "Scopus_screening_exclusion_dictionary.csv",
  row.names = FALSE
)

write.csv(
  screening_priority,
  "Scopus_screening_priority.csv",
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
    
    Reference = case_when(
      !is.na(First_author) & !is.na(Year) ~
        paste0(First_author, " et al. ", Year),
      !is.na(First_author) ~ First_author,
      TRUE ~ NA_character_
    ),
    
    manual_review_order = case_when(
      EID %in% screening_priority$EID ~ 1L,
      Auto_verdict == "" ~ 2L,
      Auto_verdict == "Exclude" ~ 3L,
      TRUE ~ 4L
    ),
    
    `Reviewer 1` = NA_character_,
    `Reviewer 2` = NA_character_,
    
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
    EID,
    PMID,
    Reference,
    Year,
    Title,
    Journal,
    Abstract,
    
    Auto_verdict,
    Auto_gate,
    Auto_reason,
    Gate4_evidence,
    
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

if (n_priority_manual == 0) {
  warning(
    "No high-priority manual records were identified. ",
    "Check the Scopus search and screening_priority criteria."
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
  cols = 1:8,
  gridExpand = TRUE
)

# Automated screening + Gate 4 evidence
addStyle(
  wb,
  "Manual screening",
  style = auto_header,
  rows = 1,
  cols = 9:12,
  gridExpand = TRUE
)

# Independent reviewers
addStyle(
  wb,
  "Manual screening",
  style = reviewer_header,
  rows = 1,
  cols = 13:14,
  gridExpand = TRUE
)

# Final adjudication
addStyle(
  wb,
  "Manual screening",
  style = final_header,
  rows = 1,
  cols = 15:18,
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
    cols = 1:18,
    gridExpand = TRUE
  )
  
  # Wrap long text columns
  addStyle(
    wb,
    "Manual screening",
    style = wrap_style,
    rows = 2:n_rows,
    cols = c(
      6,   # Title
      8,   # Abstract
      11,  # Auto_reason
      12,  # Gate4_evidence
      17,  # Final_reason
      18   # Notes
    ),
    gridExpand = TRUE,
    stack = TRUE
  )
}


# ------------------------------------------------------------
# 7. SET COLUMN WIDTHS
# ------------------------------------------------------------

# Study information
setColWidths(wb, "Manual screening", cols = 1, widths = 24)  # DOI
setColWidths(wb, "Manual screening", cols = 2, widths = 18)  # EID
setColWidths(wb, "Manual screening", cols = 3, widths = 12)  # PMID
setColWidths(wb, "Manual screening", cols = 4, widths = 22)  # Reference
setColWidths(wb, "Manual screening", cols = 5, widths = 9)   # Year
setColWidths(wb, "Manual screening", cols = 6, widths = 55)  # Title
setColWidths(wb, "Manual screening", cols = 7, widths = 28)  # Journal
setColWidths(wb, "Manual screening", cols = 8, widths = 80)  # Abstract

# Automated screening
setColWidths(wb, "Manual screening", cols = 9, widths = 16)   # Auto_verdict
setColWidths(wb, "Manual screening", cols = 10, widths = 14) # Auto_gate
setColWidths(wb, "Manual screening", cols = 11, widths = 45) # Auto_reason
setColWidths(wb, "Manual screening", cols = 12, widths = 45) # Gate4_evidence

# Independent reviewers
setColWidths(wb, "Manual screening", cols = 13:14, widths = 18)

# Final adjudication
setColWidths(wb, "Manual screening", cols = 15, widths = 16)  # Final_verdict
setColWidths(wb, "Manual screening", cols = 16, widths = 14)  # Final_gate
setColWidths(wb, "Manual screening", cols = 17, widths = 45)  # Final_reason
setColWidths(wb, "Manual screening", cols = 18, widths = 40)  # Notes


# ------------------------------------------------------------
# 8. FREEZE HEADER + IDENTIFYING COLUMNS
# ------------------------------------------------------------

freezePane(
  wb,
  sheet = "Manual screening",
  firstActiveRow = 2,
  firstActiveCol = 5
)


# ------------------------------------------------------------
# 9. ADD DATA VALIDATION DROPDOWNS
# ------------------------------------------------------------

if (n_rows > 1) {
  
  # Reviewer decisions
  dataValidation(
    wb,
    "Manual screening",
    cols = 13:14,
    rows = 2:n_rows,
    type = "list",
    value = '"Include,Exclude,Unsure"'
  )
  
  # Final verdict
  dataValidation(
    wb,
    "Manual screening",
    cols = 16,
    rows = 2:n_rows,
    type = "list",
    value = '"Include,Exclude"'
  )
  
  # Final gate
  dataValidation(
    wb,
    "Manual screening",
    cols = 16,
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
    "Auto_verdict",
    "Auto_gate",
    "Auto_reason",
    "Gate4_evidence",
    "Reviewer 1",
    "Reviewer 2",
    "Final_verdict",
    "Final_gate",
    "Final_reason",
    "Notes"
  ),
  
  Instructions = c(
    "Automated screening decision based on title, abstract and publication metadata.",
    "Automated gate assigned where there is positive evidence for exclusion.",
    "Reason assigned by the automated screening pipeline.",
    "Title/abstract evidence relevant to Gate 4 data availability. This field is for prioritisation only and does not constitute an automatic Gate 4 exclusion.",
    "Independent manual decision: Include, Exclude, or Unsure.",
    "Independent manual decision: Include, Exclude, or Unsure.",
    "Consensus decision after manual review.",
    "Eligibility gate responsible for final exclusion.",
    "Concise standardised reason for the final decision.",
    "Optional free-text comments, ambiguities, or adjudication notes."
  ),
  
  stringsAsFactors = FALSE
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
  widths = 90
)

addStyle(
  wb,
  "Instructions",
  style = study_header,
  rows = 1,
  cols = 1:2,
  gridExpand = TRUE
)

# Wrap instruction text
addStyle(
  wb,
  "Instructions",
  style = wrap_style,
  rows = 2:(nrow(instructions) + 1),
  cols = 2,
  gridExpand = TRUE
)


# ------------------------------------------------------------
# 11. SAVE WORKBOOK
# ------------------------------------------------------------

output_file <- "Scopus_screening_manual_review.xlsx"

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
# NOTE TO SELF: SCREENING HIERARCHY
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
# |                           |                | Scopus title/abstract alone.                                 |
# +---------------------------+----------------+---------------------------------------------------------------+
# | High-priority manual      | 58             | Subset of manual-review records with positive                |
# |                           |                | evidence for:                                                |
# |                           |                |   • brain tissue                                             |
# |                           |                |   • peripheral tissue                                        |
# |                           |                |   • matched/cross-tissue sampling                            |
# |                           |                |   • genome-wide DNA methylation                              |
# |                           |                | These are reviewed first.                                    |
# +---------------------------+----------------+---------------------------------------------------------------+
#
# Within the high-priority manual group, Gate 4 evidence is
# additionally used to prioritise records where the title or
# abstract suggests that downloadable CpG-level results,
# supplementary data, or relevant public data resources exist.
#
# Absence of Gate 4 evidence in the title/abstract is NOT
# sufficient for automatic exclusion.
###############################################################
#.   N Scopus records retrieved
#.   │
#.   ├── automatically excluded
#.   │
#.   └── not automatically excludable
#.   │
#.   ├── high-priority subset
#.   │    Candidate papers to screen
#.   │
#.   └── remaining records: insufficient information in title/abstract to automatically exclude


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
# CREATE INTERACTIVE SANKEY DIAGRAM OF SCREENING PIPELINE
#
# Input:
#   screening dataframe
#
# Output:
#   Scopus_screening_sankey.html
#
# Features:
#   - Wes Anderson Zissou1 colour palette
#   - exclusion counts displayed directly in node labels
#   - concise Sankey labels
#   - full screening-gate key underneath figure
#   - responsive layout designed to fit within HTML window
###############################################################

# Install once if required:
# install.packages(c(
#   "networkD3",
#   "htmlwidgets",
#   "htmltools",
#   "wesanderson"
# ))

library(networkD3)
library(htmlwidgets)
library(htmltools)
library(wesanderson)
library(dplyr)
library(stringr)


###############################################################
# 1. CHECK REQUIRED INPUTS
###############################################################

required_sankey_columns <- c(
  "EID",
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
# 2. COUNT AUTOMATED EXCLUSIONS
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
# 3. CALCULATE FLOW THROUGH EACH GATE
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

n_total_auto_excluded <- sum(
  screening$Auto_verdict == "Exclude",
  na.rm = TRUE
)


###############################################################
# 5. INTEGRITY CHECKS
###############################################################

if (n_other_manual < 0) {
  stop(
    "Sankey count error: high-priority manual records exceed ",
    "total manual-review records."
  )
}

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
# 6. CREATE SHORT EXCLUSION LABELS
#
# These are deliberately shorter than Auto_reason so that the
# diagram remains readable. The full eligibility criteria are
# given in the key below the Sankey.
###############################################################

sankey_exclusions <- sankey_exclusions %>%
  mutate(
    
    short_reason = case_when(
      
      Auto_reason == "Review article" ~
        "Review article",
      
      Auto_reason == "Meta-analysis" ~
        "Meta-analysis",
      
      Auto_reason == "Candidate-gene study" ~
        "Candidate-gene study",
      
      Auto_reason ==
        "Secondary analysis of previously published data" ~
        "Secondary analysis",
      
      Auto_reason == "Abstract only" ~
        "Abstract only",
      
      Auto_reason == "Non-human tissue/sample" ~
        "Non-human",
      
      Auto_reason ==
        "Neuroimaging phenotype; no brain-tissue DNAm" ~
        "No brain-tissue DNAm",
      
      Auto_reason ==
        "Brain and peripheral samples explicitly unmatched" ~
        "Unmatched tissues",
      
      Auto_reason ==
        "Differential methylation only / DMP-DMR analysis" ~
        "Differential methylation only",
      
      Auto_reason ==
        "Brain-tissue DNAm only; no peripheral DNAm" ~
        "Brain DNAm only",
      
      Auto_reason ==
        "Peripheral DNAm only; no brain-tissue DNAm" ~
        "Peripheral DNAm only",
      
      Auto_reason ==
        "Peripheral DNAm associated with brain-related phenotype; no brain-tissue DNAm" ~
        "Peripheral DNAm + brain phenotype",
      
      Auto_reason ==
        "Epigenetic-age/clock study without brain-tissue DNAm" ~
        "Epigenetic clock only",
      
      Auto_reason ==
        "Targeted/non-genome-wide methylation study" ~
        "Non-genome-wide DNAm",
      
      TRUE ~ Auto_reason
    ),
    
    # Label displayed in the Sankey itself
    node_label = paste0(
      Auto_gate,
      ": ",
      short_reason,
      " (n = ",
      format(n, big.mark = ","),
      ")"
    )
  )


###############################################################
# 7. CREATE CORE NODES WITH COUNTS
###############################################################

core_nodes <- tibble(
  
  internal_name = c(
    "Unique Scopus records",
    "Gate 1 passed",
    "Gate 2 passed",
    "Gate 3 passed",
    "Manual review",
    "High-priority manual",
    "Other manual review"
  ),
  
  name = c(
    paste0(
      "Unique Scopus records (n = ",
      format(n_unique, big.mark = ","),
      ")"
    ),
    
    paste0(
      "Gate 1 passed (n = ",
      format(n_after_gate1, big.mark = ","),
      ")"
    ),
    
    paste0(
      "Gate 2 passed (n = ",
      format(n_after_gate2, big.mark = ","),
      ")"
    ),
    
    paste0(
      "Gate 3 passed (n = ",
      format(n_after_gate3, big.mark = ","),
      ")"
    ),
    
    paste0(
      "Manual review (n = ",
      format(n_manual, big.mark = ","),
      ")"
    ),
    
    paste0(
      "High-priority manual (n = ",
      format(n_priority, big.mark = ","),
      ")"
    ),
    
    paste0(
      "Other manual review (n = ",
      format(n_other_manual, big.mark = ","),
      ")"
    )
  )
)


###############################################################
# 8. CREATE EXCLUSION NODES
###############################################################

reason_nodes <- sankey_exclusions %>%
  transmute(
    internal_name = paste0(
      Auto_gate,
      "::",
      Auto_reason
    ),
    
    name = node_label
  )


###############################################################
# 9. COMBINE NODES
###############################################################

sankey_nodes <- bind_rows(
  core_nodes,
  reason_nodes
)


###############################################################
# 10. CREATE LINKS
###############################################################

sankey_links <- tibble(
  source_name = character(),
  target_name = character(),
  value = numeric()
)


###############################################################
# Gate 1
###############################################################

gate1_reasons <- sankey_exclusions %>%
  filter(Auto_gate == "Gate 1")

if (nrow(gate1_reasons) > 0) {
  
  sankey_links <- bind_rows(
    
    sankey_links,
    
    gate1_reasons %>%
      transmute(
        source_name = "Unique Scopus records",
        target_name = paste0(
          Auto_gate,
          "::",
          Auto_reason
        ),
        value = n
      )
  )
}

sankey_links <- bind_rows(
  sankey_links,
  tibble(
    source_name = "Unique Scopus records",
    target_name = "Gate 1 passed",
    value = n_after_gate1
  )
)


###############################################################
# Gate 2
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
          Auto_gate,
          "::",
          Auto_reason
        ),
        value = n
      )
  )
}

sankey_links <- bind_rows(
  sankey_links,
  tibble(
    source_name = "Gate 1 passed",
    target_name = "Gate 2 passed",
    value = n_after_gate2
  )
)


###############################################################
# Gate 3
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
          Auto_gate,
          "::",
          Auto_reason
        ),
        value = n
      )
  )
}

sankey_links <- bind_rows(
  sankey_links,
  tibble(
    source_name = "Gate 2 passed",
    target_name = "Gate 3 passed",
    value = n_after_gate3
  )
)


###############################################################
# Gate 4
#
# Gate 4 is normally assessed manually in this pipeline.
# If automated Gate 4 exclusions are introduced later, they
# will automatically appear here.
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
          Auto_gate,
          "::",
          Auto_reason
        ),
        value = n
      )
  )
}


###############################################################
# Remaining records enter manual review
###############################################################

sankey_links <- bind_rows(
  sankey_links,
  tibble(
    source_name = "Gate 3 passed",
    target_name = "Manual review",
    value = n_manual
  )
)


###############################################################
# Split manual review into priority groups
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
# 11. REMOVE UNUSED NODES
###############################################################

used_node_names <- unique(
  c(
    sankey_links$source_name,
    sankey_links$target_name
  )
)

sankey_nodes <- sankey_nodes %>%
  filter(internal_name %in% used_node_names)


###############################################################
# 12. CONVERT NODE NAMES TO ZERO-BASED INDICES
###############################################################

sankey_links <- sankey_links %>%
  mutate(
    
    source = match(
      source_name,
      sankey_nodes$internal_name
    ) - 1,
    
    target = match(
      target_name,
      sankey_nodes$internal_name
    ) - 1
  )


###############################################################
# 13. ASSIGN NODE GROUPS
###############################################################

sankey_nodes <- sankey_nodes %>%
  mutate(
    
    group = case_when(
      
      internal_name == "Unique Scopus records" ~
        "Search",
      
      internal_name == "Gate 1 passed" ~
        "Passed",
      
      internal_name == "Gate 2 passed" ~
        "Passed",
      
      internal_name == "Gate 3 passed" ~
        "Passed",
      
      str_detect(internal_name, "^Gate 1::") ~
        "Gate1",
      
      str_detect(internal_name, "^Gate 2::") ~
        "Gate2",
      
      str_detect(internal_name, "^Gate 3::") ~
        "Gate3",
      
      str_detect(internal_name, "^Gate 4::") ~
        "Gate4",
      
      internal_name %in% c(
        "Manual review",
        "High-priority manual",
        "Other manual review"
      ) ~
        "Manual",
      
      TRUE ~
        "Other"
    )
  )


###############################################################
# 14. ASSIGN LINK GROUPS
###############################################################

sankey_links <- sankey_links %>%
  mutate(
    
    group = case_when(
      
      str_detect(target_name, "^Gate 1::") ~
        "Gate1",
      
      str_detect(target_name, "^Gate 2::") ~
        "Gate2",
      
      str_detect(target_name, "^Gate 3::") ~
        "Gate3",
      
      str_detect(target_name, "^Gate 4::") ~
        "Gate4",
      
      target_name %in% c(
        "Manual review",
        "High-priority manual",
        "Other manual review"
      ) ~
        "Manual",
      
      TRUE ~
        "Passed"
    )
  )


###############################################################
# 15. ZISSOU1 COLOUR PALETTE
#
# Zissou1:
#   #3B9AB2
#   #78B7C5
#   #EBCC2A
#   #E1AF00
#   #F21A00
#
# Interpolate to six colours so that each conceptual group
# can be distinguished while remaining within Zissou1.
###############################################################


colour_scale <- '
d3.scaleOrdinal()
  .domain([
    "Search",
    "Passed",
    "Gate1",
    "Gate2",
    "Gate3",
    "Gate4",
    "Manual"
  ])
  .range([
    "#3B9AB2",
    "#78B7C5",
    "#F21A00",
    "#E67E5F",
    "#F4A259",
    "#A5C2A3",
    "#BDC881"
  ])
'
    
    
  #  ),
  
#  "#3B9AB2",  # Search
#  "#78B7C5",  # Passed
#  "#91BAB6",  # Gate 1
#  "#A5C2A3",  # Gate 2
#  "#BDC881",   # Gate 3
#  "#F21A00",  # Gate 4
#  "#F4A259"	  # Manual review
#)


## other oranges to try
#"#E69F00"	
#"#E1AF00"
#"#D98C3F"	
#"#D97B29"
#"#F4A259"	
#"#E67E5F"

###############################################################
# 16. CREATE SANKEY
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
  
  colourScale = colour_scale,
  
  fontSize = 11,
  fontFamily = "Arial",
  
  nodeWidth = 18,
  nodePadding = 10,
  
  sinksRight = TRUE,
  
  width = NULL,
  height = 650
)

###############################################################
# 18. CREATE GATE KEY
###############################################################

gate_key <- tags$div(
  
  style = "
    max-width: 1350px;
    margin: 25px auto 10px auto;
    font-family: Arial, sans-serif;
  ",
  
  tags$h3(
    style = "margin-bottom: 10px;",
    "Screening gate definitions"
  ),
  
  tags$table(
    
    style = "
      width: 100%;
      border-collapse: collapse;
      font-size: 13px;
      line-height: 1.4;
    ",
    
    tags$tr(
      tags$th(
        style = "
          width: 27%;
          padding: 8px;
          text-align: left;
          border-bottom: 2px solid #444;
        ",
        "Gate"
      ),
      
      tags$th(
        style = "
          padding: 8px;
          text-align: left;
          border-bottom: 2px solid #444;
        ",
        "Eligibility question"
      )
    ),
    
    
    tags$tr(
      
      tags$td(
        style = "
          padding: 8px;
          vertical-align: top;
          border-bottom: 1px solid #ddd;
          font-weight: bold;
        ",
        "Gate 1. Appropriate study design?"
      ),
      
      tags$td(
        style = "
          padding: 8px;
          vertical-align: top;
          border-bottom: 1px solid #ddd;
        ",
        paste0(
          "Is this a primary research study on HUMAN participants ",
          "of genome-wide DNA methylation (rather than a ",
          "candidate-gene study, review, meta-analysis, conference ",
          "abstract or secondary analysis)?"
        )
      )
    ),
    
    
    tags$tr(
      
      tags$td(
        style = "
          padding: 8px;
          vertical-align: top;
          border-bottom: 1px solid #ddd;
          font-weight: bold;
        ",
        "Gate 2. Appropriate tissue sampling?"
      ),
      
      tags$td(
        style = "
          padding: 8px;
          vertical-align: top;
          border-bottom: 1px solid #ddd;
        ",
        paste0(
          "Does the study measure matched brain and peripheral ",
          "DNA methylation in the same individuals?"
        )
      )
    ),
    
    
    tags$tr(
      
      tags$td(
        style = "
          padding: 8px;
          vertical-align: top;
          border-bottom: 1px solid #ddd;
          font-weight: bold;
        ",
        "Gate 3. Appropriate analysis?"
      ),
      
      tags$td(
        style = "
          padding: 8px;
          vertical-align: top;
          border-bottom: 1px solid #ddd;
        ",
        paste0(
          "Does the study report CpG-level brain-peripheral ",
          "correlation or concordance statistics, rather than ",
          "only differential methylation (e.g. DMPs/DMRs) or ",
          "conventional EWAS analyses?"
        )
      )
    ),
    
    
    tags$tr(
      
      tags$td(
        style = "
          padding: 8px;
          vertical-align: top;
          font-weight: bold;
        ",
        "Gate 4. Appropriate data availability?"
      ),
      
      tags$td(
        style = "
          padding: 8px;
          vertical-align: top;
        ",
        paste0(
          "Are downloadable genome-wide CpG-level concordance ",
          "statistics available for extraction without ",
          "re-analysing raw data?"
        )
      )
    )
  )
)


###############################################################
# 19. CREATE TITLE AND SCREENING SUMMARY
###############################################################

screening_header <- tags$div(
  
  style = "
    max-width: 1350px;
    margin: 0 auto 5px auto;
    font-family: Arial, sans-serif;
  ",
  
  tags$h2(
    style = "margin-bottom: 5px;",
    "Scopus screening flow"
  ),
  
  tags$p(
    
    style = "
      font-size: 14px;
      margin-top: 0;
      margin-bottom: 8px;
    ",
    
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


###############################################################
# 20. ADD HEADER AND KEY TO SANKEY
###############################################################

screening_sankey <- htmlwidgets::prependContent(
  screening_sankey,
  screening_header
)

screening_sankey <- htmlwidgets::appendContent(
  screening_sankey,
  gate_key
)


###############################################################
# 21. SAVE INTERACTIVE HTML
###############################################################

sankey_output_file <- "Scopus_screening_sankey.html"

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
  "Scopus screening pipeline completed successfully\n",
  "=========================================================\n\n",

  "SEARCH RESULTS\n",
  "--------------\n",
  "Records retrieved:              ",
  format(n_scopus_retrieved, big.mark = ","), "\n",
  "Unique Scopus records:          ",
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
  "✓ Scopus_screening_records.csv\n",
  "✓ Scopus_screening_manual.csv\n",
  "✓ Scopus_screening_automatic_exclusions.csv\n",
  "✓ Scopus_screening_gate_summary.csv\n",
  "✓ Scopus_screening_exclusion_dictionary.csv\n",
  "✓ Scopus_screening_priority.csv\n\n",

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
