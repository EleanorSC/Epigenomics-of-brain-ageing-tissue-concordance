###############################################################
# Extract_GEO_screening.R
#
# Purpose:
# Search NCBI Gene Expression Omnibus (GEO) for human DNA
# methylation Series containing both brain and peripheral
# samples, and identify potentially relevant Series that are
# not linked to an indexed PubMed publication.
#
# This is a SUPPLEMENTARY study-identification workflow.
#
# The script:
#   1. searches Entrez GEO DataSets for human methylation GSEs;
#   2. requires evidence of both brain and peripheral sample
#      sources at the GEO Series level;
#   3. retrieves GEO Series metadata programmatically;
#   4. cross-references each GEO record against PubMed via ELink;
#   5. optionally inspects GEO sample-level metadata using
#      GEOquery to strengthen tissue classification and identify
#      donor/subject/sample-matching terminology;
#   6. flags GEO Series with no linked PubMed publication;
#   7. exports candidate Series for manual eligibility review.
#
# IMPORTANT:
# Evidence that a GEO Series contains brain and peripheral
# samples is NOT sufficient to establish eligibility.
# The final eligibility requirements still include:
#   - human samples;
#   - matched brain and peripheral tissue from the same people;
#   - genome-wide DNA methylation;
#   - CpG-level cross-tissue correlation/concordance statistics;
#   - downloadable results permitting extraction without
#     requiring a new analysis of raw data.
#
# A GEO Series without a linked PMID is therefore a candidate
# for manual investigation, not an automatically eligible study.
###############################################################

# Install CRAN packages once if required:
# install.packages(c(
#   "rentrez", "dplyr", "stringr", "purrr", "readr", "tibble"
# ))
#
# GEOquery is a Bioconductor package:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("GEOquery")

###############################################################
# Start timer
###############################################################

start_time <- Sys.time()

cat(
  "\n=========================================================\n",
  "XTRACE-CpG / Brain–Peripheral DNAm Concordance Project\n",
  "Supplementary GEO Search\n",
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
    sprintf("[%d/7] %s\n", stage_number, stage_name),
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
library(dplyr)
library(stringr)
library(purrr)
library(readr)
library(tibble)

has_geoquery <- requireNamespace(
  "GEOquery",
  quietly = TRUE
)

if (!has_geoquery) {
  warning(
    "GEOquery is not installed. The Entrez GEO search and ",
    "PubMed cross-reference will still run, but sample-level ",
    "metadata inspection will be skipped."
  )
}

###############################################################
# 1. Define the GEO search strategy
#
# NCBI GEO DataSets supports:
#   GSE[ETYP]       = GEO Series only
#   [ORGN]          = organism
#   [DESC]          = descriptive metadata
#   [SRC]           = Sample Source
#
# We deliberately search broadly for methylation terminology,
# but require BOTH brain and peripheral sample-source concepts.
###############################################################

methylation_terms <- c(
  "\"DNA methylation\"[DESC]",
  "methylation[DESC]",
  "methylome[DESC]",
  "methylation[All Fields]",
  "\"bisulfite sequencing\"[DESC]",
  "\"methylation profiling\"[DESC]"
)

brain_source_terms <- c(
  "brain[SRC]",
  "cortex[SRC]",
  "cortical[SRC]",
  "prefrontal[SRC]",
  "hippocampus[SRC]",
  "hippocampal[SRC]",
  "cerebellum[SRC]",
  "cerebellar[SRC]",
  "amygdala[SRC]",
  "striatum[SRC]",
  "putamen[SRC]",
  "caudate[SRC]",
  "\"nucleus accumbens\"[SRC]",
  "\"substantia nigra\"[SRC]",
  "\"superior temporal gyrus\"[SRC]",
  "\"entorhinal cortex\"[SRC]"
)

peripheral_source_terms <- c(
  "blood[SRC]",
  "\"whole blood\"[SRC]",
  "\"peripheral blood\"[SRC]",
  "leukocyte[SRC]",
  "leukocytes[SRC]",
  "PBMC[SRC]",
  "saliva[SRC]",
  "buccal[SRC]",
  "cheek[SRC]",
  "serum[SRC]",
  "plasma[SRC]",
  "placenta[SRC]",
  "\"cord blood\"[SRC]"
)

geo_query <- paste0(
  "Homo sapiens[ORGN] AND GSE[ETYP] AND ",
  "(",
  paste(methylation_terms, collapse = " OR "),
  ") AND ",
  "(",
  paste(brain_source_terms, collapse = " OR "),
  ") AND ",
  "(",
  paste(peripheral_source_terms, collapse = " OR "),
  ")"
)

cat(
  "\nGEO QUERY\n",
  "=========\n",
  geo_query,
  "\n\n",
  sep = ""
)

###############################################################
# 2. Search GEO DataSets
###############################################################

t <- stage_start(1, "Searching GEO DataSets")

geo_search <- entrez_search(
  db = "gds",
  term = geo_query,
  retmax = 10000,
  use_history = TRUE
)

cat(
  "  Candidate GEO Series returned: ",
  format(geo_search$count, big.mark = ","),
  "\n",
  sep = ""
)

if (geo_search$count == 0) {
  stop(
    "The GEO search returned zero candidate Series. ",
    "Review the GEO query before continuing."
  )
}

stage_end(
  t,
  paste0(
    format(geo_search$count, big.mark = ","),
    " candidate GEO Series identified"
  )
)

###############################################################
# 3. Retrieve GEO Series summaries
###############################################################

t <- stage_start(2, "Retrieving GEO Series metadata")

geo_summaries <- entrez_summary(
  db = "gds",
  id = geo_search$ids,
  always_return_list = TRUE
)

###############################################################
# Helper for safely extracting Entrez summary fields
###############################################################

summary_field <- function(x, candidates) {

  nms <- names(x)

  if (is.null(nms)) {
    return(NA_character_)
  }

  hit <- candidates[candidates %in% nms]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  value <- x[[hit[1]]]

  if (is.null(value) || length(value) == 0) {
    return(NA_character_)
  }

  value <- paste(
    as.character(unlist(value, use.names = FALSE)),
    collapse = "; "
  )

  if (identical(value, "")) {
    return(NA_character_)
  }

  value
}

geo_records <- map2_dfr(
  geo_summaries,
  geo_search$ids,
  function(x, uid) {

    tibble(
      GEO_UID = as.character(uid),

      GSE = summary_field(
        x,
        c(
          "accession",
          "Accession",
          "gse",
          "GSE"
        )
      ),

      Title = summary_field(
        x,
        c(
          "title",
          "Title"
        )
      ),

      Summary = summary_field(
        x,
        c(
          "summary",
          "Summary",
          "description",
          "Description"
        )
      ),

      Organism = summary_field(
        x,
        c(
          "taxon",
          "Taxon",
          "organism",
          "Organism"
        )
      ),

      GDS_type = summary_field(
        x,
        c(
          "gdstype",
          "GDS_Type",
          "type",
          "Type"
        )
      ),

      Sample_count = summary_field(
        x,
        c(
          "n_samples",
          "samples",
          "SampleCount"
        )
      ),

      GPL = summary_field(
        x,
        c(
          "gpl",
          "GPL"
        )
      ),

      GEO_release_date = summary_field(
        x,
        c(
          "PDAT",
          "pdat",
          "PublicationDate"
        )
      )
    )
  }
)

###############################################################
# Some Entrez summaries expose the Series accession in a
# slightly different structure. Recover GSE accessions from
# the complete summary text if required.
###############################################################

for (i in seq_len(nrow(geo_records))) {

  if (
    is.na(geo_records$GSE[i]) ||
    !str_detect(geo_records$GSE[i], "^GSE[0-9]+$")
  ) {

    summary_text <- paste(
      unlist(
        geo_summaries[[i]],
        use.names = FALSE
      ),
      collapse = " "
    )

    recovered_gse <- str_extract(
      summary_text,
      "\\bGSE[0-9]+\\b"
    )

    if (!is.na(recovered_gse)) {
      geo_records$GSE[i] <- recovered_gse
    }
  }
}

###############################################################
# Retrieval integrity checks
###############################################################

if (any(is.na(geo_records$GSE))) {
  warning(
    sum(is.na(geo_records$GSE)),
    " GEO record(s) are missing a recoverable GSE accession."
  )
}

stage_end(
  t,
  paste0(
    format(nrow(geo_records), big.mark = ","),
    " GEO Series summaries retrieved"
  )
)

###############################################################
# 4. Cross-reference GEO records against PubMed
#
# NCBI ELink can identify PubMed records linked to a GEO
# DataSets UID. This is the key step for identifying GEO
# deposits that do not yet have a linked indexed publication.
###############################################################

t <- stage_start(3, "Cross-referencing GEO Series against PubMed")

get_linked_pmids <- function(gds_uid) {

  result <- tryCatch(
    {
      entrez_link(
        dbfrom = "gds",
        db = "pubmed",
        id = gds_uid
      )
    },
    error = function(e) {
      warning(
        "PubMed link lookup failed for GEO UID ",
        gds_uid,
        ": ",
        conditionMessage(e)
      )
      return(NULL)
    }
  )

  if (is.null(result)) {
    return(NA_character_)
  }

  links <- result$links

  if (is.null(links) || length(links) == 0) {
    return("")
  }

  pmid_vectors <- links[
    str_detect(
      names(links),
      regex("pubmed", ignore_case = TRUE)
    )
  ]

  if (length(pmid_vectors) == 0) {
    return("")
  }

  pmids <- unique(
    as.character(
      unlist(
        pmid_vectors,
        use.names = FALSE
      )
    )
  )

  pmids <- pmids[
    str_detect(pmids, "^[0-9]+$")
  ]

  if (length(pmids) == 0) {
    return("")
  }

  paste(pmids, collapse = "; ")
}

geo_records$Linked_PMIDs <- map_chr(
  geo_records$GEO_UID,
  get_linked_pmids
)

geo_records <- geo_records %>%
  mutate(
    Has_linked_PubMed = case_when(
      is.na(Linked_PMIDs) ~ NA,
      Linked_PMIDs == "" ~ FALSE,
      TRUE ~ TRUE
    )
  )

stage_end(
  t,
  paste0(
    sum(geo_records$Has_linked_PubMed %in% TRUE, na.rm = TRUE),
    " Series linked to PubMed; ",
    sum(geo_records$Has_linked_PubMed %in% FALSE, na.rm = TRUE),
    " Series with no linked PubMed record"
  )
)

###############################################################
# 5. Inspect sample-level GEO metadata
#
# This stage is designed to distinguish:
#   - Series that genuinely contain brain samples;
#   - Series that genuinely contain peripheral samples;
#   - Series whose metadata contains donor/subject/pairing
#     terminology suggestive of matched sampling.
#
# It does NOT automatically infer that tissues came from the
# same individuals. That requires manual inspection.
###############################################################

t <- stage_start(4, "Inspecting GEO sample-level metadata")

brain_pattern <- regex(
  paste(
    c(
      "\\bbrain\\b",
      "\\bcortex\\b",
      "\\bcortical\\b",
      "prefrontal",
      "dorsolateral prefrontal",
      "entorhinal",
      "superior temporal gyrus",
      "inferior frontal gyrus",
      "hippocamp",
      "cerebell",
      "amygdala",
      "nucleus accumbens",
      "striat",
      "putamen",
      "caudate",
      "substantia nigra"
    ),
    collapse = "|"
  ),
  ignore_case = TRUE
)

peripheral_pattern <- regex(
  paste(
    c(
      "\\bwhole blood\\b",
      "\\bperipheral blood\\b",
      "\\bblood\\b",
      "\\bleukocyte",
      "\\bPBMC",
      "\\bsaliva\\b",
      "\\bbuccal\\b",
      "\\bcheek\\b",
      "\\bserum\\b",
      "\\bplasma\\b",
      "\\bplacenta\\b",
      "\\bcord blood\\b"
    ),
    collapse = "|"
  ),
  ignore_case = TRUE
)

matching_pattern <- regex(
  paste(
    c(
      "\\bmatched\\b",
      "\\bpaired\\b",
      "same individual",
      "same participant",
      "same subject",
      "same donor",
      "donor id",
      "donor_id",
      "subject id",
      "subject_id",
      "participant id",
      "participant_id",
      "individual id",
      "individual_id",
      "patient id",
      "patient_id"
    ),
    collapse = "|"
  ),
  ignore_case = TRUE
)

sample_level_results <- vector(
  "list",
  nrow(geo_records)
)

if (has_geoquery) {

  for (i in seq_len(nrow(geo_records))) {

    gse_id <- geo_records$GSE[i]

    cat(
      "  Inspecting ",
      gse_id,
      " (",
      i,
      "/",
      nrow(geo_records),
      ")...\n",
      sep = ""
    )

    if (
      is.na(gse_id) ||
      gse_id == ""
    ) {
      next
    }

    gse_object <- tryCatch(
      {
        GEOquery::getGEO(
          gse_id,
          GSEMatrix = FALSE,
          AnnotGPL = FALSE,
          getGPL = FALSE
        )
      },
      error = function(e) {
        warning(
          "Could not retrieve sample-level metadata for ",
          gse_id,
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (is.null(gse_object)) {
      next
    }

    gsm_list <- tryCatch(
      GEOquery::GSMList(gse_object),
      error = function(e) list()
    )

    if (length(gsm_list) == 0) {
      next
    }

    gsm_df <- imap_dfr(
      gsm_list,
      function(gsm, gsm_id) {

        meta <- tryCatch(
          GEOquery::Meta(gsm),
          error = function(e) list()
        )

        meta_text <- paste(
          unlist(
            meta,
            use.names = FALSE
          ),
          collapse = " | "
        )

        tibble(
          GSE = gse_id,
          GSM = gsm_id,
          Sample_metadata = meta_text,

          Is_brain = str_detect(
            meta_text,
            brain_pattern
          ),

          Is_peripheral = str_detect(
            meta_text,
            peripheral_pattern
          ),

          Has_matching_language = str_detect(
            meta_text,
            matching_pattern
          )
        )
      }
    )

    sample_level_results[[i]] <- gsm_df

    Sys.sleep(0.2)
  }

  geo_samples <- bind_rows(
    sample_level_results
  )

} else {

  geo_samples <- tibble(
    GSE = character(),
    GSM = character(),
    Sample_metadata = character(),
    Is_brain = logical(),
    Is_peripheral = logical(),
    Has_matching_language = logical()
  )
}

###############################################################
# Summarise sample-level evidence by GEO Series
###############################################################

if (nrow(geo_samples) > 0) {

  sample_summary <- geo_samples %>%
    group_by(GSE) %>%
    summarise(
      n_GSM = n(),
      n_brain_samples = sum(Is_brain, na.rm = TRUE),
      n_peripheral_samples = sum(Is_peripheral, na.rm = TRUE),
      n_samples_with_matching_language =
        sum(Has_matching_language, na.rm = TRUE),

      Sample_level_brain_and_peripheral =
        any(Is_brain, na.rm = TRUE) &
        any(Is_peripheral, na.rm = TRUE),

      Matching_language_present =
        any(Has_matching_language, na.rm = TRUE),

      .groups = "drop"
    )

  geo_records <- geo_records %>%
    left_join(
      sample_summary,
      by = "GSE"
    )

} else {

  geo_records <- geo_records %>%
    mutate(
      n_GSM = NA_integer_,
      n_brain_samples = NA_integer_,
      n_peripheral_samples = NA_integer_,
      n_samples_with_matching_language = NA_integer_,
      Sample_level_brain_and_peripheral = NA,
      Matching_language_present = NA
    )
}

stage_end(
  t,
  if (has_geoquery) {
    paste0(
      format(nrow(geo_samples), big.mark = ","),
      " GEO samples inspected"
    )
  } else {
    "Sample-level inspection skipped because GEOquery is unavailable"
  }
)

###############################################################
# 6. Candidate classification
###############################################################

t <- stage_start(5, "Classifying GEO candidates")

geo_records <- geo_records %>%
  mutate(

    Series_text = str_to_lower(
      paste(
        coalesce(Title, ""),
        coalesce(Summary, ""),
        coalesce(GDS_type, "")
      )
    ),

    Series_mentions_methylation = str_detect(
      Series_text,
      regex(
        "DNA methylation|methylation|methylome|bisulfite|epigenome",
        ignore_case = TRUE
      )
    ),

    Series_mentions_brain = str_detect(
      Series_text,
      brain_pattern
    ),

    Series_mentions_peripheral = str_detect(
      Series_text,
      peripheral_pattern
    ),

    ###########################################################
    # Main triage status
    ###########################################################

    GEO_triage = case_when(

      Has_linked_PubMed == FALSE &
        Sample_level_brain_and_peripheral %in% TRUE &
        Matching_language_present %in% TRUE ~
        "HIGH PRIORITY: brain + peripheral samples; matching terminology; no linked PMID",

      Has_linked_PubMed == FALSE &
        Sample_level_brain_and_peripheral %in% TRUE ~
        "PRIORITY: brain + peripheral samples; no linked PMID",

      Has_linked_PubMed == FALSE ~
        "No linked PMID; inspect GEO Series manually",

      Has_linked_PubMed == TRUE &
        Sample_level_brain_and_peripheral %in% TRUE ~
        "Brain + peripheral samples; linked publication exists",

      TRUE ~
        "Inspect manually"
    ),

    ###########################################################
    # Specific flag requested for the preregistered workflow:
    # candidate matched brain-peripheral methylation Series not
    # yet linked to an indexed PubMed publication.
    ###########################################################

    Candidate_unpublished_or_unlinked =
      Has_linked_PubMed == FALSE &
      (
        Sample_level_brain_and_peripheral %in% TRUE |
        (
          Series_mentions_brain &
          Series_mentions_peripheral
        )
      )
  )

geo_unlinked_candidates <- geo_records %>%
  filter(
    Candidate_unpublished_or_unlinked %in% TRUE
  ) %>%
  arrange(
    desc(Matching_language_present),
    GSE
  )

geo_high_priority <- geo_records %>%
  filter(
    Candidate_unpublished_or_unlinked %in% TRUE,
    Matching_language_present %in% TRUE
  ) %>%
  arrange(GSE)

geo_linked_candidates <- geo_records %>%
  filter(
    Has_linked_PubMed %in% TRUE
  ) %>%
  arrange(GSE)

stage_end(
  t,
  paste0(
    nrow(geo_unlinked_candidates),
    " potentially relevant GEO Series have no linked PubMed publication"
  )
)

###############################################################
# 7. Optional cross-reference against existing project outputs
#
# If your PubMed / Scopus screening files are present in the
# current working directory, linked PMIDs are compared against
# those outputs. This is supplementary bookkeeping only.
###############################################################

t <- stage_start(6, "Cross-referencing existing screening outputs")

read_existing_ids <- function(file) {

  if (!file.exists(file)) {
    return(
      tibble(
        PMID = character(),
        DOI = character()
      )
    )
  }

  x <- suppressMessages(
    readr::read_csv(
      file,
      show_col_types = FALSE
    )
  )

  tibble(
    PMID = if ("PMID" %in% names(x)) {
      as.character(x$PMID)
    } else {
      NA_character_
    },

    DOI = if ("DOI" %in% names(x)) {
      str_to_lower(
        as.character(x$DOI)
      )
    } else {
      NA_character_
    }
  )
}

existing_pubmed <- read_existing_ids(
  "PubMed_screening_records.csv"
)

existing_scopus <- read_existing_ids(
  "Scopus_screening_records.csv"
)

existing_pmids <- unique(
  na.omit(
    c(
      existing_pubmed$PMID,
      existing_scopus$PMID
    )
  )
)

geo_records <- geo_records %>%
  mutate(
    Linked_PMID_already_in_project = map_lgl(
      Linked_PMIDs,
      function(x) {

        if (
          is.na(x) ||
          x == "" ||
          length(existing_pmids) == 0
        ) {
          return(FALSE)
        }

        pmids <- str_split(
          x,
          ";\\s*"
        )[[1]]

        any(
          pmids %in% existing_pmids
        )
      }
    )
  )

geo_unlinked_candidates <- geo_records %>%
  filter(
    Candidate_unpublished_or_unlinked %in% TRUE
  ) %>%
  arrange(
    desc(Matching_language_present),
    GSE
  )

stage_end(
  t,
  paste0(
    length(existing_pmids),
    " existing PubMed/Scopus PMIDs available for cross-reference"
  )
)

###############################################################
# 8. Write outputs
###############################################################

t <- stage_start(7, "Writing GEO search outputs")

write_csv(
  geo_records,
  "GEO_brain_peripheral_methylation_all_candidates.csv"
)

write_csv(
  geo_unlinked_candidates,
  "GEO_unlinked_brain_peripheral_candidates.csv"
)

write_csv(
  geo_high_priority,
  "GEO_unlinked_high_priority_candidates.csv"
)

write_csv(
  geo_linked_candidates,
  "GEO_candidates_with_linked_publications.csv"
)

if (nrow(geo_samples) > 0) {
  write_csv(
    geo_samples,
    "GEO_candidate_sample_metadata.csv"
  )
}

###############################################################
# Create a simple manual-review sheet
###############################################################

geo_manual_review <- geo_records %>%
  transmute(
    GSE,
    Title,
    Linked_PMIDs,
    Has_linked_PubMed,
    n_GSM,
    n_brain_samples,
    n_peripheral_samples,
    Matching_language_present,
    GEO_triage,

    Gate1_human_genomewide_DNAm = "",
    Gate2_matched_brain_peripheral = "",
    Gate3_CpG_concordance_analysis = "",
    Gate4_downloadable_CpG_results = "",

    Final_verdict = "",
    Final_reason = "",
    Notes = ""
  ) %>%
  mutate(
    Candidate_order = str_detect(
      GEO_triage,
      "^HIGH PRIORITY|^PRIORITY|^No linked"
    )
  ) %>%
  arrange(
    desc(Candidate_order),
    GSE
  ) %>%
  select(-Candidate_order)

write_csv(
  geo_manual_review,
  "GEO_manual_eligibility_review.csv"
)

stage_end(
  t,
  "GEO CSV outputs written"
)

###############################################################
# Final summary
###############################################################

end_time <- Sys.time()

runtime <- as.numeric(
  difftime(
    end_time,
    start_time,
    units = "secs"
  )
)

hours <- runtime %/% 3600
minutes <- (runtime %% 3600) %/% 60
seconds <- round(runtime %% 60)

cat(
  "\n=========================================================\n",
  "GEO supplementary search completed\n",
  "=========================================================\n\n",

  "SEARCH RESULTS\n",
  "--------------\n",
  "GEO Series returned:                    ",
  format(nrow(geo_records), big.mark = ","), "\n",

  "Series with linked PubMed publication: ",
  format(
    sum(
      geo_records$Has_linked_PubMed %in% TRUE,
      na.rm = TRUE
    ),
    big.mark = ","
  ),
  "\n",

  "Series without linked PubMed record:    ",
  format(
    sum(
      geo_records$Has_linked_PubMed %in% FALSE,
      na.rm = TRUE
    ),
    big.mark = ","
  ),
  "\n",

  "Unlinked brain-peripheral candidates:   ",
  format(
    nrow(geo_unlinked_candidates),
    big.mark = ","
  ),
  "\n",

  "High-priority unlinked candidates:      ",
  format(
    nrow(geo_high_priority),
    big.mark = ","
  ),
  "\n\n",

  "OUTPUT FILES\n",
  "------------\n",
  "✓ GEO_brain_peripheral_methylation_all_candidates.csv\n",
  "✓ GEO_unlinked_brain_peripheral_candidates.csv\n",
  "✓ GEO_unlinked_high_priority_candidates.csv\n",
  "✓ GEO_candidates_with_linked_publications.csv\n",
  if (nrow(geo_samples) > 0) {
    "✓ GEO_candidate_sample_metadata.csv\n"
  } else {
    ""
  },
  "✓ GEO_manual_eligibility_review.csv\n\n",

  "IMPORTANT\n",
  "---------\n",
  "A GEO Series with no linked PMID is not assumed to be ",
  "unpublished or eligible. It is flagged for manual checking ",
  "because the corresponding publication may be unindexed, ",
  "not linked by GEO, in press, or absent.\n\n",

  "Total runtime: ",
  sprintf(
    "%02d:%02d:%02d",
    hours,
    minutes,
    seconds
  ),
  "\n",

  "Finished: ",
  format(
    end_time,
    "%Y-%m-%d %H:%M:%S"
  ),
  "\n",
  "=========================================================\n",
  sep = ""
)
