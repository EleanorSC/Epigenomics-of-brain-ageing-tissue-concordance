###############################################################
# Search_EWAS_catalog.R
#
# Purpose:
# Filter the EWAS Catalog for all eligible studies involving brain |
# peripheral tissues for screening.
###############################################################

library(readr)
library(dplyr)
library(stringr)

###############################################################
# Load EWAS Catalog metadata
###############################################################

ewas <- read_tsv(
  "ewas_catalog.tsv",
  show_col_types = FALSE
)

###############################################################
# Inspect available tissue labels
###############################################################

sort(unique(ewas$Tissue))

###############################################################
# Define tissues of interest
###############################################################

tissue_pattern <- paste(
  c(
    "brain",
    "blood",
    "whole blood",
    "peripheral blood",
    "pbmc",
    "saliva",
    "buccal"
  ),
  collapse = "|"
)

###############################################################
# Filter records
###############################################################

ewas_tissue <- ewas %>%
  filter(
    str_detect(
      string = Tissue,
      pattern = regex(tissue_pattern, ignore_case = TRUE)
    )
  )

###############################################################
# View tissues identified
###############################################################

ewas_tissue %>%
  count(Tissue, sort = TRUE)

###############################################################
# Create unique publication list
###############################################################

publications <- ewas_tissue %>%
  distinct(
    PMID,
    Author,
    Date,
    Trait,
    Tissue,
    Methylation_Array,
    N,
    .keep_all = TRUE
  ) %>%
  arrange(desc(Date))

###############################################################
# Summary
###############################################################

cat("Total EWAS records:", nrow(ewas), "\n")
cat("Filtered records:", nrow(ewas_tissue), "\n")
cat("Unique publications:", nrow(publications), "\n")

###############################################################
# Save outputs
###############################################################

write_csv(
  ewas_tissue,
  "EWAS_Catalog_Tissue_Filtered.csv"
)

write_csv(
  publications,
  "EWAS_Catalog_Publications_For_Screening.csv"
)

###############################################################
# Frequency table for reporting
###############################################################

tissue_summary <- ewas_tissue %>%
  count(Tissue, sort = TRUE)

write_csv(
  tissue_summary,
  "EWAS_Tissue_Summary.csv"
)

###############################################################
# End
###############################################################
