###############################################################
# Search_ewas_catalog.R
#
# Purpose:
# Identify EWAS Catalog publications containing analyses of
# both brain and peripheral DNA methylation.
#
# Input:
#   EWAS Catalog data downloaded 11/08/2026:
#   '"~/repos/CpG_Tissue_Corr/shiny_app/data/ewascatalog-studies.txt.gz"'
#   data/raw/ewascatalog/ewas_catalog_metadata_YYYY-MM-DD.tsv
#
# Output:
#   Candidate brain EWAS publications for manual eligibility screening
###############################################################

library(readr)
library(dplyr)
library(stringr)

###############################################################
# 1. Load EWAS Catalog metadata
###############################################################
setwd("/Users/eleanorconole/repos/CpG_Tissue_Corr/shiny_app/data")
print(getwd())

ewas <- read_tsv(
  "ewascatalog-studies.txt.gz",
  show_col_types = FALSE
)

###############################################################
# 2. Inspect tissue terminology used by the Catalog
###############################################################

tissue_counts <- ewas %>%
  count(tissue, sort = TRUE)

print(tissue_counts)

write_csv(
  tissue_counts,
  "EWAS_Catalog_tissue_labels.csv"
)

###############################################################
# 2a. DATA CLEANING OF EWAS CATALOG TISSUE TERMS
###############################################################
# n.b. separate tissue-label cleaning from study de-duplication! 

# Variants owing to use of capitals such as:
# 'Whole blood' | 'whole blood' | and 'Whole Blood'
# should be treated as the same tissue label, but they 
# are not automatically the same study.
# The same applies to Placenta/placenta, 
# Cord blood/cord blood/Cord Blood, Blood/blood, 
# Peripheral blood/Peripheral Blood, and 
# Superior temporal gyrus/superior temporal gyrus etc.

###############################################################
# 2b Preserve original tissue label and create cleaned version
###############################################################

ewas_clean <- ewas %>%
  mutate(
    Tissue_original = tissue,
    
    Tissue_clean = tissue %>%
      str_to_lower() %>%          # remove capitalisation differences
      str_squish() %>%            # remove duplicated whitespace
      str_trim()                  # remove leading/trailing whitespace
  )


#  Check whether capitalization variants represent
#  different studies / analyses

label_audit <- ewas_clean %>%
  group_by(Tissue_clean) %>%
  summarise(
    n_rows = n(),
    
    n_original_labels = n_distinct(Tissue_original),
    
    original_labels = paste(
      sort(unique(Tissue_original)),
      collapse = " | "
    ),
    
    n_studyIDs = n_distinct(study_id),
    
    n_PMIDs = n_distinct(pmid),
    
    .groups = "drop"
  ) %>%
  arrange(desc(n_rows))

# OK, this audit tells me that Capitalisation variants 
# are NOT duplicate rows of the same StudyID; they are
# (as anticiapted) inconsistent labels used across different analyses/studies.


###############################################################
# Tissue cleaning and harmonisation
###############################################################

ewas_clean <- ewas_clean %>%
  mutate(
    Has_whole_blood = str_detect(
      Tissue_clean,
      fixed("whole blood")
    )
  )




library(dplyr)
library(stringr)

ewas_clean <- ewas %>%
  mutate(
    # Preserve exactly as supplied by EWAS Catalog
    Tissue_original = tissue,
    
    # Lexical normalisation only
    Tissue_clean = Tissue_original %>%
      str_to_lower() %>%
      str_squish() %>%
      str_trim(),
    
    # Canonical display label
    Tissue_harmonised = case_when(

# PERIPHERAL-based terms      
      Tissue_original =="Whole blood" ~ "peripheral",                                                           
      Tissue_original =="whole blood"~ "peripheral",                                                           
      Tissue_original =="Cord blood"~ "peripheral",                                                            
      Tissue_original =="Placenta" ~ "peripheral",                                                             
      Tissue_original =="Whole Blood"~ "peripheral",                                                           
      Tissue_original =="Peripheral blood" ~ "peripheral",                                                     
      Tissue_original =="CD4+ T-cells"~ "peripheral",                                                          
      Tissue_original =="Leukocytes" ~ "peripheral",                                                           
                                                                     
      Tissue_original == "Saliva"~ "peripheral",                                                                
      Tissue_original == "cord blood" ~ "peripheral",                                                           
      Tissue_original == "Buccal cells"~ "peripheral",                                                          
      Tissue_original == "Skin"  ~ "peripheral",                                                                
      Tissue_original == "Airway epithelial cells"~ "peripheral",                                               
      Tissue_original == "placenta"~ "peripheral",                                                              
      Tissue_original == "Buccal"~ "peripheral",                                                                
      Tissue_original == "Monocytes"~ "peripheral",                                                             
      Tissue_original == "blood"~ "peripheral",                                                                 
      Tissue_original == "nasal epithelium"~ "peripheral",                                                      
      Tissue_original == "Adipose"~ "peripheral",                                                               
      Tissue_original == "Blood"~ "peripheral",                                                                 
      Tissue_original == "Cord Blood"~ "peripheral",                                                            
      Tissue_original == "Peripheral blood mononuclear cells"~ "peripheral",                                    
      Tissue_original == "Umbilical cord blood cells at birth and blood cells at 5y visit" ~ "peripheral",      
      Tissue_original == "Whole blood, CD4+ T-cells, CD14+ monocytes"  ~ "peripheral",                          
      Tissue_original == "Lung" ~ "peripheral",                                                                 
      Tissue_original == "Lung adenocarcinoma and lung squamous cell"~ "peripheral",                            
      Tissue_original == "Prostate tumor"  ~ "peripheral",                                                      
      Tissue_original == "Whole blood, heel prick blood spots"~ "peripheral",                                   
      Tissue_original == "B cells" ~ "peripheral",                                                              
      Tissue_original == "Breast"~ "peripheral",                                                                
      Tissue_original == "CD19+ B cells"~ "peripheral",                                                         
      Tissue_original == "Cord Tissue" ~ "peripheral",                                                          
      Tissue_original == "Liver"~ "peripheral",                                                                 
      Tissue_original == "Lymphoblasts"~ "peripheral",                                                          
      Tissue_original == "Maternal white blood cells"  ~ "peripheral",                                          
      Tissue_original == "Nasal epithelium" ~ "peripheral",                                                     
      Tissue_original == "T cells" ~ "peripheral",                                                              
      Tissue_original == "Unspecified"~ "peripheral",                                                           
      Tissue_original == "Whole blood, CD4+ T cells"   ~ "peripheral",                                          
      Tissue_original == "Whole blood, CD4+ T cells or monocytes" ~ "peripheral",                               
      Tissue_original == "nasal polyp" ~ "peripheral",                                                          
      Tissue_original == "Anal tumour" ~ "peripheral",                                                          
      Tissue_original == "Buccal cells and peripheral blood mononuclear cells" ~ "peripheral",                  
      Tissue_original == "CD14+ monocytes"  ~ "peripheral",                                                     
      Tissue_original == "CD4+ T-cells and leukocytes" ~ "peripheral",                                          
      Tissue_original == "CD4+ T-cells from cord blood"~ "peripheral",                                          
      Tissue_original == "CD4+ T-cells, leukocytes and whole blood"  ~ "peripheral",                            
      Tissue_original == "Cervical smear"    ~ "peripheral",                                                    
      Tissue_original == "Clear cell renal carcinoma tumour cells and adjacent healthy cells" ~ "peripheral",   
      Tissue_original == "Colon adenoma" ~ "peripheral",                                                        
      Tissue_original == "Colon adenomas" ~ "peripheral",                                                       
      Tissue_original == "Colon, rectum" ~ "peripheral",                                                        
      Tissue_original == "Cord blood, whole blood"  ~ "peripheral",                                             
      Tissue_original == "Fetal liver, adult liver" ~ "peripheral",                                             
      Tissue_original == "Leukocytes and CD4+ T cells" ~ "peripheral",                                          
      Tissue_original == "Leukocytes, whole blood and CD4+ T cells"  ~ "peripheral",                            
      Tissue_original == "Lung adenocarcinoma"   ~ "peripheral",                                                
      Tissue_original == "Lung cancer cells and normal cells"  ~ "peripheral",                                  
      Tissue_original == "Lung squamous cell"  ~ "peripheral",                                                  
      Tissue_original == "Lymphocytes" ~ "peripheral",                                                          
      Tissue_original == "Neutrophils"  ~ "peripheral",                                                         
      Tissue_original == "Pancreatic ductal adenocarcinoma and adjacent nontransformed pancreata"~ "peripheral",
      Tissue_original == "Peripheral Blood"     ~ "peripheral",                                                 
      Tissue_original == "Prostate cancer tissue, benign prostate tissue" ~ "peripheral",                       
      Tissue_original == "Skeletal muscle" ~ "peripheral",                                                      
      Tissue_original == "Umbilical artery" ~ "peripheral",                                                     
      Tissue_original == "Urine"   ~ "peripheral",                                                              
      Tissue_original == "White blood cells"   ~ "peripheral",                                                  
      Tissue_original == "Whole Blood, Monocytes" ~ "peripheral",                                               
      Tissue_original == "Whole blood and CD4+ T cells"  ~ "peripheral",                                        
      Tissue_original == "Whole blood and cord blood" ~ "peripheral",                                           
      Tissue_original == "Whole blood, CD4+ T cells, CD14+ monocytes"    ~ "peripheral",                        
      Tissue_original == "Whole blood, breast tissue"  ~ "peripheral",
      
      # BRAIN-based terms
      Tissue_original == "Fetal brain" ~ "brain",
      Tissue_original == "Brain" ~ "brain",
      Tissue_original == "Brain cortex" ~ "brain",
      Tissue_original == "Cerebellar tissue" ~ "brain",
      Tissue_original == "Cerebellum"~ "brain",
      Tissue_original == "Cortex"~ "brain",
      Tissue_original == "Cross-cortex"~ "brain",
      Tissue_original == "Dorsolateral prefrontal cortex"~ "brain",
      Tissue_original == "Entorhinal cortex"~ "brain",
      Tissue_original == "inferior frontal gyrus"~ "brain",
      Tissue_original == "Nucleus accumbens"~ "brain",
      
      Tissue_original == "Prefrontal cortex" ~ "brain",
      Tissue_original == "Cerebral cortex" ~ "brain",
      Tissue_original == "Superior temporal gyrus" ~ "brain",
      Tissue_original == "superior temporal gyrus" ~ "brain",
      
      
      
      TRUE ~ Tissue_original
    )
  )



tissue_counts <- ewas_clean %>%
  count(Tissue_harmonised, sort = TRUE)

tissue_counts

brain_EWAS <- ewas_clean %>%
  filter(Tissue_harmonised == "brain")

#   Candidate brain EWAS publications for manual eligibility screening

unique_brain_EWAS <- brain_EWAS %>% 
  distinct(pmid, .keep_all = TRUE)

#   We have 13 Candidate brain EWAS publications for manual screening


###############################################################
# Running a check through the peripheral EWAS to see if
# any look at brain-related phenotypes in 'trait' in case
# these EWASs contain matched brain-tissue data
###############################################################

peripheral_EWAS <- ewas_clean %>%
  filter(Tissue_harmonised == "peripheral") %>% 
  distinct(pmid, .keep_all = TRUE)


unique_peripheral_EWAS_traits <- peripheral_EWAS %>% 
  distinct(trait, .keep_all = TRUE)


unique_peripheral_EWAS_traits_brain_phenotype <- unique_peripheral_EWAS_traits$trait
###############################################################
# End
###############################################################
