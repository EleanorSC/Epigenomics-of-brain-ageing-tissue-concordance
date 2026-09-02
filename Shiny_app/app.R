
# TRACE-CpG Shiny App prototype 31/08/26 E.L.S Conole
# Prototype interactive web application for the TRACE-CpG project.
#
# This script:
#   - loads the harmonised CpG-level brain–peripheral DNAm correlation database;
#   - creates standardised brain-region and peripheral-tissue labels;
#   - creates searchable CpG- and gene-level database views;
#   - currently generates simple interactive plots of concordance distributions and
#     tissue/genomic characteristics;
#   - makes a first attempt at an interactive brain-region plot (glass-brain); and
#   - allows users to export filtered subsets of the CpG-level database.
#
# Each database row represents one CpG × brain–peripheral tissue comparison.
#
# Primary input files:
#   - data/cpg_brain_peripheral_correlation_database.csv
#   - data/gene_level_summary.csv
#
# ! NOTE:
# This is a preregistration prototype to accompany the OSF. Most elements will be
# updated post systematic screening and as the systematic TRACE-CpG project evolves.

# -------------------------------------------------------------------------
# Load Packages
# -------------------------------------------------------------------------
#install.packages("sf", type = "binary")
#install.packages("ggseg3d")

# -------------------------------------------------------------------------

library(sf)
library(ggseg3d)
library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(plotly)

# -------------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------------

setwd("/Users/eleanorconole/repos/CpG_Tissue_Corr")
print(getwd())

options(
  shiny.sanitize.errors = FALSE,
  shiny.maxRequestSize = 30 * 1024^2
)

find_database <- function() {
  candidates <- c(
    file.path("data", "cpg_brain_peripheral_correlation_database.csv"),
    "cpg_brain_peripheral_correlation_database.csv",
    file.path("..", "cpg_brain_peripheral_correlation_database.csv")
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop(
      paste0(
        "Database file not found. Checked:\n- ",
        paste(normalizePath(candidates, mustWork = FALSE), collapse = "\n- ")
      )
    )
  }
  existing[[1]]
}

database_path <- find_database()

cpg_data <- read.csv(
  database_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
) |>
  mutate(
    CpG_ID = cpg,
    Gene = gene,
    Correlation = suppressWarnings(as.numeric(corr_value)),
    Source = source,
    Brain_Region_Code = case_when(
      grepl("prefrontal|PFC", tissue_pair, ignore.case = TRUE) ~ "PFC",
      grepl("superior temporal|STG", tissue_pair, ignore.case = TRUE) ~ "STG",
      grepl("entorhinal|EC", tissue_pair, ignore.case = TRUE) ~ "EC",
      grepl("cerebell|CER", tissue_pair, ignore.case = TRUE) ~ "CER",
      TRUE ~ NA_character_
    ),
    Brain_Region = case_when(
      Brain_Region_Code == "PFC" ~ "Prefrontal cortex",
      Brain_Region_Code == "STG" ~ "Superior temporal gyrus",
      Brain_Region_Code == "EC" ~ "Entorhinal cortex",
      Brain_Region_Code == "CER" ~ "Cerebellum",
      TRUE ~ "Unclassified"
    ),
    Peripheral_Tissue = case_when(
      grepl("blood", tissue_pair, ignore.case = TRUE) ~ "Blood",
      grepl("buccal", tissue_pair, ignore.case = TRUE) ~ "Buccal",
      grepl("saliva", tissue_pair, ignore.case = TRUE) ~ "Saliva",
      TRUE ~ "Other"
    )
  )


# Database-tab data: gene-annotated rows only.
clean_gene <- function(x) {
  vapply(
    x,
    function(value) {
      if (is.na(value) || !nzchar(trimws(value))) return("")
      parts <- trimws(strsplit(value, ";")[[1]])
      parts <- parts[nzchar(parts)]
      paste(unique(parts), collapse = "; ")
    },
    character(1),
    USE.NAMES = FALSE
  )
}

cpg_db <- cpg_data |>
  mutate(
    Gene = clean_gene(Gene),
    Chromosome = ifelse(
      grepl("^chr", chr, ignore.case = TRUE),
      chr,
      paste0("chr", chr)
    )
  ) |>
  filter(!is.na(Gene), trimws(Gene) != "")

find_gene_summary <- function() {
  candidates <- c(
    file.path("data", "gene_level_summary.csv"),
    "gene_level_summary.csv",
    file.path("..", "gene_level_summary.csv")
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop(
      paste0(
        "Gene summary file not found. Checked:\n- ",
        paste(normalizePath(candidates, mustWork = FALSE), collapse = "\n- ")
      )
    )
  }
  existing[[1]]
}

gene_summary <- read.csv(
  find_gene_summary(),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gene_choices <- sort(unique(unlist(strsplit(cpg_db$Gene, "; "))))
gene_choices <- gene_choices[nzchar(gene_choices)]

chromosome_choices <- sort(
  unique(cpg_db$Chromosome),
  na.last = TRUE
)

tissue_pair_choices <- sort(unique(cpg_db$tissue_pair))
context_choices <- sort(unique(cpg_db$context))
island_choices <- sort(unique(cpg_db$island_relation))

brain_region_lookup <- c(
  PFC = "Prefrontal cortex",
  STG = "Superior temporal gyrus",
  EC = "Entorhinal cortex",
  CER = "Cerebellum"
)


# -------------------------------------------------------------------------
# Glass-brain locations
# -------------------------------------------------------------------------
#
# I need to teplace these coordinates with verified MNI coordinates before publication.
# NOTE These rows are placeholders for testing the visualisation only.

glass_brain_locations <- tibble::tribble(
  ~region_code, ~region_name,                  ~x,  ~y,  ~z,
  "PFC",        "Prefrontal cortex",           -35,  42,  28,
  "STG",        "Superior temporal gyrus",     -54, -20,   4,
  "EC",         "Entorhinal cortex",           -24, -18, -28,
  "CER",        "Cerebellum",                  -24, -58, -32
)


# -------------------------------------------------------------------------
# Interactive Brain Cortices function
# -------------------------------------------------------------------------

interactive_brain <- function() {
  div(
    class = "brain-panel",
    div(
      class = "brain-heading",
      h2("Select a brain region"),
      p("Click a region to filter the CpG concordance database.")
    ),
    div(
      class = "brain-layout",
      tags$svg(
        id = "brain-map",
        class = "brain-map",
        viewBox = "0 0 760 440",
        role = "img",
        `aria-labelledby` = "brain-title brain-description",
        tags$title(id = "brain-title", "Interactive diagram of brain regions"),
        tags$desc(
          id = "brain-description",
          paste(
            "Clickable regions include prefrontal cortex,",
            "superior temporal gyrus, entorhinal cortex and cerebellum."
          )
        ),
        tags$path(
          class = "brain-outline",
          d = paste(
            "M 113 234",
            "C 79 183, 87 115, 137 77",
            "C 181 43, 247 39, 301 62",
            "C 354 31, 430 37, 474 72",
            "C 533 73, 589 111, 606 166",
            "C 624 226, 603 282, 560 313",
            "C 525 339, 477 346, 430 334",
            "C 384 369, 316 374, 267 340",
            "C 218 355, 162 333, 135 294",
            "C 122 276, 116 255, 113 234 Z"
          )
        ),
        tags$g(
          class = "brain-region",
          `data-region` = "PFC",
          tabindex = "0",
          role = "button",
          `aria-label` = "Prefrontal cortex",
          tags$path(
            d = paste(
              "M 113 234",
              "C 79 183, 87 115, 137 77",
              "C 169 52, 208 43, 246 48",
              "C 229 92, 226 139, 238 181",
              "C 220 220, 211 260, 218 304",
              "C 183 309, 151 294, 135 270",
              "C 124 254, 118 242, 113 234 Z"
            )
          ),
          tags$text(x = "164", y = "172", class = "brain-label", "PFC")
        ),
        tags$g(
          class = "brain-region",
          `data-region` = "STG",
          tabindex = "0",
          role = "button",
          `aria-label` = "Superior temporal gyrus",
          tags$path(
            d = paste(
              "M 219 235",
              "C 275 215, 350 216, 417 236",
              "C 455 248, 493 270, 515 299",
              "C 487 327, 446 342, 405 333",
              "C 357 357, 299 355, 263 330",
              "C 239 312, 225 277, 219 235 Z"
            )
          ),
          tags$text(x = "340", y = "286", class = "brain-label", "STG")
        ),
        tags$g(
          class = "brain-region",
          `data-region` = "EC",
          tabindex = "0",
          role = "button",
          `aria-label` = "Entorhinal cortex",
          tags$path(
            d = paste(
              "M 302 282",
              "C 329 263, 370 262, 398 279",
              "C 421 293, 429 318, 414 337",
              "C 385 355, 344 356, 318 339",
              "C 299 326, 291 303, 302 282 Z"
            )
          ),
          tags$text(x = "356", y = "314", class = "brain-label brain-label-small", "EC")
        ),
        tags$path(
          class = "brain-context",
          d = paste(
            "M 241 50",
            "C 282 44, 311 51, 333 68",
            "C 372 42, 428 44, 468 73",
            "C 516 73, 564 104, 590 150",
            "C 606 180, 607 216, 596 245",
            "C 570 235, 540 238, 516 258",
            "C 483 236, 448 221, 412 211",
            "C 355 195, 294 197, 240 214",
            "C 229 167, 231 100, 241 50 Z"
          )
        ),
        tags$g(
          class = "brain-region",
          `data-region` = "CER",
          tabindex = "0",
          role = "button",
          `aria-label` = "Cerebellum",
          tags$path(
            d = paste(
              "M 502 282",
              "C 538 258, 585 263, 617 289",
              "C 650 316, 650 358, 618 383",
              "C 590 405, 544 402, 513 381",
              "C 486 361, 478 315, 502 282 Z"
            )
          ),
          tags$path(class = "cerebellar-line", d = "M 512 306 C 544 288, 588 291, 619 312"),
          tags$path(class = "cerebellar-line", d = "M 503 331 C 540 313, 593 317, 631 339"),
          tags$path(class = "cerebellar-line", d = "M 512 357 C 547 340, 590 344, 619 366"),
          tags$text(x = "565", y = "344", class = "brain-label brain-label-small", "CER")
        ),
        tags$path(
          class = "brain-context",
          d = paste(
            "M 463 318",
            "C 477 329, 491 347, 493 371",
            "L 489 417",
            "L 455 417",
            "L 451 372",
            "C 448 348, 451 330, 463 318 Z"
          )
        )
      ),
      div(
        class = "brain-information",
        div(
          class = "selected-region-card",
          span(class = "selected-region-eyebrow", "Selected region"),
          div(id = "selected-region-name", class = "selected-region-name", "All brain regions"),
          p(
            id = "selected-region-description",
            class = "selected-region-description",
            paste(
              "No regional filter is currently applied.",
              "Select a region on the brain diagram."
            )
          )
        ),
        actionButton(
          "clear_brain_region",
          "Show all regions",
          class = "btn-outline-primary"
        )
      )
    )
  )
}

ui <- page_navbar(
  id = "main_nav",
  title = "CpG Brain–Peripheral Concordance",
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#078894",
    base_font = font_google("Source Sans 3")
  ),
  header = tags$head(
    tags$link(rel = "stylesheet", href = "brain.css"),
    tags$script(src = "brain.js")
  ),
  # nav_panel(
   # "Search",
   # layout_sidebar(
    nav_panel(
      title = tagList(
        icon("magnifying-glass"),
        " Search"
      ),
      
      layout_sidebar(
      sidebar = sidebar(
        width = 310,
        open = "desktop",
        textInput(
          "query",
          "Gene symbol or CpG identifier",
          placeholder = "For example: AHRR or cg05575921"
        ),
        selectInput(
          "peripheral_tissue",
          "Peripheral tissue",
          choices = c("All tissues", sort(unique(cpg_data$Peripheral_Tissue)))
        ),
        sliderInput(
          "minimum_correlation",
          "Minimum absolute correlation",
          min = 0.60,
          max = 1.00,
          value = 0.60,
          step = 0.01
        ),
        checkboxInput(
          "gene_annotated_only",
          "Gene-annotated rows only",
          value = FALSE
        )
      ),
      div(
        class = "main-content",
        interactive_brain(),
        uiOutput("selection_summary"),
        card(
          class = "results-card plot-card",
          full_screen = TRUE,
          card_header(
            class = "results-card-header",
            div(
              class = "card-title-row",
              span("CpG concordance by tissue pair"),
              span(class = "card-subtitle", "Hover over points for CpG-level details")
            )
          ),
          card_body(
            class = "plot-card-body",
            plotlyOutput("correlation_plot", height = "340px")
          )
        ),
        card(
          class = "results-card table-card",
          full_screen = TRUE,
          card_header(
            class = "results-card-header",
            div(
              class = "card-title-row",
              span("CpG-level results"),
              span(
                class = "card-subtitle",
                "Each row is one CpG × tissue-pair comparison"
              )
            )
          ),
          card_body(
            class = "table-card-body",
            DTOutput("results_table")
          )
        )
      )
    )
  ),

  nav_panel(
    title = tagList(icon("table"), " CpG Database"),

    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        open = "desktop",

        selectizeInput(
          "db_gene",
          "Gene",
          choices = NULL,
          multiple = TRUE,
          options = list(
            placeholder = "Type to search a gene..."
          )
        ),

        selectInput(
          "db_chr",
          "Chromosome",
          choices = c("All", chromosome_choices),
          selected = "All",
          multiple = TRUE
        ),

        checkboxGroupInput(
          "db_tissue",
          "Tissue pair",
          choices = tissue_pair_choices,
          selected = tissue_pair_choices
        ),

        radioButtons(
          "db_corr_type",
          "Correlation type",
          choices = c(
            "Both",
            "Pearson r",
            "Spearman rho"
          ),
          selected = "Both"
        ),

        sliderInput(
          "db_corr_range",
          "|Correlation value|",
          min = 0.60,
          max = 1.00,
          value = c(0.60, 1.00),
          step = 0.01
        ),

        checkboxGroupInput(
          "db_context",
          "Genomic context",
          choices = context_choices,
          selected = context_choices
        ),

        checkboxGroupInput(
          "db_island",
          "CpG island relation",
          choices = island_choices,
          selected = island_choices
        ),

        actionButton(
          "db_reset",
          "Reset filters",
          icon = icon("rotate-left"),
          width = "100%"
        )
      ),

      div(
        class = "main-content",

        div(
          class = "summary-grid database-summary-grid",

          div(
            class = "metric-card",
            span(class = "metric-label", "Matching rows"),
            strong(class = "metric-value", textOutput("db_kpi_rows", inline = TRUE))
          ),

          div(
            class = "metric-card",
            span(class = "metric-label", "Unique CpGs"),
            strong(class = "metric-value", textOutput("db_kpi_cpgs", inline = TRUE))
          ),

          div(
            class = "metric-card",
            span(class = "metric-label", "Unique genes"),
            strong(class = "metric-value", textOutput("db_kpi_genes", inline = TRUE))
          )
        ),

        downloadButton(
          "download_db",
          "Download filtered rows (CSV)",
          icon = icon("download")
        ),

        card(
          class = "results-card table-card",
          full_screen = TRUE,

          card_header(
            class = "results-card-header",
            div(
              class = "card-title-row",
              span("CpG database"),
              span(
                class = "card-subtitle",
                "Gene-annotated CpG × tissue-pair observations"
              )
            )
          ),

          card_body(
            class = "table-card-body",
            DTOutput("db_table")
          )
        )
      )
    )
  ),

  nav_panel(
    title = tagList(icon("dna"), " Gene Summary"),

    div(
      class = "main-content",

      card(
        class = "results-card table-card",
        full_screen = TRUE,

        card_header(
          class = "results-card-header",
          div(
            class = "card-title-row",
            span("Gene-level summary"),
            span(
              class = "card-subtitle",
              "One row per gene; use the column filters to search"
            )
          )
        ),

        card_body(
          class = "table-card-body",
          DTOutput("gene_summary_table")
        )
      )
    )
  ),

  nav_panel(
    title = tagList(icon("chart-column"), " Visualise"),

    div(
      class = "main-content",

      div(
        class = "visualise-intro",
        h2("Visualise the filtered CpG database"),
        p(
          paste(
            "These plots update using the filters applied in the",
            "CpG Database tab."
          )
        )
      ),

      div(
        class = "visualise-grid",

        card(
          class = "results-card visual-card",
          full_screen = TRUE,
          card_header("Correlation distribution"),
          card_body(
            plotlyOutput("viz_corr_dist", height = "360px")
          )
        ),

        card(
          class = "results-card visual-card",
          full_screen = TRUE,
          card_header("Rows by tissue pair"),
          card_body(
            plotlyOutput("viz_tissue", height = "360px")
          )
        ),

        card(
          class = "results-card visual-card",
          full_screen = TRUE,
          card_header("Rows by genomic context"),
          card_body(
            plotlyOutput("viz_context", height = "360px")
          )
        ),

        card(
          class = "results-card visual-card",
          full_screen = TRUE,
          card_header("Rows by CpG-island relation"),
          card_body(
            plotlyOutput("viz_island", height = "360px")
          )
        )
      )
    )
  ),
  
  nav_panel(
    title = tagList(
      icon("brain"),
      " Glass Brain"
    ),
    
    div(
      class = "main-content",
      
      div(
        class = "visualise-intro",
        
        h2("Glass-brain view"),
        
        p(
          paste(
            "Rotate and zoom the brain to inspect the spatial locations",
            "represented in the cross-tissue resource."
          )
        )
      ),
      
      layout_columns(
        col_widths = c(8, 4),
        
        card(
          class = "results-card glass-brain-card",
          full_screen = TRUE,
          
          card_header(
            class = "results-card-header",
            
            div(
              class = "card-title-row",
              
              span("Brain sampling locations"),
              
              span(
                class = "card-subtitle",
                "Drag to rotate; scroll to zoom"
              )
            )
          ),
          
          card_body(
            class = "glass-brain-card-body",
            
          #  plotlyOutput(
          #    "glass_brain_plot",
          #    height = "620px"
          #  )
            
            ggseg3dOutput(
              "glass_brain_plot",
              height = "620px"
            )
          )
        ),
        
        card(
          class = "results-card glass-controls-card",
          
          card_header(
            class = "results-card-header",
            "Display options"
          ),
          
          card_body(
            
            selectInput(
              "glass_brain_region",
              "Brain region",
              choices = c(
                "All regions",
                setNames(
                  glass_brain_locations$region_code,
                  glass_brain_locations$region_name
                )
              ),
              selected = "All regions"
            ),
            
            sliderInput(
              "glass_brain_opacity",
              "Brain opacity",
              min = 0.05,
              max = 0.60,
              value = 0.20,
              step = 0.05
            ),
            
            sliderInput(
              "glass_point_size",
              "Marker size",
              min = 4,
              max = 18,
              value = 10,
              step = 1
            ),
            
            checkboxInput(
              "glass_show_labels",
              "Show region labels",
              value = TRUE
            ),
            
            hr(),
            
            div(
              class = "glass-brain-note",
              
              strong("Interpretation"),
              
              p(
                paste(
                  "Markers should denote verified sampling coordinates.",
                  "They do not represent genomic coordinates or the physical",
                  "location of CpG sites."
                )
              )
            )
          )
        )
      )
    )
  ),

  #nav_panel(
  #  "Methods",
  #  div(
  nav_panel(
    title = tagList(
      icon("book-open"),
      "Access Guide"
    ),
    
   # div(
   #   class = "methods-page",
   #   h2("About the database"),
   #   p(
   #     paste(
   #       "Each row represents one CpG × peripheral-tissue ×",
   #       "brain-region comparison."
   #     )
   #   ),
   #   h3("Interpretation"),
   #   p(
   #     paste(
   #       "Cross-tissue correlation does not establish that a peripheral",
   #       "signal measures the same biological process in brain.",
   #       "Concordance may reflect shared genetic regulation, cell composition,",
   #       "technical effects or stable between-person differences."
   #       
            
            div(
              class = "methods-page",
              
              h2("About the database"),
              
              div(
                class = "methods-summary-grid",
                
                div(
                  class = "methods-stat-card",
                  span(class = "methods-stat-value", "9,127"),
                  span(class = "methods-stat-label", "CpG-level observations")
                ),
                
                div(
                  class = "methods-stat-card",
                  span(class = "methods-stat-value", "815"),
                  span(class = "methods-stat-label", "Genes represented")
                ),
                
                div(
                  class = "methods-stat-card",
                  span(class = "methods-stat-value", "5"),
                  span(class = "methods-stat-label", "Primary resources")
                ),
                
                div(
                  class = "methods-stat-card",
                  span(class = "methods-stat-value", "4"),
                  span(class = "methods-stat-label", "Brain regions")
                ),
                
                div(
                  class = "methods-stat-card",
                  span(class = "methods-stat-value", "2"),
                  span(class = "methods-stat-label", "Peripheral tissues")
                ),
                
                div(
                  class = "methods-stat-card",
                  span(class = "methods-stat-value", "|r| > 0.60"),
                  span(class = "methods-stat-label", "Inclusion threshold")
                )
              ),
              
              h3("Why this resource was developed"),
              
              p(
                paste(
                  "This database brings together published evidence on cross-tissue",
                  "DNA methylation concordance between human brain tissue and",
                  "accessible peripheral tissues. Its purpose is to allow researchers",
                  "to search for a gene or CpG site and determine whether methylation",
                  "measured in blood or buccal tissue has previously shown strong",
                  "correlation with methylation in a brain region."
                )
              ),
              
              p(
                paste(
                  "Brain tissue is rarely accessible in living cohorts. Studies of",
                  "cognitive ageing, dementia and neuropsychiatric disease therefore",
                  "often rely on blood, saliva or buccal samples as surrogate tissues.",
                  "However, DNA methylation is strongly tissue-specific, and a",
                  "peripheral methylation signal should not be assumed to reflect brain",
                  "methylation without CpG-level evidence of concordance."
                )
              ),
              
              p(
                paste(
                  "The relevant evidence was previously distributed across several",
                  "non-interoperable journal supplements and interactive web tools.",
                  "The aim was to convert that scattered evidence into a single,",
                  "structured and gene-searchable resource while retaining clear",
                  "provenance for every observation."
                )
              ),
              
              hr(),
              
              h3("Primary resources"),
              
              p(
                paste(
                  "Five published resources were identified as the principal human",
                  "brain–peripheral DNA methylation concordance resources."
                )
              ),
              
              tags$ul(
                tags$li(
                  strong("Exeter Blood–Brain DNA Methylation Comparison Tool: "),
                  "blood compared with four postmortem brain regions, based on ",
                  tags$em("Hannon et al. (2015)")
                ),
                
                tags$li(
                  strong("BECon: "),
                  "blood compared with three Brodmann areas, based on ",
                  tags$em("Edgar et al. (2017)")
                ),
                
                tags$li(
                  strong("IMAGE-CpG: "),
                  "brain compared with blood, saliva and buccal tissue in a living ",
                  "epilepsy-surgery cohort, based on ",
                  tags$em("Braun et al. (2019)")
                ),
                
                tags$li(
                  strong("Buccal–Brain Correlation Map: "),
                  "buccal tissue compared with prefrontal cortex, based on ",
                  tags$em("Sommerer et al. (2022)")
                ),
                
                tags$li(
                  strong("AMAZE-CpG: "),
                  "brain compared with blood, saliva and buccal tissue in an ",
                  "independent Japanese cohort, based on ",
                  tags$em("Nishitani et al. (2023)")
                )
              ),
              
              hr(),
              
              h3("Data availability assessment"),
              
              p(
                paste(
                  "Each resource was checked to determine whether the underlying",
                  "CpG-level correlation results were available as a bulk-downloadable",
                  "supplementary dataset rather than only through an interactive",
                  "single-gene lookup interface."
                )
              ),
              
              p(
                paste(
                  "Bulk supplementary data were available for the Exeter",
                  "blood–brain resource and the Sommerer buccal–brain resource.",
                  "Equivalent downloadable correlation files were not identified for",
                  "IMAGE-CpG, BECon or AMAZE-CpG. These three resources are therefore",
                  "represented in the database access guide and in individually",
                  "verified literature examples, but not as bulk-extracted rows in the",
                  "main CpG-level table."
                )
              ),
              
              hr(),
              
              h3("Bulk data extraction"),
              
              tags$ol(
                tags$li(
                  "The shared CpG annotation manifest was loaded from Hannon et al. ",
                  "(2015), including chromosome, genomic position, gene annotation, ",
                  "UCSC RefGene group and CpG-island relation."
                ),
                
                tags$li(
                  "Blood–brain correlations were extracted from the Hannon ",
                  "supplementary tables covering four brain regions."
                ),
                
                tags$li(
                  "Buccal–brain correlations were extracted from the Sommerer ",
                  "supplementary table covering buccal tissue and prefrontal cortex."
                ),
                
                tags$li(
                  "Every CpG × tissue-pair observation with an absolute Pearson ",
                  "correlation or Spearman correlation greater than 0.60 was retained."
                ),
                
                tags$li(
                  "Genomic annotations were harmonised across the two source datasets."
                ),
                
                tags$li(
                  "The extracted rows were combined into one master table and assigned ",
                  "consistent tissue-pair labels, correlation-type labels, source ",
                  "citations and database-tool metadata."
                )
              ),
              
              hr(),
              
              h3("Inclusion threshold"),
              
              p(
                paste(
                  "The database retains CpG-level observations satisfying",
                  "|r| > 0.60 or |ρ| > 0.60. The threshold is applied to the",
                  "absolute value of the reported coefficient, so both strong positive",
                  "and strong negative cross-tissue correlations are eligible."
                )
              ),
              
              p(
                paste(
                  "The threshold identifies comparatively strong concordance but should",
                  "not be treated as a universal biological cut-off. The app allows",
                  "users to raise the minimum absolute-correlation threshold when a",
                  "more stringent search is required."
                )
              ),
              
              hr(),
              
              h3("Data harmonisation and annotation"),
              
              p(
                paste(
                  "Each retained record represents one CpG × peripheral-tissue ×",
                  "brain-region comparison. The same CpG can therefore appear in",
                  "multiple rows when it was evaluated across several brain regions or",
                  "tissue pairs."
                )
              ),
              
              tags$ul(
                tags$li("Illumina CpG probe identifier"),
                tags$li("Chromosome and genomic coordinate"),
                tags$li("Gene annotation"),
                tags$li("Promoter, gene-body or intergenic context"),
                tags$li("CpG-island relation"),
                tags$li("Peripheral tissue"),
                tags$li("Brain region"),
                tags$li("Correlation metric and coefficient"),
                tags$li("Source publication"),
                tags$li("Source URL and associated database tool")
              ),
              
              hr(),
              
              h3("Database outputs"),
              
              p(
                paste(
                  "The two bulk source datasets yielded 7,252 blood–brain observations",
                  "and 1,875 buccal–brain observations. Together they form a master",
                  "database of 9,127 CpG × tissue-pair rows."
                )
              ),
              
              p(
                paste(
                  "The CpG-level table was also collapsed to a one-row-per-gene summary",
                  "covering 815 genes. The gene-level summary reports whether a gene",
                  "contains one or more qualifying CpGs and enables rapid screening",
                  "before examining individual CpG records."
                )
              ),
              
              tags$ul(
                tags$li(
                  strong("CpG-level database: "),
                  "9,127 cross-tissue observations."
                ),
                
                tags$li(
                  strong("Gene-level summary: "),
                  "815 genes represented by at least one qualifying CpG."
                ),
                
                tags$li(
                  strong("Curated literature examples: "),
                  "individually verified CpGs reported in peer-reviewed studies that ",
                  "queried interactive concordance tools."
                ),
                
                tags$li(
                  strong("Database access guide: "),
                  "links and query instructions for all five primary resources."
                )
              ),
              
              hr(),
              
              h3("Literature-curated examples"),
              
              p(
                paste(
                  "A separate literature pass identified published studies that had",
                  "queried interactive resources such as BECon or the Exeter tool for",
                  "specific CpGs. These examples were verified individually and are",
                  "kept separate from the bulk-extracted master database."
                )
              ),
              
              p(
                paste(
                  "The curated examples include CpGs associated with AHRR, MAD1L1,",
                  "HEXDC, CDC42BPB, LDHC and the MIR886 cluster. They are presented as",
                  "illustrative literature evidence and are not merged into the",
                  "automatically extracted dataset."
                )
              ),
              
              hr(),
              
              h3("How to interpret a result"),
              
              p(
                paste(
                  "A high cross-tissue correlation indicates that between-person",
                  "variation in methylation at a particular CpG was similar across the",
                  "specified peripheral and brain tissues in the source study."
                )
              ),
              
              p(
                paste(
                  "It does not establish that the peripheral and brain measurements",
                  "reflect the same cellular mechanism, that changes in one tissue",
                  "cause changes in the other, or that the CpG is functionally relevant",
                  "to a neurological phenotype."
                )
              ),
              
              p(
                paste(
                  "Peripheral methylation should therefore be interpreted as a",
                  "potential surrogate marker rather than as a direct measurement of",
                  "brain methylation."
                )
              ),
              
              hr(),
              
              h3("Important sources of concordance"),
              
              p(
                paste(
                  "Strong cross-tissue correlations can arise for several reasons,",
                  "including shared developmental regulation, stable between-person",
                  "differences, genetic control of methylation, cell-composition",
                  "patterns and technical characteristics of the methylation assay."
                )
              ),
              
              p(
                paste(
                  "In particular, some very high correlations may reflect methylation",
                  "quantitative trait loci (mQTLs), where the same inherited genetic",
                  "variant influences methylation across multiple tissues. Such a CpG",
                  "may be a reliable cross-tissue marker without necessarily indexing",
                  "a shared environmentally responsive biological process."
                )
              ),
              
              hr(),
              
              h3("Coverage limitations"),
              
              tags$ul(
                tags$li(
                  "No bulk-downloadable saliva–brain correlation dataset was ",
                  "identified. Saliva–brain queries currently require gene-by-gene ",
                  "lookup through IMAGE-CpG or AMAZE-CpG."
                ),
                
                tags$li(
                  "Bulk CpG-level rows in this app derive from Hannon et al. and ",
                  "Sommerer et al.; the other resources are not represented as complete ",
                  "downloadable datasets."
                ),
                
                tags$li(
                  "The included studies differ in cohort composition, sample size, ",
                  "brain region, tissue source, methylation platform and correlation ",
                  "metric."
                ),
                
                tags$li(
                  "Most available evidence is based on bulk tissue, so observed ",
                  "correlations can be influenced by differences in cellular ",
                  "composition."
                ),
                
                tags$li(
                  "The absence of a CpG from the database does not establish absence of ",
                  "cross-tissue concordance. It may instead indicate that the CpG was ",
                  "not measured, not reported, did not exceed the chosen threshold or ",
                  "is only queryable through an external interactive tool."
                ),
                
                tags$li(
                  "The database reports published associations and cannot be used to ",
                  "infer causal relationships between peripheral and brain methylation."
                )
              ),
              
              hr(),
              
              h3("Reproducibility"),
              
              p(
                paste(
                  "The database was generated through a scripted workflow. The source",
                  "supplementary files are read and processed by",
                  "`database_extract_source_data.py`, which produces the standardised",
                  "CpG-level and gene-level CSV files."
                )
              ),
              
              p(
                paste(
                  "The final formatted workbook is generated separately by",
                  "`build_cpg_database.py`. This separation preserves a clear distinction",
                  "between primary-data extraction and workbook presentation."
                )
              ),
              
              p(
                paste(
                  "No manually curated CpG values are inserted into the bulk-extracted",
                  "master table. Literature-curated examples remain explicitly labelled",
                  "and separate."
                )
              ),
              
              hr(),
              
              h3("Recommended use"),
              
              p(
                paste(
                  "The resource is intended for study design, candidate-CpG appraisal,",
                  "interpretation of peripheral methylation findings and identification",
                  "of loci requiring direct cross-tissue validation."
                )
              ),
              
              p(
                paste(
                  "Researchers should cite both this database and the original source",
                  "publication from which each correlation was obtained."
                )
              )
            )
          )
          ###
          
#          
#        )
#      )
#    )
#  )
)

server <- function(input, output, session) {

  updateSelectizeInput(
    session,
    "db_gene",
    choices = gene_choices,
    server = TRUE
  )

  observeEvent(input$db_reset, {
    updateSelectizeInput(
      session,
      "db_gene",
      selected = character(0),
      choices = gene_choices,
      server = TRUE
    )
    updateSelectInput(session, "db_chr", selected = "All")
    updateCheckboxGroupInput(
      session,
      "db_tissue",
      selected = tissue_pair_choices
    )
    updateRadioButtons(session, "db_corr_type", selected = "Both")
    updateSliderInput(
      session,
      "db_corr_range",
      value = c(0.60, 1.00)
    )
    updateCheckboxGroupInput(
      session,
      "db_context",
      selected = context_choices
    )
    updateCheckboxGroupInput(
      session,
      "db_island",
      selected = island_choices
    )
  })

  db_filtered <- reactive({
    dat <- cpg_db

    if (length(input$db_gene) > 0) {
      selected_genes <- toupper(input$db_gene)
      dat <- dat |>
        filter(
          vapply(
            strsplit(toupper(Gene), "; "),
            function(x) any(x %in% selected_genes),
            logical(1)
          )
        )
    }

    if (
      !is.null(input$db_chr) &&
      length(input$db_chr) > 0 &&
      !("All" %in% input$db_chr)
    ) {
      dat <- dat |> filter(Chromosome %in% input$db_chr)
    }

    if (is.null(input$db_tissue) || length(input$db_tissue) == 0) {
      dat <- dat[0, ]
    } else {
      dat <- dat |> filter(tissue_pair %in% input$db_tissue)
    }

    if (
      !is.null(input$db_corr_type) &&
      input$db_corr_type != "Both"
    ) {
      dat <- dat |> filter(corr_type == input$db_corr_type)
    }

    dat <- dat |>
      filter(
        abs(Correlation) >= input$db_corr_range[1],
        abs(Correlation) <= input$db_corr_range[2]
      )

    if (is.null(input$db_context) || length(input$db_context) == 0) {
      dat <- dat[0, ]
    } else {
      dat <- dat |> filter(context %in% input$db_context)
    }

    if (is.null(input$db_island) || length(input$db_island) == 0) {
      dat <- dat[0, ]
    } else {
      dat <- dat |> filter(island_relation %in% input$db_island)
    }

    dat
  })

  output$db_kpi_rows <- renderText(
    format(nrow(db_filtered()), big.mark = ",")
  )

  output$db_kpi_cpgs <- renderText(
    format(n_distinct(db_filtered()$CpG_ID), big.mark = ",")
  )

  output$db_kpi_genes <- renderText({
    genes <- unlist(strsplit(db_filtered()$Gene, "; "))
    genes <- genes[nzchar(genes)]
    format(length(unique(genes)), big.mark = ",")
  })

  output$db_table <- renderDT({
    dat <- db_filtered()

    validate(
      need(nrow(dat) > 0, "No rows match the current database filters.")
    )

    table_data <- dat |>
      transmute(
        `CpG ID` = CpG_ID,
        Chromosome = Chromosome,
        `Position (hg19)` = pos,
        Gene = Gene,
        `Genomic context` = context,
        `Island relation` = island_relation,
        `Tissue pair` = tissue_pair,
        `Correlation type` = corr_type,
        `Correlation value` = Correlation,
        `Database/tool` = database_tool,
        Source = Source,
        `Source URL` = source_url
      )

    datatable(
      table_data,
      rownames = FALSE,
      filter = "top",
      extensions = c("Buttons", "FixedHeader"),
      class = "stripe hover compact",
      options = list(
        pageLength = 20,
        lengthMenu = c(10, 20, 50, 100),
        deferRender = TRUE,
        processing = TRUE,
        scrollX = TRUE,
        fixedHeader = TRUE,
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel")
      )
    ) |>
      formatRound("Correlation value", digits = 3)
  }, server = TRUE)

  output$download_db <- downloadHandler(
    filename = function() {
      paste0("cpg_database_filtered_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(db_filtered(), file, row.names = FALSE)
    }
  )

  output$gene_summary_table <- renderDT({
    display_data <- gene_summary

    names(display_data) <- gsub("_", " ", names(display_data))
    names(display_data) <- tools::toTitleCase(names(display_data))

    datatable(
      display_data,
      rownames = FALSE,
      filter = "top",
      extensions = "FixedHeader",
      class = "stripe hover compact",
      options = list(
        pageLength = 20,
        lengthMenu = c(10, 20, 50, 100),
        deferRender = TRUE,
        processing = TRUE,
        scrollX = TRUE,
        fixedHeader = TRUE
      )
    )
  }, server = TRUE)

  output$viz_corr_dist <- renderPlotly({
    dat <- db_filtered()

    validate(
      need(nrow(dat) > 0, "No rows match the current database filters.")
    )

    plot_ly(
      dat,
      x = ~abs(Correlation),
      color = ~corr_type,
      type = "histogram",
      xbins = list(start = 0.60, end = 1.00, size = 0.02)
    ) |>
      layout(
        barmode = "stack",
        xaxis = list(title = "|Correlation value|"),
        yaxis = list(title = "Count"),
        legend = list(orientation = "h", y = -0.2),
        margin = list(l = 60, r = 20, t = 20, b = 80)
      ) |>
      config(displaylogo = FALSE, responsive = TRUE)
  })

  output$viz_tissue <- renderPlotly({
    counts <- db_filtered() |>
      count(tissue_pair, sort = TRUE) |>
      arrange(n)

    validate(
      need(nrow(counts) > 0, "No rows match the current database filters.")
    )

    plot_ly(
      counts,
      x = ~n,
      y = ~factor(tissue_pair, levels = tissue_pair),
      type = "bar",
      orientation = "h"
    ) |>
      layout(
        xaxis = list(title = "Count"),
        yaxis = list(title = "", automargin = TRUE),
        margin = list(l = 210, r = 20, t = 20, b = 55)
      ) |>
      config(displaylogo = FALSE, responsive = TRUE)
  })

  output$viz_context <- renderPlotly({
    counts <- db_filtered() |>
      count(context, sort = TRUE) |>
      arrange(n)

    validate(
      need(nrow(counts) > 0, "No rows match the current database filters.")
    )

    plot_ly(
      counts,
      x = ~n,
      y = ~factor(context, levels = context),
      type = "bar",
      orientation = "h"
    ) |>
      layout(
        xaxis = list(title = "Count"),
        yaxis = list(title = "", automargin = TRUE),
        margin = list(l = 180, r = 20, t = 20, b = 55)
      ) |>
      config(displaylogo = FALSE, responsive = TRUE)
  })

# output$viz_island <- renderPlotly({
#   counts <- db_filtered() |>
#     count(island_relation, sort = TRUE) |>
#     arrange(n)

#   validate(
#     need(nrow(counts) > 0, "No rows match the current database filters.")
#   )

#   plot_ly(
#     counts,
#     x = ~n,
#     y = ~factor(island_relation, levels = island_relation),
#     type = "bar",
#     orientation = "h"
#   ) |>
#     layout(
#       xaxis = list(title = "Count"),
#       yaxis = list(title = "", automargin = TRUE),
#       margin = list(l = 160, r = 20, t = 20, b = 55)
#     ) |>
#     config(displaylogo = FALSE, responsive = TRUE)
# })
# 
  output$viz_island <- renderPlotly({
    counts <- db_filtered() |>
      count(island_relation, sort = TRUE) |>
      arrange(n)
    
    validate(
      need(nrow(counts) > 0, "No rows match the current database filters.")
    )
    
    plot_ly(
      counts,
      x = ~n,
      y = ~factor(island_relation, levels = island_relation),
      type = "bar",
      orientation = "h"
    ) |>
      layout(
        xaxis = list(title = "Count"),
        yaxis = list(title = "", automargin = TRUE),
        margin = list(l = 160, r = 20, t = 20, b = 55)
      ) |>
      config(displaylogo = FALSE, responsive = TRUE)
  })
  
  # -------------------------------------------------------------------------
  # Glass-brain plot
  # -------------------------------------------------------------------------
  
 # output$glass_brain_plot <- renderPlotly({
 #   
 #   locations <- glass_brain_locations
 #   
 #   if (
 #     !is.null(input$glass_brain_region) &&
 #     input$glass_brain_region != "All regions"
 #   ) {
 #     locations <- locations |>
 #       filter(region_code == input$glass_brain_region)
 #   }
 #   
 #   validate(
 #     need(
 #       nrow(locations) > 0,
 #       "No brain locations match the selected filters."
 #     )
 #   )
 #   
 #   marker_text <- paste0(
 #     "<b>", locations$region_name, "</b>",
 #     "<br>MNI x: ", locations$x,
 #     "<br>MNI y: ", locations$y,
 #     "<br>MNI z: ", locations$z
 #   )
 #   
 #   marker_labels <- if (isTRUE(input$glass_show_labels)) {
 #     locations$region_code
 #   } else {
 #     rep("", nrow(locations))
 #   }
 #   
 #   brain <- ggseg3d(
 #     atlas = aseg(),
 #     na_alpha = 0
 #   ) |>
 #     add_glassbrain(
 #       hemisphere = c("left", "right"),
 #       colour = "#B9C4C7",
 #       opacity = input$glass_brain_opacity
 #     )
 #   
 #   brain |>
 #     add_trace(
 #       data = locations,
 #       x = ~x,
 #       y = ~y,
 #       z = ~z,
 #       type = "scatter3d",
 #       
 #       mode = if (isTRUE(input$glass_show_labels)) {
 #         "markers+text"
 #       } else {
 #         "markers"
 #       },
 #       
 #       text = marker_labels,
 #       textposition = "top center",
 #       hovertext = marker_text,
 #       hoverinfo = "text",
 #       
 #       marker = list(
 #         size = input$glass_point_size,
 #         color = "#078894",
 #         opacity = 0.95,
 #         line = list(
 #           color = "#FFFFFF",
 #           width = 2
 #         )
 #       ),
 #       
 #       name = "Sampling locations",
 #       inherit = FALSE
 #     ) |>
 #     layout(
 #       showlegend = FALSE,
 #       
 #       scene = list(
 #         xaxis = list(
 #           visible = FALSE,
 #           title = ""
 #         ),
 #         
 #         yaxis = list(
 #           visible = FALSE,
 #           title = ""
 #         ),
 #         
 #         zaxis = list(
 #           visible = FALSE,
 #           title = ""
 #         ),
 #         
 #         aspectmode = "data",
 #         
 #         camera = list(
 #           eye = list(
 #             x = 1.55,
 #             y = 0.10,
 #             z = 0.35
 #           )
 #         )
 #       ),
 #       
 #       margin = list(
 #         l = 0,
 #         r = 0,
 #         t = 0,
 #         b = 0
 #       ),
 #       
 #       paper_bgcolor = "rgba(0,0,0,0)",
 #       plot_bgcolor = "rgba(0,0,0,0)"
 #     ) |>
 #     config(
 #       responsive = TRUE,
 #       displaylogo = FALSE
 #     )
 # })
  
  output$glass_brain_plot <- renderGgseg3d({
    
    ggseg3d(
      atlas = aseg(),
      hemisphere = c("left", "right"),
      na_alpha = 0
    ) |>
      add_glassbrain(
        hemisphere = c("left", "right"),
        colour = "#B9C4C7",
        opacity = input$glass_brain_opacity
      ) |>
      set_background("white") |>
      set_legend(FALSE) |>
      pan_camera("right lateral")
  })
  
  selected_brain_region <- reactiveVal(NULL)

  #selected_brain_region <- reactiveVal(NULL)

  query_debounced <- debounce(reactive(input$query), millis = 350)
  correlation_debounced <- debounce(
    reactive(input$minimum_correlation),
    millis = 250
  )

  observeEvent(input$brain_region, {
    req(input$brain_region)
    selected_brain_region(input$brain_region)
  })

  observeEvent(input$clear_brain_region, {
    selected_brain_region(NULL)
    session$sendCustomMessage("clear-brain-selection", list())
  })

  filtered_data <- reactive({
    dat <- cpg_data
    query <- trimws(query_debounced())

    if (nzchar(query)) {
      if (grepl("^cg[0-9]+$", query, ignore.case = TRUE)) {
        dat <- dat |> filter(tolower(CpG_ID) == tolower(query))
      } else {
        dat <- dat |>
          filter(
            !is.na(Gene),
            grepl(
              paste0("(^|;)", toupper(query), "(;|$)"),
              toupper(Gene),
              perl = TRUE
            )
          )
      }
    }

    region <- selected_brain_region()
    if (!is.null(region)) {
      dat <- dat |> filter(Brain_Region_Code == region)
    }

    if (
      !is.null(input$peripheral_tissue) &&
      input$peripheral_tissue != "All tissues"
    ) {
      dat <- dat |> filter(Peripheral_Tissue == input$peripheral_tissue)
    }

    if (isTRUE(input$gene_annotated_only)) {
      dat <- dat |> filter(!is.na(Gene), trimws(Gene) != "")
    }

    dat |>
      filter(
        !is.na(Correlation),
        abs(Correlation) >= correlation_debounced()
      )
  })

  output$selection_summary <- renderUI({
    dat <- filtered_data()
    region <- selected_brain_region()
    region_name <- if (is.null(region)) {
      "All brain regions"
    } else {
      unname(brain_region_lookup[[region]])
    }

    div(
      class = "summary-grid",
      div(
        class = "metric-card",
        span(class = "metric-label", "Brain region"),
        strong(class = "metric-value metric-value-text", region_name)
      ),
      div(
        class = "metric-card",
        span(class = "metric-label", "Matching rows"),
        strong(class = "metric-value", format(nrow(dat), big.mark = ","))
      ),
      div(
        class = "metric-card",
        span(class = "metric-label", "Unique CpGs"),
        strong(
          class = "metric-value",
          format(n_distinct(dat$CpG_ID), big.mark = ",")
        )
      ),
      div(
        class = "metric-card",
        span(class = "metric-label", "Strongest |correlation|"),
        strong(
          class = "metric-value",
          if (nrow(dat) == 0) "—" else sprintf("%.2f", max(abs(dat$Correlation), na.rm = TRUE))
        )
      )
    )
  })

  output$correlation_plot <- renderPlotly({
    dat <- filtered_data()
    validate(need(nrow(dat) > 0, "No CpGs match the current search and filters."))

    plot_data <- dat |>
      transmute(
        CpG = CpG_ID,
        Gene = ifelse(is.na(Gene) | trimws(Gene) == "", "Unannotated", Gene),
        Tissue_Pair = tissue_pair,
        Correlation = Correlation,
        Correlation_Type = corr_type,
        Source = Source,
        Hover = paste0(
          "<b>", CpG_ID, "</b>",
          "<br>Gene: ",
          ifelse(is.na(Gene) | trimws(Gene) == "", "Unannotated", Gene),
          "<br>Tissue pair: ", tissue_pair,
          "<br>Correlation: ", sprintf("%.3f", Correlation),
          "<br>Type: ", corr_type,
          "<br>Source: ", Source
        )
      )

    if (nrow(plot_data) > 4000) {
      plot_data <- plot_data |>
        slice_max(order_by = abs(Correlation), n = 4000, with_ties = FALSE)
    }

    plot_ly(
      data = plot_data,
      x = ~Correlation,
      y = ~Tissue_Pair,
      type = "scattergl",
      mode = "markers",
      text = ~Hover,
      hoverinfo = "text",
      marker = list(size = 6, opacity = 0.65)
    ) |>
      layout(
        autosize = TRUE,
        hovermode = "closest",
        showlegend = FALSE,
        xaxis = list(
          title = "Cross-tissue correlation",
          range = c(-1, 1),
          tickvals = c(-1, -0.5, 0, 0.5, 1),
          zeroline = TRUE
        ),
        yaxis = list(
          title = "",
          automargin = TRUE,
          categoryorder = "array",
          categoryarray = rev(c(
            "Blood vs Prefrontal cortex",
            "Blood vs Superior temporal gyrus",
            "Blood vs Entorhinal cortex",
            "Blood vs Cerebellum",
            "Buccal vs Prefrontal cortex"
          ))
        ),
        margin = list(l = 210, r = 30, t = 12, b = 60)
      ) |>
      config(
        responsive = TRUE,
        displaylogo = FALSE,
        modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d")
      )
  })

  output$results_table <- renderDT({
    dat <- filtered_data()
    validate(need(nrow(dat) > 0, "No records match the current search and filters."))

    table_data <- dat |>
      transmute(
        CpG = CpG_ID,
        Chromosome = chr,
        Position = pos,
        Gene = ifelse(is.na(Gene) | trimws(Gene) == "", "—", Gene),
        Context = context,
        `CpG island relation` = island_relation,
        `Tissue pair` = tissue_pair,
        `Correlation type` = corr_type,
        Correlation = Correlation,
        Source = Source,
        Tool = database_tool,
        `Source link` = source_url
      )

    datatable(
      table_data,
      rownames = FALSE,
      filter = "top",
      extensions = c("Buttons", "FixedHeader"),
      class = "stripe hover compact",
      options = list(
        pageLength = 20,
        lengthMenu = c(10, 20, 50, 100),
        deferRender = TRUE,
        processing = TRUE,
        scrollX = TRUE,
        fixedHeader = TRUE,
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel")
      )
    ) |>
      formatRound(columns = "Correlation", digits = 3)
  }, server = TRUE)
}

shinyApp(ui = ui, server = server)
