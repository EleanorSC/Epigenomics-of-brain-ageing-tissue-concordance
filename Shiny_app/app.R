#!/usr/bin/env Rscript
# ============================================================================
# CpG Brain-Peripheral Tissue Correlation Explorer
#
# Interactive replacement for CpG_Brain_Peripheral_Correlation_Database.xlsx.
# Reads the same two source CSVs used to build that workbook and reproduces
# its five sheets (README, CpG_Database, Gene_Summary, Additional_Curated_
# Examples, Database_Access_Guide) as linked, filterable, searchable tabs,
# plus a new Visualize tab with interactive plots (Zissou1 palette).
#
# Data provenance:
#   - Blood-brain rows: Hannon E, Lunnon K, Schalkwyk L, Mill J. Epigenetics
#     2015;10(11):1024-1032. doi:10.1080/15592294.2015.1100786 (PMC4844197)
#   - Buccal-brain rows: Sommerer Y, Ohlei O, Dobricic V, et al. Clinical
#     Epigenetics 2022;14:118. doi:10.1186/s13148-022-01357-w
#
# Run locally:  shiny::runApp("app.R")
# ============================================================================
setwd("/Users/eleanorconole/repos/CpG_Concordance")

library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(plotly)
library(scales)

# ---------------------------------------------------------------------------
# Palette (Wes Anderson "Zissou1", matching the static figures already
# delivered for this database)
# ---------------------------------------------------------------------------
ZISSOU1 <- c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")
ACCENT  <- "#01696F"  # Nexus Hydra Teal, used for chrome/UI

# ---------------------------------------------------------------------------
# Data loading — same two source CSVs used by build_cpg_database.py
# ---------------------------------------------------------------------------
DATA_DIR <- if (dir.exists("data")) "data" else "."

clean_gene <- function(raw_gene) {
  # Illumina manifest gene fields often repeat the same symbol per-probe
  # transcript (e.g. "SDK1;SDK1") or list several distinct genes. Dedupe
  # while preserving order, join with "; ". Vectorized over a character vector.
  vapply(raw_gene, function(x) {
    if (is.na(x) || !nzchar(trimws(x))) return("")
    parts <- trimws(strsplit(x, ";")[[1]])
    parts <- parts[nzchar(parts)]
    paste(unique(parts), collapse = "; ")
  }, character(1), USE.NAMES = FALSE)
}

raw <- read.csv(file.path(DATA_DIR, "cpg_brain_peripheral_correlation_database.csv"),
                 stringsAsFactors = FALSE, colClasses = "character")

cpg_all <- raw %>%
  mutate(
    chr_num     = chr,
    chr_label   = ifelse(nzchar(trimws(chr)), paste0("chr", chr), "n.a."),
    pos         = suppressWarnings(as.integer(pos)),
    gene_clean  = clean_gene(gene),
    context     = ifelse(nzchar(trimws(context)), context, "n.a."),
    island_relation = ifelse(nzchar(trimws(island_relation)), island_relation, "n.a."),
    corr_value  = suppressWarnings(as.numeric(corr_value)),
    has_gene    = nzchar(gene_clean)
  )

# CpG_Database sheet equivalent: gene-annotated rows only (2,937 of 9,127),
# exactly mirroring build_cpg_database.py so this app is a faithful
# replacement for that workbook's primary filterable sheet.
cpg_db <- cpg_all %>%
  filter(has_gene) %>%
  transmute(
    cpg_id   = cpg,
    chr      = chr_label,
    chr_num  = suppressWarnings(as.integer(chr_num)),
    position = pos,
    gene     = gene_clean,
    context  = context,
    island_relation = island_relation,
    tissue_pair = tissue_pair,
    corr_type   = corr_type,
    corr_value  = corr_value,
    database_tool = database_tool,
    source      = source,
    source_url  = source_url
  )

gene_summary <- read.csv(file.path(DATA_DIR, "gene_level_summary.csv"),
                          stringsAsFactors = FALSE) %>%
  transmute(
    gene = gene,
    chr  = ifelse(nzchar(trimws(as.character(chr))), paste0("chr", chr), "n.a."),
    n_qualifying_cpgs = as.integer(n_qualifying_cpgs),
    max_abs_corr = as.numeric(max_abs_corr),
    tissue_pairs = tissue_pairs,
    contexts = contexts
  )

# Expanded single-gene lookup for the searchable selectize (a gene column may
# list several genes separated by "; ") — used to drive the gene filter.
gene_choices <- cpg_db %>%
  pull(gene) %>%
  strsplit("; ") %>%
  unlist() %>%
  unique() %>%
  sort()

tissue_pair_choices <- sort(unique(cpg_db$tissue_pair))
context_choices      <- sort(unique(cpg_db$context))
island_choices        <- sort(unique(cpg_db$island_relation))
chr_choices <- cpg_db %>%
  distinct(chr, chr_num) %>%
  arrange(chr_num) %>%
  pull(chr)

# ---------------------------------------------------------------------------
# Additional_Curated_Examples sheet equivalent — 12 individually verified
# CpGs from secondary papers that queried BECon / the Essex tool directly,
# not part of the two bulk-downloaded supplementary tables above.
# ---------------------------------------------------------------------------
curated_examples <- tibble::tribble(
  ~cpg_id, ~chr_position, ~gene, ~context, ~tissue_pair, ~correlation, ~source, ~source_url, ~notes,
  "cg07093428", "chr11 (exact position n.a. in source)", "LDHC", "n.a. in source",
  "Blood vs Brain (mean Brodmann Areas 7, 10, 20)", "Spearman r = 0.62 (via BECon Tool)",
  "Dahrendorff et al. 2025, BMC Medical Genomics 18:181, Table 3",
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC12613596/table/Tab3/",
  "Identified in a TMS treatment-response DMR analysis; cross-referenced against BECon",

  "cg11821245", "chr11 (exact position n.a. in source)", "LDHC", "n.a. in source",
  "Blood vs Brain (mean BA7/10/20)", "Spearman r = 0.62 (via BECon Tool)",
  "Dahrendorff et al. 2025, BMC Medical Genomics 18:181, Table 3",
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC12613596/table/Tab3/", "",

  "cg00124993", "chr5 (exact position n.a. in source)", "MIR886", "n.a. in source",
  "Blood vs Brain (mean BA7/10/20)", "Spearman r = 0.72 (via BECon Tool)",
  "Dahrendorff et al. 2025, BMC Medical Genomics 18:181, Table 3",
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC12613596/table/Tab3/",
  "Strongest of a 5-CpG MIR886 cluster reported in this DMR analysis",

  "cg06536614", "chr5 (exact position n.a. in source)", "MIR886", "n.a. in source",
  "Blood vs Brain (mean BA7/10/20)", "Spearman r = 0.69 (via BECon Tool)",
  "Dahrendorff et al. 2025, BMC Medical Genomics 18:181, Table 3",
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC12613596/table/Tab3/", "",

  "cg08658272", "chr5 (exact position n.a. in source)", "MIR886", "n.a. in source",
  "Blood vs Brain (mean BA7/10/20)", "Spearman r = 0.66 (via BECon Tool)",
  "Dahrendorff et al. 2025, BMC Medical Genomics 18:181, Table 3",
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC12613596/table/Tab3/", "",

  "cg18894933", "chr5 (exact position n.a. in source)", "MIR886", "n.a. in source",
  "Blood vs Brain (mean BA7/10/20)", "Spearman r = 0.63 (via BECon Tool)",
  "Dahrendorff et al. 2025, BMC Medical Genomics 18:181, Table 3",
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC12613596/table/Tab3/", "",

  "cg12028246", "chr5 (exact position n.a. in source)", "MIR886", "n.a. in source",
  "Blood vs Brain (mean BA7/10/20)", "Spearman r = 0.61 (via BECon Tool)",
  "Dahrendorff et al. 2025, BMC Medical Genomics 18:181, Table 3",
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC12613596/table/Tab3/", "",

  "cg23576855", "chr5:373,299 (hg19)", "AHRR", "n.a. in source",
  "Blood vs Prefrontal Cortex (PFC)", "r = 0.91 (via Essex Blood-Brain Tool)",
  "PGC PTSD Epigenetics Workgroup et al. 2024, Genome Medicine 16:147",
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC11275670/",
  "Well-known smoking-associated AHRR locus; queried live via the Essex tool, not present in the bulk-downloaded Supplementary Tables 6/7 subset",

  "cg05656210", "chr5:141,660,565 (hg19)", "Intergenic", "Intergenic",
  "Blood vs PFC/EC/STG/CER (all 4 regions)", "r >= 0.93 for all 4 regions",
  "Snijders et al. 2020, Clinical Epigenetics 12:11",
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC6958602/",
  "PTSD-associated CpG; per-region exact values in a supplementary table not independently retrievable",

  "cg12169700", "chr7:1,923,695 (hg19)", "MAD1L1", "n.a. in source",
  "Blood vs PFC/EC/STG/CER (all 4 regions)", "r >= 0.93 for all 4 regions",
  "Snijders et al. 2020, Clinical Epigenetics 12:11",
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC6958602/", "",

  "cg20756026", "chr17:80,394,529 (hg19)", "HEXDC", "n.a. in source",
  "Blood vs PFC/EC/STG/CER (all 4 regions)", "r >= 0.93 for all 4 regions",
  "Snijders et al. 2020, Clinical Epigenetics 12:11",
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC6958602/", "",

  "cg04987734", "chr14:103,415,873 (hg19)", "CDC42BPB", "n.a. in source",
  "Blood vs Brodmann Area 7", "r = 0.81 (via BECon Tool)",
  "Jovanova et al. 2018, JAMA Psychiatry",
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC6142917/",
  "From a multiethnic depression-methylation meta-analysis"
)

# ---------------------------------------------------------------------------
# Database_Access_Guide sheet equivalent — direct links to the 5 primary
# brain-peripheral tissue methylation correlation resources.
# ---------------------------------------------------------------------------
access_guide <- tibble::tribble(
  ~tool, ~url, ~tissues, ~sample_size, ~how_to_query, ~citation,
  "IMAGE-CpG", "http://han-lab.org/methylation/default/imageCpG",
  "Live brain (surgical resection gray matter) vs blood, saliva, AND buccal — the only bulk-queryable resource covering all three peripheral tissues plus brain in the same living individuals",
  "27 epilepsy patients; Illumina 450K (12 subjects: brain/blood/saliva) and EPIC (21 subjects: brain/blood/saliva/buccal)",
  "Enter a gene name, CpG ID, or genomic coordinate range into the web search tool to view Pearson correlation coefficients between brain and each peripheral tissue at individual CpGs",
  "Braun PR et al. 2019, Transl Psychiatry 9:47, doi:10.1038/s41398-019-0376-y",

  "BECon (Blood-Brain Epigenetic Concordance)", "https://redgar598.shinyapps.io/BECon/",
  "Postmortem brain (Brodmann Areas 7, 10, 20) vs whole blood",
  "16 individuals; Illumina 450K array",
  "Search by gene symbol or CpG ID in the Shiny app to view Spearman correlation coefficients for each brain region, plus CpG variability and cell-composition-effect metrics",
  "Edgar RD et al. 2017, Transl Psychiatry 7:e1187, doi:10.1038/tp.2017.171",

  "Blood Brain DNA Methylation Comparison Tool", "https://epigenetics.essex.ac.uk/bloodbrain/",
  "Postmortem brain (prefrontal cortex, entorhinal cortex, superior temporal gyrus, cerebellum) vs whole blood",
  "71-75 individuals; Illumina 450K array",
  "Query by probe ID (cg-number) or gene name via the web interface for correlation coefficients (r) and scatterplots per brain region. This tool's underlying Supplementary Tables were used to build the CpG Database / Gene Summary tabs in this app",
  "Hannon E et al. 2015, Epigenetics 10(11):1024-1032, doi:10.1080/15592294.2015.1100786",

  "Buccal-Brain Correlation Map (MADRC)", "http://www.liga.uni-luebeck.de/buccal_brain_correlation_results/",
  "Postmortem prefrontal cortex vs buccal (oral) epithelial cells",
  "120 paired samples; Illumina EPIC array",
  "Look up a CpG or gene of interest for correlation statistics and location plots. This resource's Supplementary Tables were also used to build the CpG Database / Gene Summary tabs in this app",
  "Sommerer Y et al. 2022, Clinical Epigenetics 14:118, doi:10.1186/s13148-022-01357-w",

  "AMAZE-CpG", "https://snishit-amaze-cpg.web.app/",
  "Living brain tissue vs blood, saliva, AND buccal epithelial cells in an independent Japanese/Asian-ancestry cohort",
  "19 Japanese subjects; Illumina EPIC array",
  "Query individual CpGs or genes to retrieve Spearman rho for each brain-peripheral tissue pair; useful as an independent-population replication check, especially for saliva-brain correlations",
  "Nishitani S et al. 2023, Transl Psychiatry 13:72, doi:10.1038/s41398-023-02370-0"
)

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
app_theme <- bs_theme(
  version = 5,
  bg = "#F7F6F2", fg = "#28251D",
  primary = ACCENT, secondary = "#7A7974",
  base_font = font_google("Inter"),
  heading_font = font_google("Inter")
)

ui <- navbarPage(
  title = "CpG Brain\u2013Peripheral Tissue Correlation Explorer",
  theme = app_theme,
  id = "nav",

  # ---- About tab ---------------------------------------------------------
  tabPanel("About", icon = icon("info-circle"),
    fluidRow(column(10, offset = 1,
      div(style = "padding: 24px 0;",
        h2("Human CpG Sites with Validated Brain\u2013Peripheral Tissue Methylation Correlations (|r| or |rho| > 0.6)"),
        p(em("A search-friendly database for identifying surrogate peripheral-tissue markers of brain DNA methylation")),
        hr(),
        h4("What this app contains"),
        p("This app compiles individual CpG sites (Illumina 450K/EPIC probe IDs) for which a validated correlation ",
          "coefficient (|r| or |rho| > 0.6) has been reported between DNA methylation in brain gray matter and a ",
          "peripheral surrogate tissue (whole blood or buccal epithelium), using paired/matched samples from the ",
          "same individuals. It is an interactive replacement for the equivalent Excel workbook, with the same ",
          "underlying data."),
        h4("Tabs in this app"),
        tags$ul(
          tags$li(strong("CpG Database"), sprintf(" \u2014 the primary filterable table: %s gene-annotated CpG \u00d7 tissue-pair rows (chromosome, position, gene, promoter/intragenic context, CpG-island relation, tissue pair, correlation value, source).", comma(nrow(cpg_db)))),
          tags$li(strong("Gene Summary"), sprintf(" \u2014 one row per gene (%s genes), showing how many qualifying CpGs exist for that gene, the strongest correlation found, and which tissue pairs/contexts are covered. Click a row here to jump straight to that gene's CpGs in the CpG Database tab.", comma(nrow(gene_summary)))),
          tags$li(strong("Visualize"), " \u2014 interactive charts (Zissou1 palette) that update live as you apply filters on the CpG Database tab."),
          tags$li(strong("Curated Examples"), " \u2014 12 individually verified high-correlation CpGs from other primary cross-tissue methylation resources (BECon, IMAGE-CpG-related literature, Essex/Exeter tool as cited in secondary EWAS papers) that could not be bulk-downloaded but are documented with exact values in peer-reviewed text, including SALIVA-brain examples not covered by the two bulk-downloaded sources."),
          tags$li(strong("Access Guide"), " \u2014 direct links and usage notes for all five primary brain-peripheral tissue methylation correlation resources, so you can look up any additional target gene not already listed here.")
        ),
        h4("Data provenance (CpG Database + Gene Summary tabs)"),
        tags$ul(
          tags$li(HTML("<strong>Blood-brain rows</strong> (tissue pair = \u201cBlood vs [PFC/EC/STG/CER] (brain)\u201d): extracted from Supplementary Tables 2, 6 and 7 of Hannon E, Lunnon K, Schalkwyk L, Mill J. \u201cInterindividual methylomic variation across blood, cortex, and cerebellum.\u201d <em>Epigenetics</em> 2015;10(11):1024-1032. doi:10.1080/15592294.2015.1100786. <a href='https://pmc.ncbi.nlm.nih.gov/articles/PMC4844197/' target='_blank'>PMC4844197</a>. Pearson r between blood methylation and each of 4 postmortem brain regions across 71-75 matched individuals, Illumina 450K array.")),
          tags$li(HTML("<strong>Buccal-brain rows</strong> (tissue pair = \u201cBuccal vs Prefrontal Cortex (brain)\u201d): extracted from Supplementary Tables 2 and 6 of Sommerer Y, Ohlei O, Dobricic V, et al. \u201cA correlation map of genome-wide DNA methylation patterns between paired human brain and buccal samples.\u201d <em>Clinical Epigenetics</em> 2022;14:118. <a href='https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-022-01357-w' target='_blank'>doi:10.1186/s13148-022-01357-w</a>. Spearman rho between buccal and prefrontal cortex methylation across 120 matched postmortem samples, Illumina EPIC array.")),
          tags$li("Genomic context classification: derived from the Illumina UCSC_RefGene_Group manifest annotation. \u201cPromoter-associated\u201d = TSS200, TSS1500, 5'UTR, or 1stExon. \u201cGene body/intragenic\u201d = Body or 3'UTR. \u201cUnannotated/Intergenic\u201d = no RefGene annotation in the manifest.")
        ),
        h4("Important methodological note"),
        p("A meaningful fraction of very high (r > 0.9) blood-brain correlations reflect strong genetic (SNP/mQTL) ",
          "control of methylation at that CpG rather than a brain-specific biological process \u2014 this is a known, ",
          "well-documented feature of cross-tissue methylation correlation and does not reduce the validity of these ",
          "CpGs as statistically reliable peripheral surrogates, per Hannon et al. 2015 and Braun et al. 2019."),
        h4("How to use this app"),
        tags$ol(
          tags$li("Go to Gene Summary and search the gene column to check if your gene(s) of interest already have a validated surrogate CpG."),
          tags$li("Click that row (or switch to CpG Database and filter by gene name) to see the exact CpG ID(s), tissue pair(s), and correlation value(s) to select from."),
          tags$li("If your gene is not listed, or you need saliva-brain or additional buccal/blood coverage, check Curated Examples, then consult the live tools in Access Guide directly (IMAGE-CpG and AMAZE-CpG in particular include saliva as a tissue)."),
          tags$li("Always cite the original source study/database (Source / Database-Tool columns in CpG Database) when using a value in your own research, not this compilation app.")
        )
      )
    ))
  ),

  # ---- CpG Database tab ---------------------------------------------------
  tabPanel("CpG Database", icon = icon("table"),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        selectizeInput("f_gene", "Gene", choices = NULL, multiple = TRUE,
                        options = list(placeholder = "Type to search a gene...")),
        selectInput("f_chr", "Chromosome", choices = c("All", chr_choices),
                    selected = "All", multiple = TRUE),
        checkboxGroupInput("f_tissue", "Tissue pair", choices = tissue_pair_choices,
                            selected = tissue_pair_choices),
        radioButtons("f_corr_type", "Correlation type",
                     choices = c("Both" = "Both", "Pearson r" = "Pearson r", "Spearman rho" = "Spearman rho"),
                     selected = "Both"),
        sliderInput("f_corr_range", "|Correlation value|", min = 0.6, max = 1.0,
                    value = c(0.6, 1.0), step = 0.01),
        checkboxGroupInput("f_context", "Genomic context", choices = context_choices,
                            selected = context_choices),
        checkboxGroupInput("f_island", "CpG island relation", choices = island_choices,
                            selected = island_choices),
        actionButton("reset_filters", "Reset filters", icon = icon("rotate-left"), width = "100%")
      ),
      mainPanel(
        width = 9,
        fluidRow(
          column(4, wellPanel(h4(textOutput("kpi_rows")), p("Matching rows"))),
          column(4, wellPanel(h4(textOutput("kpi_cpgs")), p("Unique CpGs"))),
          column(4, wellPanel(h4(textOutput("kpi_genes")), p("Unique genes")))
        ),
        downloadButton("download_cpg", "Download filtered rows (CSV)"),
        br(), br(),
        DTOutput("cpg_table")
      )
    )
  ),

  # ---- Gene Summary tab ----------------------------------------------------
  tabPanel("Gene Summary", icon = icon("dna"),
    fluidRow(column(12,
      p("One row per gene. Click a row to jump to its CpGs in the CpG Database tab."),
      DTOutput("gene_table")
    ))
  ),

  # ---- Visualize tab --------------------------------------------------------
  tabPanel("Visualize", icon = icon("chart-column"),
    p("These charts reflect the current filters applied on the CpG Database tab."),
    fluidRow(
      column(6, plotlyOutput("plot_corr_dist", height = "380px")),
      column(6, plotlyOutput("plot_tissue_bar", height = "380px"))
    ),
    br(),
    fluidRow(
      column(6, plotlyOutput("plot_context_bar", height = "380px")),
      column(6, plotlyOutput("plot_island_bar", height = "380px"))
    )
  ),

  # ---- Curated Examples tab -------------------------------------------------
  tabPanel("Curated Examples", icon = icon("star"),
    p("12 individually verified high-correlation CpGs from secondary papers that queried BECon or the Essex ",
      "Blood-Brain Tool directly. They are not part of the two bulk supplementary tables used for the CpG ",
      "Database / Gene Summary tabs, so are listed separately here. IMAGE-CpG and AMAZE-CpG (which uniquely ",
      "cover saliva-brain and buccal-brain comparisons and an independent Japanese cohort) could not be ",
      "bulk-downloaded or queried programmatically for arbitrary genes \u2014 see the Access Guide tab to query ",
      "them directly for your own target gene set, including saliva-brain correlations not represented ",
      "anywhere else in this app."),
    DTOutput("curated_table")
  ),

  # ---- Access Guide tab -----------------------------------------------------
  tabPanel("Access Guide", icon = icon("link"),
    p("Direct links to the five primary brain-peripheral tissue methylation correlation resources."),
    DTOutput("access_table")
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  updateSelectizeInput(session, "f_gene", choices = gene_choices, server = TRUE)

  observeEvent(input$reset_filters, {
    updateSelectizeInput(session, "f_gene", selected = character(0), choices = gene_choices, server = TRUE)
    updateSelectInput(session, "f_chr", selected = "All")
    updateCheckboxGroupInput(session, "f_tissue", selected = tissue_pair_choices)
    updateRadioButtons(session, "f_corr_type", selected = "Both")
    updateSliderInput(session, "f_corr_range", value = c(0.6, 1.0))
    updateCheckboxGroupInput(session, "f_context", selected = context_choices)
    updateCheckboxGroupInput(session, "f_island", selected = island_choices)
  })

  filtered_cpg <- reactive({
    d <- cpg_db

    if (length(input$f_gene) > 0) {
      pattern <- paste(str_replace_all(input$f_gene, "([\\W])", "\\\\\\1"), collapse = "|")
      d <- d[grepl(pattern, d$gene), ]
    }
    if (!is.null(input$f_chr) && !("All" %in% input$f_chr)) {
      d <- d[d$chr %in% input$f_chr, ]
    }
    if (!is.null(input$f_tissue)) {
      d <- d[d$tissue_pair %in% input$f_tissue, ]
    } else {
      d <- d[0, ]
    }
    if (!is.null(input$f_corr_type) && input$f_corr_type != "Both") {
      d <- d[d$corr_type == input$f_corr_type, ]
    }
    d <- d[abs(d$corr_value) >= input$f_corr_range[1] & abs(d$corr_value) <= input$f_corr_range[2], ]
    if (!is.null(input$f_context)) {
      d <- d[d$context %in% input$f_context, ]
    } else {
      d <- d[0, ]
    }
    if (!is.null(input$f_island)) {
      d <- d[d$island_relation %in% input$f_island, ]
    } else {
      d <- d[0, ]
    }
    d
  })

  output$kpi_rows  <- renderText(comma(nrow(filtered_cpg())))
  output$kpi_cpgs  <- renderText(comma(length(unique(filtered_cpg()$cpg_id))))
  output$kpi_genes <- renderText({
    genes <- filtered_cpg()$gene
    comma(length(unique(unlist(strsplit(genes, "; ")))))
  })

  output$cpg_table <- renderDT({
    d <- filtered_cpg() %>%
      transmute(
        `CpG ID` = cpg_id, Chromosome = chr, `Position (hg19)` = position,
        Gene = gene, `Genomic Context` = context, `Island Relation` = island_relation,
        `Tissue Pair` = tissue_pair, `Correlation Type` = corr_type,
        `Correlation Value` = round(corr_value, 3),
        `Database/Tool` = database_tool, `Source` = source,
        `Source URL` = ifelse(nzchar(source_url),
                               sprintf('<a href="%s" target="_blank">link</a>', trimws(strsplit(source_url, ";")[[1]][1])),
                               "")
      )
    datatable(d, escape = FALSE, rownames = FALSE, filter = "top",
              extensions = "Buttons",
              options = list(pageLength = 15, scrollX = TRUE,
                              dom = "Bfrtip", buttons = c("copy")))
  })

  output$download_cpg <- downloadHandler(
    filename = function() sprintf("cpg_database_filtered_%s.csv", Sys.Date()),
    content = function(file) write.csv(filtered_cpg(), file, row.names = FALSE)
  )

  # ---- Gene Summary table + click-through to CpG Database ----------------
  output$gene_table <- renderDT({
    datatable(gene_summary %>%
                rename(Gene = gene, Chromosome = chr,
                       `# Qualifying CpGs (r>0.6)` = n_qualifying_cpgs,
                       `Max |Correlation|` = max_abs_corr,
                       `Tissue Pairs Covered` = tissue_pairs,
                       `Genomic Contexts Covered` = contexts) %>%
                mutate(`Max |Correlation|` = round(`Max |Correlation|`, 3)) %>%
                arrange(desc(`# Qualifying CpGs (r>0.6)`)),
              rownames = FALSE, selection = "single", filter = "top",
              options = list(pageLength = 15, scrollX = TRUE))
  })

  observeEvent(input$gene_table_rows_selected, {
    sel <- input$gene_table_rows_selected
    if (length(sel) == 1) {
      g <- gene_summary %>% arrange(desc(n_qualifying_cpgs)) %>% slice(sel) %>% pull(gene)
      updateSelectizeInput(session, "f_gene", choices = gene_choices, selected = g, server = TRUE)
      updateTabsetPanel(session, "nav", selected = "CpG Database")
    }
  })

  # ---- Visualize tab: reactive plots on the same filtered data ------------
  # NOTE: these use plotly's native plot_ly() API rather than ggplot2 +
  # ggplotly(), which currently throws "subscript out of bounds" for
  # histogram/bar geoms with the ggplot2 4.x / plotly 4.10.4 combination
  # installed in this environment (verified via isolated reproduction).
  output$plot_corr_dist <- renderPlotly({
    d <- filtered_cpg()
    validate(need(nrow(d) > 0, "No rows match the current filters."))
    plot_ly(d, x = ~abs(corr_value), color = ~corr_type,
            colors = c("Pearson r" = ZISSOU1[1], "Spearman rho" = ZISSOU1[5]),
            type = "histogram",
            xbins = list(start = 0.6, end = 1.0, size = 0.02)) %>%
      layout(title = "Distribution of correlation values",
             xaxis = list(title = "|Correlation value|"),
             yaxis = list(title = "Count"),
             barmode = "stack",
             legend = list(orientation = "h", y = -0.2))
  })

  output$plot_tissue_bar <- renderPlotly({
    d <- filtered_cpg()
    validate(need(nrow(d) > 0, "No rows match the current filters."))
    counts <- d %>% count(tissue_pair, sort = TRUE) %>% arrange(n)
    plot_ly(counts, x = ~n, y = ~factor(tissue_pair, levels = tissue_pair),
            type = "bar", orientation = "h",
            marker = list(color = rep(ZISSOU1, length.out = nrow(counts)))) %>%
      layout(title = "Rows by tissue pair",
             xaxis = list(title = "Count"), yaxis = list(title = ""))
  })

  output$plot_context_bar <- renderPlotly({
    d <- filtered_cpg()
    validate(need(nrow(d) > 0, "No rows match the current filters."))
    counts <- d %>% count(context, sort = TRUE) %>% arrange(n)
    plot_ly(counts, x = ~n, y = ~factor(context, levels = context),
            type = "bar", orientation = "h",
            marker = list(color = rep(ZISSOU1, length.out = nrow(counts)))) %>%
      layout(title = "Rows by genomic context",
             xaxis = list(title = "Count"), yaxis = list(title = ""))
  })

  output$plot_island_bar <- renderPlotly({
    d <- filtered_cpg()
    validate(need(nrow(d) > 0, "No rows match the current filters."))
    counts <- d %>% count(island_relation, sort = TRUE) %>% arrange(n)
    plot_ly(counts, x = ~n, y = ~factor(island_relation, levels = island_relation),
            type = "bar", orientation = "h",
            marker = list(color = rep(ZISSOU1, length.out = nrow(counts)))) %>%
      layout(title = "Rows by CpG island relation",
             xaxis = list(title = "Count"), yaxis = list(title = ""))
  })

  # ---- Curated Examples tab -------------------------------------------------
  output$curated_table <- renderDT({
    d <- curated_examples %>%
      transmute(
        `CpG ID` = cpg_id, `Chromosome:Position` = chr_position, Gene = gene,
        `Genomic Context` = context, `Tissue Pair` = tissue_pair,
        `Correlation (type, value)` = correlation,
        `Source Study` = source,
        `Source URL` = sprintf('<a href="%s" target="_blank">link</a>', source_url),
        Notes = notes
      )
    datatable(d, escape = FALSE, rownames = FALSE,
              options = list(pageLength = 12, scrollX = TRUE, dom = "t"))
  })

  # ---- Access Guide tab -----------------------------------------------------
  output$access_table <- renderDT({
    d <- access_guide %>%
      transmute(
        `Database/Tool` = tool,
        URL = sprintf('<a href="%s" target="_blank">%s</a>', url, url),
        `Tissues Covered` = tissues,
        `Sample Size / Array` = sample_size,
        `How to Query for Your Gene` = how_to_query,
        `Primary Citation` = citation
      )
    datatable(d, escape = FALSE, rownames = FALSE,
              options = list(pageLength = 5, scrollX = TRUE, dom = "t"))
  })
}

shinyApp(ui, server)
