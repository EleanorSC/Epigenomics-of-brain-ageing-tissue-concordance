
# Prototype TRACE-CpG Shiny App Prototype 31/08/26

Interactive web application for navigating the **TRACE-CpG**: **T**issue **R**eference-database for **A**nalysing **C**oncordant-**T**issue **E**pigenetics

The app provides a searchable UI to CpG-level DNA methylation (DNAm) correlations  derived from matched human brain and peripheral tissue samples.

## What the prototype app does

The app allows users to:

- search for individual CpG sites or genes;
- filter results by brain region and peripheral tissue;
- filter by correlation magnitude and correlation metric;
- filter by chromosome, genomic context and CpG-island relation;
- browse CpG-level concordance estimates in an interactive table;
- view gene-level summaries;
- visualise correlation distributions and database characteristics;
- explore the brain regions represented in the database; and
- download filtered CpG-level results as CSV files.

Each row in the CpG-level database represents a **CpG × brain–peripheral tissue comparison**. A single CpG may therefore occur in multiple rows where concordance has been estimated across different brain regions or tissue pairs.

## Input data

The application currently reads two primary files:

- `data/cpg_brain_peripheral_correlation_database.csv` — harmonised CpG-level concordance database
- `data/gene_level_summary.csv` — gene-level summary derived from the CpG database

## Application structure / tabs

The current prototype contains the following main sections:

- **Search** — search by CpG or gene and interactively explore brain–peripheral concordance.
- **CpG Database** — filter, browse and download CpG-level observations.
- **Gene Summary** — browse gene-level summaries.
- **Visualise** — interactively visualise correlation distributions, tissue pairs, genomic context and CpG-island relationships.
- **Glass Brain** — visualise the brain regions represented in the resource.
- **Access Guide** — information about the database, source resources, interpretation and limitations.

## Running the app

The application is written in R using Shiny.

Required packages include:

```r
library(sf)
library(ggseg3d)
library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(plotly)
