#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
})

#settings
input_file  <- "X:/imb-beligr2/MS/outNPLOC4ganz/alphafold_iptm_ranking.csv"
output_dir  <- "X:/imb-beligr2/MS/outNPLOC4ganz"
iptm_cutoff <- 0.5

#read alphafold_iptm_ranking.csv
first_line <- readLines(input_file, n = 1)
sep <- if (str_count(first_line, ";") > str_count(first_line, ",")) ";" else ","

iptm <- read_delim(input_file, delim = sep, show_col_types = FALSE)

#assign and filter by "NPLOC4_full"
iptm_parsed <- iptm %>%
  mutate(
    dir_base = basename(directory),
    OTUD5_peptide_index = as.integer(
      coalesce(
        str_match(dir_base, "OTUD5_peptid([0-9]+)")[, 2],
        str_match(dir_base, "OTUD5_peptide([0-9]+)")[, 2]
      )
    ),
    NPLOC4_full = str_detect(dir_base, "NPLOC4_full")
  ) %>%
  filter(!is.na(OTUD5_peptide_index), NPLOC4_full)

#output
summary_df <- iptm_parsed %>%
  group_by(OTUD5_peptide_index) %>%
  summarise(
    n_models = n(),
    max_iptm  = max(iptm, na.rm = TRUE),
    mean_iptm = mean(iptm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(max_iptm >= iptm_cutoff) %>%
  arrange(desc(max_iptm))

out_file <- file.path(output_dir, "summary_otudPeptide_NPLOC4full.csv")

#output
write.table(
  summary_df,
  file      = out_file,
  sep       = ";",
  dec       = ".",
  row.names = FALSE,
  quote     = FALSE
)
