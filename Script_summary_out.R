#load packages
library(dplyr)
library(stringr)
library(readr)

#settings
input_file  <- "Z:/IMB/Beli/alphafold/Moritz/out/alphafold_iptm_ranking_out.csv"
output_dir  <- "Z:/IMB/Beli/alphafold/Moritz/out"
iptm_cutoff <- 0.5


#input, define separator
first_line <- readLines(input_file, n = 1)
sep <- if (str_count(first_line, ";") > str_count(first_line, ",")) ";" else ","

iptm <- read_delim(input_file, delim = sep, show_col_types = FALSE)


#zuordnen (iptm zu directory name)
iptm_parsed <- iptm %>%
  mutate(
    dir_base = basename(directory),
    
#NPLOC4 domain
    NPLOC4_domain_index = as.integer(str_match(dir_base, "NPLOC4_domain([0-9]+)")[, 2]),
    
#OTUD5 domain
    OTUD5_domain_index  = as.integer(str_match(dir_base, "OTUD5_domain([0-9]+)")[, 2]),

#OTUD5 peptid(e)
    OTUD5_peptide_index = as.integer(
      coalesce(
        str_match(dir_base, "OTUD5_peptid([0-9]+)")[, 2],
        str_match(dir_base, "OTUD5_peptide([0-9]+)")[, 2]
      )
    ),
    
#NPLOC4 peptid(e)
    NPLOC4_peptide_index = as.integer(
      coalesce(
        str_match(dir_base, "NPLOC4_peptid([0-9]+)")[, 2],
        str_match(dir_base, "NPLOC4_peptide([0-9]+)")[, 2]
      )
    )
  )


#sort and filter iptm>0.5 based on original fasta name
df1 <- iptm_parsed %>%
  filter(!is.na(NPLOC4_domain_index),
         !is.na(OTUD5_peptide_index)) %>%
  group_by(NPLOC4_domain_index, OTUD5_peptide_index) %>%
  summarise(
    n_models = n(),
    max_iptm = max(iptm, na.rm = TRUE),
    mean_iptm = mean(iptm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(max_iptm >= iptm_cutoff) %>%
  arrange(desc(max_iptm))

write.table(
  df1,
  file = file.path(output_dir, "summary_nplocDomain_otudPeptide.csv"),
  sep = ";",
  dec = ".",
  row.names = FALSE,
  quote = FALSE
)


df2 <- iptm_parsed %>%
  filter(!is.na(OTUD5_domain_index),
         !is.na(NPLOC4_peptide_index)) %>%
  group_by(OTUD5_domain_index, NPLOC4_peptide_index) %>%
  summarise(
    n_models = n(),
    max_iptm = max(iptm, na.rm = TRUE),
    mean_iptm = mean(iptm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(max_iptm >= iptm_cutoff) %>%
  arrange(desc(max_iptm))

write.table(
  df2,
  file = file.path(output_dir, "summary_otudDomain_nplocPeptide.csv"),
  sep = ";",
  dec = ".",
  row.names = FALSE,
  quote = FALSE
)