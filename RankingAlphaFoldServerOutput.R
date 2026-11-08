#!/usr/bin/env Rscript

#Rank AlphaFold Server outputs using *summary_confidences*.json.

#Paackages installed?
suppressWarnings({
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Bitte 'jsonlite' installieren.")
  if (!requireNamespace("dplyr",   quietly = TRUE)) stop("Bitte 'dplyr' installieren.")
  if (!requireNamespace("readr",   quietly = TRUE)) stop("Bitte 'readr' installieren.")
  if (!requireNamespace("purrr",   quietly = TRUE)) stop("Bitte 'purrr' installieren.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Bitte 'stringr' installieren.")
})

#load libraries
library(jsonlite)
library(dplyr)
library(readr)
library(purrr)
library(stringr)

#define paths manually!
base_dir   <- "C:/Users/morit/OneDrive - JGU/Desktop/Arbeitsgruppe_Beli/R_Beli_Lab/Results/AlphaFoldServerResults/Erster Ansatz"
output_dir <- "C:/Users/morit/OneDrive - JGU/Desktop/Arbeitsgruppe_Beli/R_Beli_Lab/Results/AlphaFoldServerResults/Erster Ansatz"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#read json/text
`%||%` <- function(a, b) if (!is.null(a)) a else b
num     <- function(x) suppressWarnings(as.numeric(x))

read_json_safely <- function(path) tryCatch(fromJSON(path, simplifyVector = TRUE), error = function(e) NULL)
read_text_safely <- function(path) tryCatch(paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n"), error = function(e) "")

#searches for any "keys"
pluck_any <- function(x, keys) {
  if (is.null(x)) return(NULL)
  for (k in keys) {
    v <- tryCatch(x[[k]], error = function(e) NULL)
    if (!is.null(v)) return(v)
  }
  if (is.list(x)) for (el in x) {
    v <- pluck_any(el, keys)
    if (!is.null(v)) return(v)
  }
  NULL
}

#raw text if not json
extract_number_regex <- function(text, key_regex) {
  if (!nzchar(text)) return(NA_real_)
  pattern <- paste0('(?:"|\\b)', key_regex, '(?:"|\\b)\\s*:\\s*(-?\\d+(?:\\.\\d+)?(?:[eE][+-]?\\d+)?)')
  m <- regexpr(pattern, text, perl = TRUE)
  if (m[1] == -1) return(NA_real_)
  val <- sub("^.*:\\s*", "", regmatches(text, m))
  num(val)
}

#extract iptm from json
extract_conf_record_from_file <- function(f) {
  j <- read_json_safely(f)
  
#search from json
  rs  <- pluck_any(j, c("ranking_score", "rankingScore", "ranking.score"))
  ip  <- pluck_any(j, c("iptm", "iPTM"))
  pt  <- pluck_any(j, c("ptm"))
  mpl <- pluck_any(j, c("mean_plddt", "meanPlddt", "mean.plddt"))
  
#raw text if needed
  if (is.null(rs)) {
    txt <- read_text_safely(f)
    rs  <- extract_number_regex(txt, "(ranking[_\\.]?score|rankingScore)")
#Na -> return
    if (is.na(rs)) return(tibble())  
    if (is.null(ip))  ip  <- extract_number_regex(txt, "i?ptm")
    if (is.null(pt))  pt  <- extract_number_regex(txt, "ptm")
    if (is.null(mpl)) mpl <- extract_number_regex(txt, "(mean[_\\.]?plddt|meanPlddt)")
  }
  
  tibble(
    model        = tools::file_path_sans_ext(basename(f)),
    ranking_score= num(rs),
    iptm         = num(ip),
    ptm          = num(pt),
    mean_plddt   = num(mpl),
    directory    = normalizePath(dirname(f)),
    file         = basename(f),
    source       = paste0("file:", normalizePath(f))
  )
}

#Find files
all_json <- list.files(base_dir, pattern = "\\.json$", full.names = TRUE, recursive = TRUE)
summary_files <- all_json[grepl("summary_confidences", basename(all_json), ignore.case = TRUE, perl = FALSE)]

message("Gefundene Dateien: ", length(summary_files), " JSONs; 0 ZIPs (nur *summary_confidences*.json)")

if (!length(summary_files)) {
  message("Keine *summary_confidences*.json gefunden.")
  return(invisible(NULL))
}

#Read in files, collect iptm
df_all <- purrr::map_dfr(summary_files, function(f) {
  res <- extract_conf_record_from_file(f)
  if (nrow(res)) res else tibble()
})

message("nr of ranking_score files", nrow(df_all))

if (!nrow(df_all)) {
  message("did not find summary_confidences.json with ranking_score in table")
  return(invisible(NULL))
}

#cleans model names
clean_model <- function(x) {
  x %>%
    str_replace("^result[_-]*", "") %>%
    str_replace("^ranking[_-]*", "") %>%
    str_replace("^model[_-]*", "") %>%
    str_replace("\\.pkl$", "") %>%
    str_squish()
}

#create ranking
ranked <- df_all %>%
  mutate(model = clean_model(as.character(model))) %>%
  filter(!is.na(ranking_score)) %>%
  arrange(desc(ranking_score)) %>%
  mutate(rank = row_number()) %>%
  select(rank, ranking_score, iptm, ptm, mean_plddt, model, directory, source, file)

#cancels if table is empty
if (!nrow(ranked)) {
  message("table empty after filtering")
  return(invisible(NULL))
}

#saves ranking in file
out_csv <- file.path(output_dir, "alphafoldserver_ranking_by_ranking_score.csv")
readr::write_csv(ranked, out_csv, na = "")

cat(sprintf("Rangliste (nach ranking_score) gespeichert unter: %s\n", out_csv))
print(utils::head(ranked, 20), row.names = FALSE)

invisible(ranked)
