#!/usr/bin/env Rscript

#to run this script: cmd line: system('Rscript "Analyze_co_IP_hits_with_4_args.R" "X:/imb-beligr2/MS"')

suppressWarnings({
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Bitte 'jsonlite' installieren.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Bitte 'dplyr' installieren.")
  if (!requireNamespace("readr", quietly = TRUE)) stop("Bitte 'readr' installieren.")
  if (!requireNamespace("purrr", quietly = TRUE)) stop("Bitte 'purrr' installieren.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Bitte 'stringr' installieren.")
})

library(jsonlite)
library(dplyr)
library(readr)
library(purrr)
library(stringr)

#read in arguments
args <- commandArgs(trailingOnly = TRUE)

base_dir <- args[1]
csv_output <- args[2]
fasta_dir <- base_dir

cat("Alphafold directory: ", base_dir, "\n")
cat("CSV output file: ", csv_output, "\n")

#read in files
fastas <- list.files(fasta_dir, pattern = "\\.fasta$", full.names = TRUE)

#create output directory if missing
dir.create(dirname(csv_output), recursive = TRUE, showWarnings = FALSE)

#Helper functions
`%||%` <- function(x, y) if (is.null(x)) y else x

read_json_safely <- function(path) {
  tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE), error = function(e) NULL)
}

collect_from_result_json <- function(dir) {
  files <- list.files(dir, pattern = "^result_.*\\.json$", full.names = TRUE, recursive = FALSE)
  if (length(files) == 0) return(tibble())
  purrr::map_dfr(files, function(f) {
    j <- read_json_safely(f)
    if (is.null(j)) return(tibble())
    iptm <- j[["iptm"]]
    ptm  <- j[["ptm"]]
    mpl  <- j[["mean_plddt"]]
    tibble(
      directory  = normalizePath(dirname(f), winslash = "/"),
      file       = basename(f),
      model      = tools::file_path_sans_ext(basename(f)),
      iptm       = suppressWarnings(as.numeric(iptm)),
      ptm        = suppressWarnings(as.numeric(ptm)),
      mean_plddt = suppressWarnings(as.numeric(mpl)),
      source     = "result_json"
    )
  })
}

collect_from_ranking_debug <- function(dir) {
  f <- file.path(dir, "ranking_debug.json")
  if (!file.exists(f)) return(tibble())
  j <- read_json_safely(f)
  if (is.null(j)) return(tibble())
  
  df1 <- try({
    if (!is.null(j$iptm) && is.list(j$iptm)) {
      tibble(
        directory  = normalizePath(dir, winslash = "/"),
        file       = "ranking_debug.json",
        model      = names(j$iptm),
        iptm       = as.numeric(unlist(j$iptm)),
        ptm        = NA_real_,
        mean_plddt = NA_real_,
        source     = "ranking_debug"
      )
    } else tibble()
  }, silent = TRUE)
  
  df2 <- try({
    mc <- j$model_confidences
    if (!is.null(mc) && (is.data.frame(mc) || is.matrix(mc))) {
      cols    <- tolower(colnames(mc))
      iptmcol <- which(cols == "iptm")
      ptmcol  <- which(cols == "ptm")
      mplcol  <- which(cols %in% c("mean_plddt","mean_plddt_0","mean_plddt_1"))
      if (length(iptmcol) || length(ptmcol) || length(mplcol)) {
        tibble(
          directory  = normalizePath(dir, winslash = "/"),
          file       = "ranking_debug.json",
          model      = rownames(mc) %||% as.character(seq_len(nrow(mc))),
          iptm       = if (length(iptmcol)) as.numeric(mc[[iptmcol[1]]]) else NA_real_,
          ptm        = if (length(ptmcol))  as.numeric(mc[[ptmcol[1]]])  else NA_real_,
          mean_plddt = if (length(mplcol))  as.numeric(mc[[mplcol[1]]])  else NA_real_,
          source     = "ranking_debug"
        )
      } else tibble()
    } else tibble()
  }, silent = TRUE)
  
  bind_rows(if (inherits(df1, "try-error")) tibble() else df1,
            if (inherits(df2, "try-error")) tibble() else df2) |>
    distinct(directory, model, .keep_all = TRUE)
}

#Ordner finden
candidate_dirs <- unique(dirname(list.files(
  base_dir,
  pattern = "ranking_debug\\.json|ranked_.*\\.pdb|result_.*\\.json",
  full.names = TRUE, recursive = TRUE
)))

if (length(candidate_dirs) == 0) {
  message("Keine AlphaFold-Output-Dateien gefunden unter: ", base_dir)
  
  invisible(return(NULL))
}

df <- purrr::map_dfr(candidate_dirs, function(d) {
  bind_rows(
    collect_from_result_json(d),
    collect_from_ranking_debug(d)
  )
})

if (nrow(df) == 0 || all(is.na(df$iptm))) {
  message("Es wurden Dateien gefunden, aber keine ipTM-Werte extrahiert.")
  invisible(return(NULL))
}

#create ranking
ranked <- df |>
  mutate(iptm = suppressWarnings(as.numeric(iptm))) |>
  filter(!is.na(iptm)) |>
  arrange(desc(iptm)) |>
  mutate(rank = dplyr::row_number()) |>
  select(rank, iptm, ptm, mean_plddt, model, directory, source, file)

#highscore
highscore <- ranked |>
  group_by(directory) |>                        
  arrange(desc(iptm), .by_group = TRUE) |>       
  slice(1L) |>                                   
  ungroup() |>
  arrange(desc(iptm)) |>
  mutate(global_rank = dplyr::row_number())


#save ranking as csv file
ranking_csv <- file.path(csv_output, "alphafold_iptm_ranking.csv")
readr::write_csv(ranked, ranking_csv, na = "")

cat(sprintf("Rangliste gespeichert unter: %s\n", ranking_csv))
print(utils::head(ranked, 10), row.names = FALSE)

#save highscore as csv file
highscore_csv <- file.path(csv_output, "alphafold_iptm_highscore.csv")
readr::write_csv(highscore, highscore_csv, na = "")

cat(sprintf("Highscore gespeichert unter: %s\n", highscore_csv))
print(utils::head(highscore, 10), row.names = FALSE)