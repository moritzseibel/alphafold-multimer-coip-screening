#load packages
library(jsonlite)
library(dplyr)
library(readr)
library(purrr)
library(stringr)

#define dirs
base_dir   <- "X:/imb-beligr2/MS/outNPLOC4ganz"
output_dir <- "X:/imb-beligr2/MS/outNPLOC4ganz"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#helper functions
`%||%` <- function(x, y) if (is.null(x)) y else x

read_json_safely <- function(path) {
  tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE), error = function(e) NULL)
}

#collects data from result.json
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

#collects data from ranking_debug.json
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

#find the relevant files
candidate_dirs <- unique(dirname(list.files(
  base_dir,
  pattern = "ranking_debug\\.json|ranked_.*\\.pdb|result_.*\\.json",
  full.names = TRUE, recursive = TRUE
)))

if (length(candidate_dirs) == 0) {
  message("no AFM output data found", base_dir)
#quits script if not output data found
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

#create ranked list
ranked <- df |>
  mutate(iptm = suppressWarnings(as.numeric(iptm))) |>
  filter(!is.na(iptm)) |>
  arrange(desc(iptm)) |>
  mutate(rank = dplyr::row_number()) |>
  select(rank, iptm, ptm, mean_plddt, model, directory, source, file)

#create csv, print head
out_csv <- file.path(output_dir, "alphafold_iptm_ranking.csv")
readr::write_csv(ranked, out_csv, na = "")

cat(sprintf("Rangliste gespeichert unter: %s\n", out_csv))
print(utils::head(ranked, 10), row.names = FALSE)
