#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 00-make_probe_key.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Build the SomaScan probe -> gene symbol lookup used by
#'   01-WSSmod_prep.Rmd from the probe list used to construct the WSS modules
#'
#' Author: Robert B Butler III
#'
#' Date Created: 2026-11-09
#'
#' Copyright (c) 2026
#' Email: rrbutler@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'   Internal only (needs the probe list); collaborators receive the output.
#'   Uses the same one-probe-per-gene choice as 09-WSS_construction.R so
#'   WSSmod::calculate_WSS() reproduces the original module scores.
#'
#' ---------------------------
library(optparse)
library(dplyr)
library(readr)
library(stringr)

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "../../dataProt_SS-205063.probe_list.csv",
    help = "Probe list with ix and GENE columns [default: %default]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "soma_probe_to_symbol_key.csv",
    help = "Output key [default: %default]"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

#' SomaScan SeqId from a probe name
#'
#' @param ix Probe names like `CRYBB2.10000.28.3`.
#' @return SeqIds like `seq.10000.28` (the last field is the assay version).
probe_to_seq <- function(ix) {
  str_replace(ix, ".*\\.(\\d+)\\.(\\d+)\\.\\d+$", "seq.\\1.\\2")
}

key <- read_csv(opt$input, show_col_types = FALSE) %>%
  transmute(ix, seq = probe_to_seq(ix), Symbol = GENE) %>%
  filter(!is.na(Symbol), Symbol != "")

# a multi-gene "complex" probe repeats once per constituent gene
stopifnot(!anyDuplicated(key$Symbol), !anyDuplicated(key[c("ix", "Symbol")]))
write_csv(key, opt$output)
message(sprintf("Wrote %d probes to %s", nrow(key), opt$output))
