#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 03-Protein_QC.R
##
## Version: 0.0.10
##
## Purpose of script: Prepare LnPlasmaPTau data for gene set analyses
##
## Author: Robert B Butler III, Michael P Brown
##
## Date Created: 2022-12-21
##
## Copyright (c) 2022
## Email: rrbutler@stanford.edu, mpb@stanford.edu
##
## ---------------------------
##
## Notes:
##   v0.0.10 - convertGeneSymbol() no longer uses biomaRt (broken upstream,
##     see Huber-group-EMBL/biomaRt#161); now reads Ensembl release 115's
##     FTP-published entrez xref table + GTF directly, cached locally on
##     first use. Release 115 (published 2025-09-03) is the closest FTP
##     snapshot to the sep2025.archive.ensembl.org host this previously
##     pinned to, and resolving locally from flat files sidesteps both the
##     archive-browsing breakage and remart's lack of an entrezgene_id
##     filter (checked; not supported as of Huber-group-EMBL/remart@1f7f140).
##
##   Post linear regression, filters based on protein qualities, dropping
##   non-human and non-protein markers. Then matches to Ensembl gene IDs.
##
##   For spots detecting gene complexes it splits the genes to duplicate rows,
##   then for each gene, it selects the row with the highest effect size by
##   absolute value of beta. Interesting genes should be checked to see if they
##   were detected across multiple spots.
##
##   Finally, prints some plots, and calculates some lists of genes for
##   downstream analysis, including a gprofiler enrichment check using whole
##   SomaLogic set as custom background
##
##   REMEMBER, WHEN COMMENTING IN CODE, DO NOT REFERENCE INDIVIDUAL INFORMATION:
##   SAMPLES/PIDNS, DATES OR ANY INDIVIDUAL LEVEL INFO
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
library(data.table)
library(ggplot2)
library(ggpubr)
library(gprofiler2)
library(patchwork)
library(qqman)
library(optparse)
library(tools)
library(Cairo)


# Nameset helper function
nameset_default <- function(input) {
  file_path_sans_ext(basename(input))
}

# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = ".",
    help = "File path to DE input file [default: %default]"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = ".",
    help = "File path to output directory [default: %default]"
  ),
  make_option(
    c("-r", "--ref"),
    type = "character",
    default = "..",
    help = "File path to directory of reference data [default: %default]"
  ),
  make_option(
    c("-p", "--plot-col"),
    dest = "plot_col",
    type = "character",
    default = "beta",
    help = "For plots, which effect size column  [default: %default]"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# Nameset default case
if (is.null(args_list$nameset)) {
  args_list$nameset <- nameset_default(args_list$input)
}

# read in fields
fields <- c(
  "input",
  "output",
  "ref",
  "plot_col",
  "nameset"
)

# check for values
for (field in fields) {
  if (!field %in% names(args_list)) {
    cat(str_c('Argument "', field, '" must be provided\n'), fill = TRUE)
    stop(parse_args(OptionParser(option_list = option_list), args = "--help"),
      call. = FALSE
    )
  }

  assign(field, args_list[[field]])
  print(paste0(field, ": ", get(field)))
}


# Functions -------------------------------------
#' Download (if not already cached) an Ensembl release's gene-to-Entrez
#' xref table and gene annotation GTF from the Ensembl FTP site
#'
#' @param ref_dir directory to cache the downloaded files in
#' @param release Ensembl release number
#'
#' @return list(entrez = path, gtf = path) to the two cached files
ensembl_ftp_cache <- function(ref_dir, release = 115) {
  dir.create(ref_dir, showWarnings = FALSE, recursive = TRUE)
  base <- sprintf("https://ftp.ensembl.org/pub/release-%d", release)
  entrez_f <- file.path(ref_dir, sprintf("Homo_sapiens.GRCh38.%d.entrez.tsv.gz", release))
  gtf_f <- file.path(ref_dir, sprintf("Homo_sapiens.GRCh38.%d.gtf.gz", release))
  if (!file.exists(entrez_f)) {
    download.file(
      sprintf("%s/tsv/homo_sapiens/Homo_sapiens.GRCh38.%d.entrez.tsv.gz", base, release),
      entrez_f, quiet = TRUE
    )
  }
  if (!file.exists(gtf_f)) {
    download.file(
      sprintf("%s/gtf/homo_sapiens/Homo_sapiens.GRCh38.%d.gtf.gz", base, release),
      gtf_f, quiet = TRUE
    )
  }
  list(entrez = entrez_f, gtf = gtf_f)
}

#' for a given column of NCBI IDs, replace them with ensembl gene ids
#' includes a filter for autosomal genes that are not pseudo or TEC
#'
#' @param dt data table to replace names
#' @param colname string, name of the column containing NCBI gene IDs
#' @param symname string, name of the column containing gene Symbols
#' @param keep.symbols boolean, retain the symbols column?
#'
#' @return dt with a GENE column
convertGeneSymbol <- function(dt, colname, symname, keep.symbols = FALSE) {
  print(paste(length(unique(dt[[colname]])), "genes to look up"))

  ff <- ensembl_ftp_cache(file.path(dirname(ref), "reference/ensembl_115"))

  # entrez xref table (file is transcript-level; dedup to gene-level pairs)
  entrez_dt <- fread(cmd = sprintf("zcat %s", shQuote(ff$entrez)), sep = "\t")
  entrez_dt <- entrez_dt[db_name == "EntrezGene"]
  entrez_dt <- unique(entrez_dt[, .(gene_stable_id, xref)])
  setnames(entrez_dt, c("gene_stable_id", "xref"), c("ensembl_gene_id", "entrezgene_id"))
  entrez_dt[, entrezgene_id := as.integer(entrezgene_id)]

  # GTF gene-feature rows only (prefiltered with awk before fread, rather
  # than parsing the whole transcript/exon-level file)
  gtf_dt <- fread(
    cmd = sprintf("zcat %s | awk -F'\\t' '$3==\"gene\"'", shQuote(ff$gtf)),
    header = FALSE, sep = "\t", quote = ""
  )
  setnames(gtf_dt, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))
  gtf_dt[, gene_id := sub('.*gene_id "([^"]+)".*', "\\1", attribute)]
  gtf_dt[, gene_biotype := sub('.*gene_biotype "([^"]+)".*', "\\1", attribute)]
  gtf_dt[, gene_name := ifelse(grepl("gene_name", attribute), sub('.*gene_name "([^"]+)".*', "\\1", attribute), NA_character_)]
  gtf_dt <- gtf_dt[, .(gene_id, external_gene_name = gene_name, gene_biotype, chromosome_name = seqname)]

  # match with ensembl_gene_id across both
  x <- merge(entrez_dt, gtf_dt, by.x = "ensembl_gene_id", by.y = "gene_id")
  x <- x[entrezgene_id %in% unique(dt[[colname]])]
  x <- x[!is.na(external_gene_name) & external_gene_name != ""]
  print(paste(nrow(x), "total matches returned"))
  # ensembl transcript types to keep, excluding pseudogenes and TECs
  keep_categ <- grep("pseudo", unique(x$gene_biotype), value = TRUE, invert = TRUE)
  keep_categ <- grep("TEC", keep_categ, value = TRUE, invert = TRUE)
  x <- x[gene_biotype %in% keep_categ]
  print(paste(nrow(x), "genes correct biotype"))
  # drop duplicates on patches and scaffolds if exist on main chromosomes
  x[, chr_length := nchar(chromosome_name)]
  x <- unique(
    x[order(entrezgene_id, external_gene_name, chr_length)],
    by = c("entrezgene_id", "external_gene_name")
  )
  setnames(x,
    old = c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
    new = c("GENE", "NCBI_ID", "Symbol")
  )
  # replace gene names
  dt <- merge(
    x[, .(GENE, NCBI_ID, Symbol)], 
    dt, 
    all.y = TRUE, 
    by.x = c("NCBI_ID", "Symbol"),
    by.y = c(colname, symname)
  )

  # protein complexes will not match by symbol, need to update them separately
  # convert Symbol to list and match by NCBI and any symbol in list
  cols <- names(dt)[names(dt) != "Symbol"]
  dt[, Symbol := as.list(strsplit(Symbol, "\\|"))]
  dt[UniProt %like% "\\|", GENE := mapply(function(z, y) {
    return(x[NCBI_ID == z & Symbol %in% y, GENE])
  }, NCBI_ID, Symbol), by = .I]
  # collapse Symbol back to a character
  dt[, Symbol := paste(unlist(Symbol), collapse = "|"), by = .I]
  dt$Symbol <- unlist(dt$Symbol) # for the life of me this unlist will not work inside the line above

  # write to file
  print(paste(length(unique(dt[["GENE"]])), "unique genes with ensembl ids"))
  if (keep.symbols == FALSE) dt[, NCBI_ID := NULL]
  return(unique(dt))
}

#' for a string of ens ids: explode, look up symbols, return symbols as string
#'
#' @param string comma separated ens_ids
#' @param dds DESeq Dataset with annotations in mcol
#' @param col string, colname to grab
#'
#' @return string of concatenated, comma separated symbols
get_mcols <- function(str, dds, col) {
  require("stringr")
  y <- dds[[col]][match(unlist(str_split(str, ",")), dds[["GENE"]])]
  return(paste0(y, collapse = ","))
}


#' Filter Extreme Coefficient Estimates from Logistic Regression
#'
#' Identifies numerically unstable or biologically implausible coefficient 
#' estimates that arise from complete or quasi-complete separation in logistic 
#' regression models. Vectorized version for use with data.frames/data.tables.
#'
#' @param beta Numeric vector. The coefficient estimates (log odds ratios) from models.
#' @param se Numeric vector. The standard errors of the coefficient estimates.
#' @param pval Numeric vector. The p-values associated with the coefficient estimates.
#'
#' @return Logical vector. Returns TRUE for estimates that are extreme and should 
#'   be filtered out, FALSE for estimates that should be retained.
#'
#' @details
#' The function identifies three patterns of problematic estimates:
#' \describe{
#'   \item{Pattern 1: Astronomical values}{Beta or SE > 1000, indicating 
#'     numerical instability from separation}
#'   \item{Pattern 2: Extreme beta/SE ratios}{Beta > 100 with beta/SE ratio 
#'     > 10^6, characteristic of complete separation}
#'   \item{Pattern 3: Numerical underflow}{P-value exactly 0 with large beta 
#'     (> 50), indicating p-value underflow from extreme z-statistics}
#' }
#'
#' @note This function handles NA values gracefully by returning FALSE (keep) 
#'   for any observation where beta, se, or pval is NA.
#'
#' @examples
#' # Example with data.table
#' library(data.table)
#' results <- data.table(
#'   protein = c("A", "B", "C"),
#'   beta = c(1e14, 13.5, 2.1),
#'   se = c(1e6, 4.3, 0.52),
#'   pval = c(0, 0.0017, 0.00006)
#' )
#' 
#' results[, exclude := filter_extreme_vec(beta, se, pval)]
#' results[, beta := ifelse(exclude, NA_real_, beta)]
#'
#' # Example with dplyr
#' library(dplyr)
#' results %>%
#'   mutate(
#'     exclude = filter_extreme_vec(beta, se, pval),
#'     beta = ifelse(exclude, NA, beta),
#'     se = ifelse(exclude, NA, se),
#'     pval = ifelse(exclude, NA, pval)
#'   )
#'
#' @export
filter_extreme <- function(beta, se, pval) {
  # Handle NA values
  result <- rep(FALSE, length(beta))
  
  # Pattern 1: Astronomical values
  pattern1 <- !is.na(beta) & !is.na(se) & (abs(beta) > 1000 | se > 1000)
  
  # Pattern 2: Extreme beta/SE ratios
  pattern2 <- !is.na(beta) & !is.na(se) & 
              abs(beta) > 100 & (abs(beta) / se) > 10^6
  
  # Pattern 3: P-value underflow
  pattern3 <- !is.na(pval) & !is.na(beta) & 
              pval == 0 & abs(beta) > 50
  
  # Combine patterns
  result <- pattern1 | pattern2 | pattern3
  
  # Ensure we return FALSE for NAs (keep them, don't filter)
  result[is.na(result)] <- FALSE
  
  return(result)
}


#' Make a volcano with labels & cutoffs, including separated proteins
#'
#' @param x data.table, results 
#' @param titl string, label for graph
#' @param sub string, subtitle for graph
#' @param thresh numeric, cutoff for adjusted p-val
#' @param lfc numeric, log2 fold change cutoff for significance
#' @param lfc_col string, effect_size column
#'
#' @returns ggplot object (patchwork combined plot)
my_volcano <- function(x, titl = NULL, sub = NULL, thresh = 0.05, lfc = 0, lfc_col = "beta") {
  require("ggplot2")
  require("EnhancedVolcano")
  require("patchwork")
  require("ggrepel")
  
  # Split data into separated and non-separated
  x_separated <- x[flag %in% c("separated"), ]
  x_main <- x[!flag %in% c("separated", "failed") & !is.na(padj), ]
  
  # Limit p-values for main plot
  x_main[padj < 1e-150, padj := 1e-150]
  
  # Color by up/downregulated for main plot
  keyvals <- ifelse(
    x_main[[lfc_col]] < -lfc & x_main[["padj"]] < thresh, 
    "#2d004b",
    ifelse(x_main[[lfc_col]] > lfc & x_main[["padj"]] < thresh, "#7f3b08", "grey45")
  )
  names(keyvals)[keyvals == "#7f3b08"] <- "upregulated"
  names(keyvals)[keyvals == "grey45"] <- "not significant"
  names(keyvals)[keyvals == "#2d004b"] <- "downregulated"
  
  # Main volcano plot
  p1 <- EnhancedVolcano(x_main,
    x = lfc_col,
    y = "pval",
    lab = x_main[["Target"]],
    xlim = c(min(x_main[[lfc_col]], na.rm = TRUE) * 1.1, max(x_main[[lfc_col]], na.rm = TRUE) * 1.1),
    ylim = c(0, max(-log10(x_main[["pval"]]), na.rm = TRUE) * 1.25),
    selectLab = x_main[abs(get(lfc_col)) > lfc & padj < thresh, Target],
    labFace = "bold",
    drawConnectors = TRUE,
    widthConnectors = 0.25,
    max.overlaps = 40,
    # maxoverlapsConnectors = 100,  # Changed: allow all connectors
    arrowheads = FALSE, 
    xlab = lfc_col,
    title = NULL,
    subtitle = NULL,
    pCutoff = NA,
    FCcutoff = lfc,
    cutoffLineCol = "light gray",
    pointSize = 1.5,
    labSize = 3.5,
    colCustom = keyvals,
    colAlpha = 1,
    legendPosition = 'none',
    legendLabSize = 11,
    legendIconSize = 5,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    border = "full",
    borderWidth = 0.5,
    borderColour = "black"
  ) + theme(
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )
  
  # If no separated proteins, return main plot only
  if (nrow(x_separated) == 0) {
    return(p1 + 
           ggtitle(titl, subtitle = sub) + 
           theme(plot.title = element_text(hjust = 0.5),
                 plot.subtitle = element_text(hjust = 0.5)))
  }
  
  # Get x-axis limits from main plot
  x_limits <- layer_scales(p1)$x$range$range
  
  # Create separated proteins plot
  x_separated[, y_position := 1]
  x_separated[, y_jitter := y_position + rnorm(.N, 0, 0.15)]
  
  # Color separated points by direction of Cohen's d
  x_separated[, color := ifelse(get(lfc_col) > 0, "#7f3b08", "#2d004b")]
  
  p2 <- ggplot(x_separated, aes(x = .data[[lfc_col]], y = y_jitter)) +
    geom_point(aes(color = color), size = 2, alpha = 0.7) +
    geom_text_repel(
      aes(label = Target),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.3,
      point.padding = 0.3,
      segment.size = 0.2,
      min.segment.length = 0,
      fontface = "bold"
    ) +
    scale_color_identity() +
    scale_x_continuous(limits = x_limits) +  # Match x-axis limits
    scale_y_continuous(
      breaks = 1,
      labels = "Quasi-complete\nseparation",
      limits = c(0.5, 1.5)
    ) +
    labs(
      x = lfc_col,
      y = NULL
    ) +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 10, hjust = 0.5, face = "italic"),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_blank(),  # Hide x-axis text (shared with bottom)
      axis.title.x = element_blank(), # Hide x-axis title (shared with bottom)
      axis.ticks.x = element_blank(), # Hide x-axis ticks
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(t = 5, r = 5, b = 0, l = 5),
      panel.background = element_rect(fill = "white"),  # White background
      plot.background = element_rect(fill = "white")    # White background
    )
  
  # Combine plots
  combined <- p2 / p1 + 
    plot_layout(heights = c(1, 4)) +
    plot_annotation(
      title = titl,
      subtitle = sub,
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11)
      )
    )
  
  return(combined)
}

#' Plot gprofiler enrichment results for multiple gene sets
#'
#' Creates three gostplot-style Manhattan plots showing pathway enrichment for 
#' different gene sets from a gprofiler analysis.
#'
#' @param g List. The gprofiler object returned from gost() containing a 
#'   'result' component with enrichment results.
#' @param n_labels Integer. Number of max.overlaps to allow for labels. Higher is more.
#'   Default is 5.
#'
#' @return A named list of three ggplot objects:
#'   \describe{
#'     \item{all}{Plot for all significant proteins}
#'     \item{updown}{Plot comparing upregulated vs downregulated proteins}
#'     \item{topbottom}{Plot comparing top 20 vs bottom 20 proteins}
#'   }
#'
#' @details
#' The function applies the following filters to enrichment results:
#' \itemize{
#'   \item Term size < 1500
#'   \item Intersection size > 3
#' }
#'
#' @examples
#' \dontrun{
#' plots <- plot_gprofiler_sets(g, n_labels = 10)
#' plots$all
#' plots$updown
#' plots$topbottom
#' }
#'
#' @export
plot_gprofiler_sets <- function(g, n_labels = 5) {
  require("gprofiler2")
  require("ggplot2")
  require("ggrepel")
  require("data.table")
  
  # Filter results
  dt <- as.data.table(g$result)[
    term_size < 1500 & 
    intersection_size > 3
  ]
  
  # Create filtered g objects for each plot
  g_all <- g
  g_all$result <- dt[query == "all"]
  
  g_updown <- g
  g_updown$result <- dt[query %in% c("up", "down")][order(-query, p_value)]
  
  g_topbottom <- g
  g_topbottom$result <- dt[query %in% c("top20", "bottom20")][order(-query, p_value)]
  
  # Plot 1: all
  p1 <- gostplot(g_all, interactive = FALSE) +
    geom_text_repel(
      mapping = aes(label = term_name), 
      family = "Arial", 
      # nudge_y = 4,
      force = 100, 
      force_pull = 0.0001,
      max.overlaps = 5,
      ylim = c(10, NA)
    )
  
  # Plot 2: up/down
  p2 <- gostplot(g_updown, interactive = FALSE) +
    geom_text_repel(
      mapping = aes(label = term_name), 
      family = "Arial", 
      # nudge_y = 4,
      force = 100, 
      force_pull = 0.0001,
      max.overlaps = 5,
      ylim = c(10, NA)
    )
  
  # Plot 3: top20/bottom20
  p3 <- gostplot(g_topbottom, interactive = FALSE) +
    geom_text_repel(
      mapping = aes(label = term_name), 
      family = "Arial", 
      # nudge_y = 4,
      force = 100, 
      force_pull = 0.0001,
      max.overlaps = 5,
      ylim = c(10, NA)
    )
  
  return(list(all = p1, updown = p2, topbottom = p3))
}

# Main ------------------------------------------
# Load assocaition table
a <- fread(input)
a <- a[ix != ""]
# a$variable = gsub("^X\\.", "\\.", a$variable)

# Apply filtering for additional quasi-separation
a[, exclude := filter_extreme(beta, se, pval)]
a[exclude == TRUE, flag := fcase(
  pval < 0.03, "separated",
  pval > 0.97, "failed", 
  default = "separated"
)]
a[, exclude := NULL]  # Remove the exclude column

# handle pvals for quasi-complete separation that would be NAs
a[flag == "separated" & is.na(pval), pval := 0]

# padjust
a[, padj := p.adjust(pval, method = "fdr")]
a[, pbonf := p.adjust(pval, method = "bonferroni")]

# map to reference key
b <- fread(ref)
b <- b[Organism == "Human" & Type == "Protein"]
b[, c("Organism", "Units", "Type", "EntrezGeneSymbol_v2") := NULL]

# Source - https://stackoverflow.com/q/78808739
# Posted by dragon951, modified by community. See post 'Timeline' for change history
# Retrieved 2026-03-13, License - CC BY-SA 4.0

# Explode complex columns with multiple genes/prot entries
mycol <- c("EntrezGeneID", "EntrezGeneSymbol", "UniProt")
b <- b[, lapply(.SD, \(x) unlist(tstrsplit(x, "\\|"))), 
  by = c(names(b)[!names(b) %in% mycol]), 
  .SDcols = mycol
]

# replication cohort probes are SomaScan SeqIds ("seq.10000.28"), so key the
# reference the same way (SYMBOL.10000.28.3 -> seq.10000.28)
if (all(startsWith(a$ix, "seq."))) {
  b[, Key_2 := sub(".*\\.(\\d+)\\.(\\d+)\\.\\d+$", "seq.\\1.\\2", Key_2)]
}
c <- merge(a, b, by.x = "ix", by.y = "Key_2")

# clean up (total set of human protein probes)
c <- c[EntrezGeneSymbol != ""]
c$EntrezGeneID <- as.integer(c$EntrezGeneID)

# drop lower duplicates, use pval then magnitude beta
# c[duplicated(EntrezGeneSymbol) | duplicated(EntrezGeneSymbol, fromLast=T)][order(EntrezGeneSymbol, -beta)]
c <- c[order(pval, -abs(beta))][, head(.SD, 1), by = "EntrezGeneID"]

# lookup NCBI IDs for matching
c <- convertGeneSymbol(c, 
  colname = "EntrezGeneID", 
  symname = "EntrezGeneSymbol", 
  keep.symbols = TRUE
)

fwrite(c, file = file.path(output, paste(nameset, "csv", sep = ".")))

# drop failed for the remainder
c <- c[flag != "failed" | is.na(flag)]

# plot volcanoes -----------------------------------------------
# use cohen's d & val, jitter quasi box at the top
ggsave(
  plot = my_volcano(
    c,
    titl = sprintf("Protein associations with %s", basename(output)),
    sub = NULL,
    thresh = 0.05,
    lfc = 0,
    lfc_col = plot_col
  ),
  filename = paste(nameset, "volcano.pdf", sep = "."),
  path = output, 
  width = 8, 
  height = 8, 
  device = cairo_pdf
)

# plot histograms ----------------------------------------------
# all betas
xmax <- max(abs(c[pval < 0.1, get(plot_col)]))
p1 <- ggplot(c, aes(x = .data[[plot_col]])) +
  geom_density(fill = "#999999") +
  ggtitle("All proteins") +
  theme_classic() +
  xlim(-xmax, xmax) +
  theme(plot.title = element_text(hjust = 0.5))

# calculate genomic inflation
lambda.c <- median(qchisq(c$pval, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
# plot QQ
p3 <- ggplot(c, aes(sample = -log10(pval))) +
  stat_qq(distribution = qexp, dparams = list(rate = log(10))) +
  geom_abline(col = "red", intercept = 0, slope = 1) +
  ggtitle(sprintf("QQ (Lambda: %.3f)", lambda.c)) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# check for significance before writing genesets, pathway enrich, etc
if (nrow(c[padj < 0.05]) < 2) {
  # write plots
  ggsave(
    plot = p1 + p3,
    filename = paste(nameset, paste0(plot_col, "_distributions.pdf"), sep = "."),
    path = output,
    width = 10,
    height = 6,
    device = cairo_pdf
  )

  # exit early
  print("fewer than 2 genes are padj < 0.05, exiting...")
  stop(call. = FALSE)
}

# Bonferroni significant betas
c_cut <- c[padj < 0.05]
# make density distribution for quantiles
dens <- density(c_cut[[plot_col]])
qs <- quantile(c_cut[[plot_col]], c(0.2, 0.8))
dt <- data.frame(x = dens$x, y = dens$y)
dt$quant <- factor(
  findInterval(dt$x, c(min(qs[[1]], -.Machine$double.xmin), 0, qs[[2]]))
)
# plot
p2 <- ggplot(dt, aes(x, y)) +
  geom_line() +
  labs(x = plot_col, y = "density", title = "Proteins with Padj < 0.05") +
  geom_ribbon(aes(ymin = 0, ymax = y, fill = quant)) +
  scale_fill_brewer(
    type = "div", 
    palette = "RdYlBu", 
    direction = -1, 
    name = "",
    labels = c("Bottom 10%", "Less than 0", "Greater than 0", "Top 10%")
  ) +
  xlim(-xmax, xmax) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    legend.position.inside = c(0.9, 0.5)
  )

# write plots
ggsave(
  plot = p2 +
    inset_element(p1, left = 0.05, top = 1, bottom = 0.6, right = 0.3) +
    inset_element(p3, left = 0.7, top = 1, bottom = 0.6, right = 0.95),
  filename = paste(nameset, paste0(plot_col, "_distributions.pdf"), sep = "."),
  path = output, 
  width = 10, 
  height = 6, 
  device = cairo_pdf
)

# write gene lists
vars <- list(p = plot_col)
sets <- list(
  all = c[!is.na(GENE) & padj < 0.05, .(GENE, p), env = vars
          ][order(-abs(p)), env = vars],
  up = c[!is.na(GENE) & padj < 0.05 & p > 0, .(GENE, p), env = vars
         ][order(-p), env = vars],
  down = c[!is.na(GENE) & padj < 0.05 & p < 0, .(GENE, p), env = vars
           ][order(p), env = vars],
  top20 = c[!is.na(GENE) & padj < 0.05 & p > qs[[2]], .(GENE, p), env = vars
            ][order(-p), env = vars],
  bottom20 = c[!is.na(GENE) & padj < 0.05 & p < qs[[1]], .(GENE, p), env = vars
               ][order(p), env = vars]
)
# lapply(seq(sets), function(i) {
  # fwrite(sets[[i]],
    # file = file.path(output, paste(nameset, names(sets)[i], "txt", sep = ".")),
    # sep = "\t", 
    # col.names = FALSE
  # )
# })

# gprofiler enrichment
g <- gost(
  query = sets, 
  organism = "hsapiens", 
  ordered_query = FALSE, 
  evcodes = TRUE,
  custom_bg = unique(c[!is.na(GENE), GENE]),
  sources = c("GO", "KEGG", "REAC", "WP")
)
# write to file
if (is.list(g$result)) {
  # get symbols for terms
  g$result$symbols <- lapply(g$result$intersection, get_mcols, c, "Symbol")
  # write to file
  fwrite(g$result, file = file.path(output, paste0(nameset, ".gprofiler.csv")))
}

# plot with labels
plots <- plot_gprofiler_sets(g, n_labels = 10)

CairoPDF(
  file = file.path(output, paste(nameset, "gprofiler_plots.pdf", sep=".")), 
  onefile = TRUE, 
  paper = "USr", 
  width = 11, 
  height = 8.5,
)
plots$all
plots$updown
plots$topbottom
dev.off()


