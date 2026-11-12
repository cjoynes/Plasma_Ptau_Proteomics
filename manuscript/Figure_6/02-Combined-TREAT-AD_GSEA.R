#!/usr/bin/env Rscript
#' ---------------------------
#'
#' Script name: 03-Combined-TREAT-AD_GSEA.R
#'
#' Version: 0.0.1
#'
#' Purpose of script: Comparison of GSEA correlations of TREAT-AD in proteomic associations 
#'
#' Author: Robert B Butler III, Michael Peter Brown
#'
#' Date Created: 2024-11-18
#'
#' Copyright (c) 2024
#' Email: rrbutler@stanford.edu, mpb@stanford.edu
#'
#' ---------------------------
#'
#' Notes:
#'
#'
#' ---------------------------
#' load up the packages we will need:  (uncomment as required)
library(optparse)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(data.table)
library(stringr)
library(dplyr)
library(patchwork)
library(gt)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(Cairo)
library(extrafont)
# library(GeneOverlap)
# library(tools)

# Nameset helper function
nameset_default <- function(input){
  file_path_sans_ext(basename(input))
}

# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "List of input datasets comma separated"
  ),
  make_option(
    c("--inpath"),
    type = "character",
    default = "../../WSS_prot",
    help = "Path to input datasets folder [default: %default]"
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
    default = "../../module_scores/STRING.INT.WSS.terms.csv",
    help = "File path to module term names [default: %default]"
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    default = "dataProt_SS-205063.PC.related",
    help = "filename prefix of input project [default: %default]"
  ),
  make_option(
    c("-s", "--suffix"),
    type = "character",
    default = "all_bios",
    help = "output suffix to add to nameset [default: %default]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))

# Nameset default case
if(is.null(args_list$nameset)){
  args_list$nameset <- nameset_default(args_list$input)
}

# read in fields
fields <- c(
  "input",
  "inpath",
  "output",
  "ref",
  "nameset",
  "suffix"
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


# input variables
input <- strsplit(input, ",")[[1]]
input <- c(t(outer(input, c("", ".pred"), paste0)))

# Functions ---------------------------

#' for a string of ens ids: explode, look up symbols, return symbols as string
#'
#' @param str string of comma separated ens_ids
#' @param result data.table with ens_ids in GENE col
#' @param col colname to grab (i.e. Symbol)
#'
#' @return concatenated string of lookup terms 
get_mcols <- function(str, result, col) {
  require("stringr")
  y <- result[[col]][match(unlist(str_split(str, ",")), result$GENE)]
  return(paste0(y, collapse = ","))
}



# Main ------------------------------------------
# DE file locations
fh <- file.path(inpath, input,
  paste(nameset, input, "GSEA_Biodomain.rds", sep = ".")
)
result <- lapply(fh, readRDS)
setattr(result, 'names', input)

# Extract relevant columns from each gseaResult object and convert to data.table
result_dt <- lapply(names(result), function(module_name) {
  
  # Extract the result data.frame from the gseaResult object
  df <- result[[module_name]]@result
  
  # Convert to data.table and select relevant columns
  dt <- as.data.table(df)[, .(
    Module = module_name,
    BioDomain = label,  # or use 'Biodomain' column if that's what you want
    ONTOLOGY = ONTOLOGY,
    ID = ID,
    Description = Description,
    NES = NES,
    padj = p.adjust
  )]
  
  return(dt)
})
result_dt <- rbindlist(result_dt)

# # if predicted scores remove
# result_dt[, Module := gsub(".pred", "", Module)]
# input <- gsub(".pred", "", input)

# Figure the mean NES * proportion significant -------------------------
scores <- result_dt[, .(
  n_total = .N,
  n_significant = sum(padj < 0.05, na.rm = TRUE),
  proportion_sig = sum(padj < 0.05, na.rm = TRUE) / .N,
  mean_NES_sig = ifelse(sum(padj < 0.05, na.rm = TRUE) > 0,
                        mean(NES[padj < 0.05], na.rm = TRUE),
                        NA_real_),
  BioDomain_Score = ifelse(sum(padj < 0.05, na.rm = TRUE) > 0,
                           mean(NES[padj < 0.05], na.rm = TRUE) * 
                             sum(padj < 0.05, na.rm = TRUE) / .N,
                           0)
), by = .(Module, BioDomain)][order(Module, -BioDomain_Score)]

# # Reshape to wide format (Modules as rows, BioDomains as columns)
# dcast(scores, Module ~ BioDomain, value.var = "BioDomain_Score", fill = 0)

# Plotting for characteristics ------------------------------------
p1 <- ggplot(scores, aes(x = BioDomain_Score, y = Module, fill = Module)) +
  geom_density_ridges(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "BioDomain Score Distribution by Module",
    x = "BioDomain Score",
    y = "Module"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

p2 <- scores[, .(
  Mean = mean(BioDomain_Score, na.rm = TRUE),
  Median = median(BioDomain_Score, na.rm = TRUE),
  SD = sd(BioDomain_Score, na.rm = TRUE),
  Min = min(BioDomain_Score, na.rm = TRUE),
  Max = max(BioDomain_Score, na.rm = TRUE),
  N_BioDomains = .N,
  N_Positive = sum(BioDomain_Score > 0),
  N_Negative = sum(BioDomain_Score < 0)
), by = Module] |>
  gt() |>
  fmt_number(columns = c(Mean, Median, SD, Min, Max), decimals = 3) |>
  tab_options(
    table.font.size = 10,
    data_row.padding = px(2),
    column_labels.font.weight = "bold"
  )


ggsave(
  filename = file.path(output, 
    paste(nameset, suffix, "Distribution_BioDomain_Scores.pdf", sep = ".")
  ),
  plot = p1 + p2,
  width = 14
)

# Heatmap -------------------------------------------------------------
# load term names 
names <- fread(ref)
names[, group := tstrsplit(name, "\\.", keep = c(1))]
names[, group := { 
  group %>%
    str_remove_all(., '_|Plasma') %>%
    str_replace_all(., 'AB140P', 'A\u03b240') %>%
    str_replace_all(., 'AB142P', 'A\u03b242') %>%
    str_replace_all(., 'PTau', 'pTau') %>%
    str_replace_all(., 'ABRatio', 'A\u03b242/A\u03b240') %>%
    str_replace_all(., 'AB42Ratio', '/A\u03b242')
}]
g_names <- c("Merged", sort(unique(names[group != "Merged", group])))
names[, group := factor(group, levels = g_names)]

# add in pred names
pred_names <- copy(names)[, name := paste0(name, ".pred")]
pred_names[, term := "Predicted"]
names <- rbindlist(list(names, pred_names))
names <- names[match(input, name)]

# add to scores
scores <- merge(scores, names, by.x = 'Module', by.y = 'name', all.x = TRUE)
b_names <- c(sort(levels(scores$BioDomain))[-13], "none")
scores$BioDomain <- factor(scores$BioDomain, levels = b_names)
scores$Module <- factor(scores$Module, levels = input)
setorder(scores, group, Module, BioDomain)
fwrite(scores, 
  file = file.path(output, 
    paste(nameset, suffix, "BioDomain_Scores.csv", sep = ".")
  )
)

# Convert to matrix
score_matrix <- as.matrix(
  x = dcast(scores, BioDomain ~ Module, value.var = "BioDomain_Score", fill = 0),
  rownames = "BioDomain"
)
prop_matrix <- as.matrix(
  x = dcast(scores, BioDomain ~ Module, value.var = "proportion_sig", fill = 0),
  rownames = "BioDomain"
)

# Create color function
col_fun <- colorRamp2(
  c(min(score_matrix, na.rm = TRUE), 0, max(score_matrix, na.rm = TRUE)), 
  c("#164B6E", "white", "#85070C")
)

# # Create column annotation
# col_ha <- HeatmapAnnotation(
  # Biomarker = names$group[match(colnames(score_matrix), names$name)],
  # col = list(
    # Biomarker = setNames(
      # object = brewer.pal(length(levels(names$group)), "Dark2"), 
      # nm = levels(names$group))
  # ),
  # show_legend = TRUE,
  # annotation_name_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"),
  # annotation_name_side = "left"
# )

# Calculate number of unique GO terms per BioDomain
gos <- unique(result_dt[, .(ID, BioDomain)])[, .N, by = BioDomain]

# Match to the row order of your heatmap
gos_ord <- gos$N[match(rownames(score_matrix), gos$BioDomain)]

# Create right-side row annotation with GO term counts
row_ha <- rowAnnotation(
  `# GO Terms` = anno_text(
    gos_ord,
    gp = gpar(fontsize = 9, fontface = "bold", fontfamily = "Arial")
  ),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"),
  annotation_name_rot = 0
)

# Create heatmap
p <- Heatmap(
  score_matrix,
  name = "BioDomain\nScore",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_order = input,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10, fontface = 'bold', fontfamily = 'Arial'),
  row_names_side = 'left',
  column_names_gp = gpar(fontsize = 9, fontface = 'bold', fontfamily = 'Arial'),
  column_labels = structure(names$term, names = names$name),
  column_names_side = 'top',
  column_names_rot = 45,
  column_split = rep(1:(length(input)/2), each = 2),
  column_gap = unit(3, "mm"),          # Adjust spacing
  border_gp = gpar(col = "black", lwd = 2),  # Dark border line
  column_title = "TREAT-AD BioDomain Enrichment Scores by WSS Module",
  # top_annotation = col_ha,
  right_annotation = row_ha,
  heatmap_legend_param = list(
    title = "Score\nmeanNES * prop(FDR < 0.05)",
    direction = "horizontal",
    title_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"),
    labels_gp = gpar(fontsize = 9, fontfamily = "Arial")
  )
)

CairoPDF(
  file = file.path(output, 
    paste(nameset, suffix, "BioDomain_Scores_Heatmap.pdf", sep = ".")
  ), 
  width = 13, 
  height = 8
)
draw(p, heatmap_legend_side = "bottom")
dev.off()


# Possibly focus on sub-domains ---------------------------------------


# Driver domains/proteins ---------------------------------------------
# p1 <- gt(
  # unique(
    # result_dt[Module == "Merged.M8" 
      # & !BioDomain %in% c("none", "Immune Response [IR]") 
      # & padj < 0.05
    # ][order(-NES)][1:13, .(ID, Description, BioDomain, NES)]
  # )
# ) |> 
  # tab_header(title = md("**Merged.M8 - Synaptic vesicle|SNARE complex**")) |>
  # cols_align(
    # align = "right",
    # columns = c(BioDomain, NES)
  # )

# ggsave(
  # plot = wrap_table(p1), 
  # file = file.path(output, 
    # paste(nameset, suffix, "BioDomain_Top_Terms.png", sep = ".")
  # ),
  # width = 10,
  # height = 4
# )
