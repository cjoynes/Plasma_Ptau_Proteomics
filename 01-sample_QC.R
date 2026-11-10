#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 01-sample_QC.R
##
## Version: 0.0.8
##
## Purpose of script: Prepare LnPlasmaPTau sample data for linear regression
##
## Author: Robert B Butler III, Michael P Brown, Edward N Wilson III
##
## Date Created: 2023-07-04
##
## Copyright (c) 2023
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   As noted in 01a-Examining_SomaLogic_QC.R, it is really important to use the
##   data that is scale to the external reference. When comparing a protein 
##   across samples, it is also important to use the inverse rank normalized
##   data, as that is then scaled approriately to actually compare the lmm betas
##
##   REMEMBER, WHEN COMMENTING IN CODE, DO NOT REFERENCE INDIVIDUAL INFORMATION:
##   SAMPLES/PIDNS, DATES OR ANY INDIVIDUAL LEVEL INFO
##
## ---------------------------

## load up the packages we will need:  (uncomment as required)
options(future.cmdargs = character(0)) # some dependency loads future
library(data.table)
library(RNOmni)
library(dplyr)     # for data wrangling
library(ggplot2)   # for awesome plotting
library(patchwork)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(matrixStats)
library(sva)
# Modeling packages
# library(earth)     # for fitting MARS models
library(pls)     # for fitting PLS models
library(caret)    # for tuning models
library(ranger)   # a c++ implementation of random forest 
library(vip)       # for variable importance
library(optparse)

# Nameset helper function
nameset_default <- function(input){
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
    c("-p", "--prots"),
    type = "character",
    default = "../prot_names.txt",
    help = "File path to prots file [default: %default]"
  ),
  make_option(
    c("-v", "--covs"),
    type = "character",
    default = "",
    help = "Terminal input of covariates in comma-separated format [default: %default]"
  ),
  make_option(
    c("-m", "--pcs"),
    type = "character",
    default = "genomic/mypcair.PCs.csv",
    help = "File path to PCs csv [default: %default]"
  ),
  make_option(
    c("-e", "--eigengenes"),
    type = "character",
    default = "",
    help = "File path to eigengenes csv [default: skip]"
  ),
  make_option(
    c("-c", "--outcome"),
    type = "character",
    help = "String specifying the column in data containing the outcome variable [required]"     
  ),
  make_option(
    c("--skip-qc"),
    action = "store_false",
    default = TRUE,
    dest = "qc",
    help = "Skip sample qc? [default: Run QC]" 
  ),
  make_option(
    c("-n", "--nameset"),
    type = "character",
    help = "Name of project [default: Nameset of input]"
  )
)

args_list <- parse_args(OptionParser(option_list = option_list))


# read in fields
fields <- c(
  "input",
  "output",
  "prots",
  "covs",
  "pcs",
  "eigengenes",
  "outcome",
  "qc",
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

# Nameset default case
if(is.null(args_list$nameset)){
  args_list$nameset <- nameset_default(args_list$input)
}

# check outpath status
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

# variables
covs <- strsplit(covs,",")[[1]]
prots <- readLines(prots)
# first_prot = prots[1]                                             
# last_prot = prots[length(prots)]

qc_cols <- c(
  "PIDN", "SampleID", "ADRC_ID", "Visit", "Visit_Consensus", "ConnectivityZscore",
  "Diagnosis_group", "Date.of.draw", "DS1V"
)

# Main -----------------------
# So, the final INTnorm dataset:
a <- fread(input)

# Map in PCs assumes your primary data table exists 
if (pcs != "") {
  print(sprintf("Reading PCs for inclusion from %s", pcs))
  x <- fread(pcs, colClasses = c(PIDN = "character"))
  a <- merge(a, x, by = "PIDN")
  covs <- c(covs, names(a)[names(a) %like% "^PC\\d{1,2}$"]) # flexible for PCs
}

# adding in eigengenes if specified
if (eigengenes != "") {
  print(sprintf("Reading eigengenes for inclusion from %s", eigengenes))
  x <- fread(eigengenes)
  a <- merge(a, x, by = "SampleID")
}


# Initial QC filters
a = a[!is.na(get(outcome))]
cols = c(qc_cols, covs)


# Sample outlier analysis --------------------------------------------------
if (qc) {
  # check ConnectivityZscore for extreme outliers:
  a[abs(ConnectivityZscore) > 3, ..cols]
  
  # make dummy SampID
  a[, SampID := paste0("Samp", .I)]
  
  # plots connectivity histogram
  p1 = ggplot(a, aes(x=ConnectivityZscore)) + geom_histogram()
  ggsave(
    paste(nameset, "ConnectivityZscore.pdf", sep = "."), 
    p1, 
    useDingbats = F, 
    path = output
  )

  # One sample had an extreme Z-score, less than -23. Participant's baseline visit
  # was in a normal range, and there was no other metadata explaining the outlier.
  # The rest included some samples with z in the +/- 3-4 range, but not notable

  # fetch the protein columns for plotting
  cols = c("SampID", prots) 
  b = a[, ..cols]

  # Poisson Distance Pairwise
  poisd = PoissonDistance(as.matrix(b, rownames="SampID"))
  samplePoisDistMatrix = as.matrix(poisd$dd)
  rownames(samplePoisDistMatrix) = b[, SampID]
  colnames(samplePoisDistMatrix) = b[, SampID]

  # set colors and draw
  colors = colorRampPalette(brewer.pal(11, "RdYlBu"))(255)
  pheatmap(samplePoisDistMatrix, col=colors, clustering_distance_rows=poisd$dd,
           clustering_distance_cols=poisd$dd, width=110, height=85, border_color=NA,
           filename=file.path(output, paste(nameset, "Poisson_distance.pdf", sep=".")))

  # Some groupings that could be interesting, one major outlier is the same one
  # from z-score connectivity, final check PCA

  # PCA
  # Function ---------------------------
  # PCA top X genes by variance
  # c_by = column to color by
  # c_lab = Label for color legend
  # c_type = "disc" discrete or "cont" continuous coloring
  # num = number of genes, default is all
  my_pca_plot = function(vst, c_by, c_lab, c_type="disc", num=1000000) {
    rv = rowVars(vst)
    top_genes = order(rv, decreasing=T)[seq_len(min(num, length(rv)))]
    pca = prcomp(t(vst[top_genes, ]))
    percent_pc = round((pca$sdev^2 / sum(pca$sdev^2))*100, 2)
    if (c_type == "disc") {
      col_by = as.factor(a[[c_by]])
    } else if (c_type == "cont") {
      col_by = a[[c_by]]
    } else {
      print("Incorrect value for c_type!")
      quit("no", 1)
    }

    # PC1 v PC2
    ggPCA12 = ggplot(data.frame(pca$x), aes(x=PC1, y=PC2, col=col_by)) +
      labs(x=sprintf("PC1 - %.2f%%", percent_pc[1]),
           y=sprintf("PC2 - %.2f%%", percent_pc[2]), col=c_lab) +
      geom_point(size=3, alpha=0.5) +
      geom_text_repel(aes(label=rownames(data.frame(pca$x)))) +
      ggtitle(sprintf("top %s genes", length(top_genes))) + theme_classic() + {
        if(c_type == "disc") scale_color_brewer(type="qual", palette="Paired") 
        else scale_color_distiller(type="div", palette="RdYlBu")
      }

    # PC1 v PC3
    ggPCA13 = ggplot(data.frame(pca$x), aes(x=PC1, y=PC3, col=col_by)) +
      labs(x=sprintf("PC1 - %.2f%%", percent_pc[1]),
           y=sprintf("PC3 - %.2f%%", percent_pc[3]), col=c_lab) +
      geom_point(size=3, alpha=0.5) +
      geom_text_repel(aes(label=rownames(data.frame(pca$x)))) +
      ggtitle(sprintf("top %s genes", length(top_genes))) + theme_classic() + {
        if(c_type == "disc") scale_color_brewer(type="qual", palette="Paired") 
        else scale_color_distiller(type="div", palette="RdYlBu")
      }
    plot(ggPCA12)
    plot(ggPCA13)
  }

  # Plot PCAs
  pdf(file=file.path(output, paste(nameset, "pca_plots.pdf", sep=".")), onefile=T, paper="USr", 
      width=11, height=8.5)
  my_pca_plot(t(as.matrix(b, rownames="SampID")), c_by="Diagnosis_group", 
              c_lab="Diagnosis_group", c_type="disc")
  my_pca_plot(t(as.matrix(b, rownames="SampID")), c_by="Diagnosis_group", 
              c_lab="Top 500 Diagnosis_group", c_type="disc", num=500)
  dev.off()

  # Same outlier sample in PCs as well, same across three different metrics, remove
  # drop outlier samples
  drop = a[abs(ConnectivityZscore) > 10, SampleID]
  a = a[!SampleID %in% drop]
  
  # drop dummy names
  a[, SampID := NULL]
  print("Finished QC...")
}

# Setup for baseline visit only analysis ---------------------------
# # # drop all samples with DS1V != 0
# # a = a[DS1V == 0]
# Correction, keep only first visit where plasma available
a = a[order(DS1V)][, head(.SD, 1), by = "PIDN"]

# checking for missed duplicates or bad visit labeling
a[duplicated(PIDN), ..qc_cols]

# Inverse rank normalization (INT) to z-scale
if (eigengenes != "") { # don't INT normalize WSS/eigengenes
  cols <- prots
} else {
  cols <- c(outcome, prots)
}

a[, (cols) := lapply(.SD, RankNorm), .SDcols=cols]

# handling incomplete cases
cols <- c("SampleID", outcome, covs)
drops <- a[!complete.cases(a[, ..cols]), SampleID]
a <- a[!SampleID %in% drops]

# generate surrogate variables ---------------------------------------
# matrix of proteins
cols <- c("SampleID", prots)
prot <- t(as.matrix(a[, ..cols], rownames = "SampleID"))

# matrix of phenotypes to regress
cols <- c("SampleID", outcome, covs)
pheno <- data.frame(a[, ..cols], row.names = "SampleID")
all(colnames(prot) == rownames(pheno))

# defining models with static set of 6 SVs
mod <- model.matrix( # this will include all covs and the outcome
  as.formula(paste("~", paste(names(pheno), collapse=" + "))), 
  data = pheno
) 
mod0 <- model.matrix( # this will drop the outcome
  as.formula(
    paste("~", paste(names(pheno)[names(pheno) != outcome], collapse = " + "))
  ), 
  data = pheno
)
svobj <- sva(prot, mod, mod0, n.sv = 6)

# extract and add to a
svs <- data.table(svobj$sv)
setnames(svs, gsub("V", "sv", names(svs)))
a <- cbind(a, svs)

# add to covs 
covs <- c(covs, names(svs))
cols <- c("SampleID", outcome, covs)

# GLM based feature selection -----------------------------------------
# examined potential covariates (use pheno data.frame above)
# drop missing values for feature importance testing
pheno <- data.frame(a[, ..cols], row.names = "SampleID")
pheno <- pheno[complete.cases(pheno),] 

# include all possible main effects
model = glm(as.formula(paste(outcome, "~", ".")), data = pheno) 

# simple training
set.seed(123)
cv_model <- train(
  as.formula(paste(outcome, "~", ".")), 
  data = pheno, 
  method = "glm",
  trControl = trainControl(method = "cv", number = 5)
)

# summary(cv_model$finalModel)

# Partial Least Squares feature scoring -------------------------------
# perform 10-fold cross validation on a PLS model tuning the 
# number of principal components to use as predictors from 1-30

# # custom model npcs
# npcs = 39
# model_pls = mvr(outcome ~ ., data=b, ncomp=npcs) 

set.seed(123)
cv_model_pls <- train(
  as.formula(paste(outcome, "~", ".")), 
  data = pheno, 
  method = "pls",
  trControl = trainControl(method = "cv", number = 5),
  preProcess = c("zv", "center", "scale"),
  tuneLength = 15
)

# model with lowest RMSE
cv_model_pls$bestTune
## ncomp
##   1

# results for model with lowest RMSE
cv_model_pls$results %>%
  dplyr::filter(ncomp == pull(cv_model_pls$bestTune))
##  ncomp      RMSE Rsquared       MAE     RMSESD RsquaredSD     MAESD
##      1 0.9399164 0.133929 0.7375119 0.09900189 0.07164967 0.0701038
# plot cross-validated RMSE
# ggplot(cv_model_pls)
#' Variable importance dot plot
#'
#' Replacement for `vip::vip(geom = "point")`: vip >= 0.5 draws with base
#' graphics instead of returning a ggplot, so build the plot from `vip::vi()`,
#' whose values are identical across vip versions.
#'
#' @param object A fitted model supported by `vip::vi()`.
#' @param num_features Number of top features to show.
#' @return A ggplot with the most important feature at the top.
vip_plot <- function(object, num_features = 50) {
  vip::vi(object) %>%
    slice_max(Importance, n = num_features, with_ties = FALSE) %>%
    mutate(Variable = reorder(Variable, Importance)) %>%
    ggplot(aes(x = Importance, y = Variable)) +
    geom_point() +
    ylab(NULL)
}

p1 <- vip_plot(cv_model) + 
  ggtitle("GLM - 5-fold Cross Validation") + 
  geom_hline(yintercept = min(25, floor(length(cv_model$coefnames) / 2)))
p2 <- vip_plot(cv_model_pls) +
  ggtitle(paste("Partial Least Squares - Best Fit ", 
                cv_model_pls$bestTune$ncomp, "PCs")) + 
  geom_hline(yintercept = min(25, floor(length(cv_model_pls$coefnames) / 2)))


# Random Forest feature scoring ---------------------------------------
# number of features
n_features <- ncol(pheno) - 1

# train a default random forest model
rf1 <- ranger(
  as.formula(paste(outcome, "~", ".")), 
  data = pheno,
  mtry = floor(n_features / 3),
  respect.unordered.factors = "order"
)

default_rmse <- sqrt(rf1$prediction.error)

# create hyperparameter grid
hyper_grid <- expand.grid(
  mtry = floor(n_features * c(.05, .15, .25, .333, .4)),
  min.node.size = c(1, 3, 5, 10), 
  replace = c(TRUE, FALSE),                               
  sample.fraction = c(.5, .63, .8),                       
  rmse = NA                                               
)

# execute full cartesian grid search
for(i in seq_len(nrow(hyper_grid))) {
  # fit model for ith hyperparameter combination
  fit.rf <- ranger(
    formula = as.formula(paste(outcome, "~", ".")), 
    data = pheno,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$min.node.size[i],
    replace         = hyper_grid$replace[i],
    sample.fraction = hyper_grid$sample.fraction[i],
    verbose         = FALSE,
    respect.unordered.factors = 'order',
  )
  # export OOB error 
  hyper_grid$rmse[i] <- sqrt(fit.rf$prediction.error)
}

# assess top 10 models
top = hyper_grid %>%
  arrange(rmse) %>%
  mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
  head(10)
# mtry min.node.size replace sample.fraction      rmse perc_gain
# 22             5   FALSE            0.63 0.9101213 1.2430475

# re-run model with impurity-based variable importance
rf_impurity <- ranger(
  formula = as.formula(paste(outcome, "~", ".")), 
  data = pheno,
  num.trees = 2000,
  mtry = top$mtry[1],
  min.node.size = top$min.node.size[1],
  sample.fraction = top$sample.fraction[1],
  replace = top$replace[1],
  importance = "impurity",
  respect.unordered.factors = "order",
  verbose = FALSE,
)

# re-run model with permutation-based variable importance
rf_permutation <- ranger(
  formula = as.formula(paste(outcome, "~", ".")), 
  data = pheno,
  num.trees = 2000,
  mtry = top$mtry[1],
  min.node.size = top$min.node.size[1],
  sample.fraction = top$sample.fraction[1],
  replace = top$replace[1],
  importance = "permutation",
  respect.unordered.factors = "order",
  verbose = FALSE,
)

p3 <- vip_plot(rf_impurity) + 
  ggtitle("Random Forest - Impurity") + 
  geom_hline(yintercept = min(25, floor(n_features/2)))
p4 <- vip_plot(rf_permutation) +
  ggtitle("Random Forest - Permutation") + 
  geom_hline(yintercept = min(25, floor(n_features/2)))

pdf(
  file = file.path(output, paste(nameset, "feature_importance.pdf", sep=".")), 
  width = 11, 
  height = 11
)
p1 + p2 + p3 + p4
dev.off()

# Post-covariate determination ----------------------------------------------
cols = c(qc_cols, covs, outcome, prots)
fwrite(
  a[, ..cols], 
  file = file.path(output, paste(nameset, "INTnorm.csv.gz", sep="."))
)



# Testing -------------------------------------------------------------------
# a = data.table(SampleID=seq(6), 
               # PIDN=c(rep("A",2),rep("B",3),rep("C",1)), 
               # Visit=c(seq(2), seq(3), 1), 
               # Date.of.draw=c("3/12/2018", "4/23/2019", "12/28/2016", "2/14/2018",
                      # "7/3/2019", "4/9/2012"))
# a[20:30, .(SampleID, PIDN, Visit, Date.of.draw, Visit_Consensus, DS1V)]
# ---------------------------------------------------------------------------
# p <- a %>% ggplot(aes(outcome, CSRP3.9171.11.3))
# alpha = 0.5
# p + 
  # geom_point(aes(col = PlateId), alpha = alpha) + 
  # geom_smooth(aes(col = PlateId), method = 'glm') +
  # theme(legend.background = element_rect(fill = "transparent"),
        # legend.justification = c(0, 1),
        # legend.position = c(0, 1))

# # # # MARS based feature scoring ---------------------------------------
# # # # examined potential covariates
# # # cols = c(names(a)[2:(grep(first_prot, names(a)) - 1)])
# # # cols = cols[! cols %in% c("Diagnosis_group", "PIDN", "Visit", "DS1V")]
# # # b = a[, ..cols]
# # # b = b[complete.cases(b)] # drop missing values for feature importance testing

# # # # earth quick feature scoring to examine covariates
# # # fit <- earth(outcome ~ ., data=b, degree=2)
# # # evimp(fit, trim=FALSE) # trim=FALSE to show unused variables
# # # p0 = ~plot(evimp(fit), main="Earth (MARS) Model Importance")

# # # # create a tuning grid
# # # hyper_grid <- expand.grid(
  # # # degree = 1:2, 
  # # # nprune = seq(5, 200, length.out = 10) %>% floor()
# # # )

# # # # Cross-validated model
# # # set.seed(123)  # for reproducibility
# # # cv_mars <- train(
  # # # x = subset(b, select = -outcome),
  # # # y = b[, outcome],
  # # # method = "earth",
  # # # metric = "RMSE",
  # # # trControl = trainControl(method = "cv", number = 10),
  # # # tuneGrid = hyper_grid
# # # )

# # # # View results
# # # cv_mars$bestTune
# # # ##    nprune degree
# # # ## 16     56      2

# # # cv_mars$results %>%
  # # # filter(nprune == cv_mars$bestTune$nprune, degree == cv_mars$bestTune$degree)
# # # ##   degree nprune    RMSE  Rsquared      MAE   RMSESD RsquaredSD    MAESD
# # # ## 1      2     56 26817.1 0.8838914 16439.15 11683.73 0.09785945 1678.672

# # # ggplot(cv_mars)

