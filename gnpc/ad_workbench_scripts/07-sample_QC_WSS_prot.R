#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 07-sample_QC_WSS_prot.R
##
## Version: 0.0.1
##
## Purpose of script: Prepare sample data for linear regression
##
## Author: Robert B Butler III, Michael P Brown, Edward N Wilson III
##
## Date Created: 2026-03-06
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
library(data.table)
library(RNOmni)
library(dplyr)     # for data wrangling
library(ggplot2)   # for awesome plotting
library(patchwork)
library(ggrepel)
library(RColorBrewer)
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
    default = "GNPC.imputed.csv",
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
    default = "",
    help = "File path to PCs csv [default: %default]"
  ),
  make_option(
    c("-e", "--eigengenes"),
    type = "character",
    default = "celltype_scores/GNPC.STRING.INT.WSS.csv",
    help = "File path to eigengenes csv [default: skip]"
  ),
  make_option(
    c("-c", "--outcome"),
    type = "character",
    help = "String specifying the column in data containing the outcome variable [required]"     
  ),
  make_option(
    c("--no-svs"),
    action = "store_false",
    dest = "svs",
    default = TRUE,
    help = "skip surrogate variable calculation? [default: no]"
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
  "svs",
  "nameset"
)

# Nameset default case
if(is.null(args_list$nameset)){
  args_list$nameset <- nameset_default(args_list$input)
}

# check outpath status
if (!dir.exists(args_list$output)) {
  dir.create(args_list$output, recursive = TRUE)
}

# Start logging
log_file <- file.path(args_list$output, 
                      paste(args_list$nameset, "INTnorm.log", sep=".")
)
sink(log_file, split = TRUE)

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

# global variables
covs <- strsplit(covs,",")[[1]]
orig_outcome <- outcome
outcome <- make.names(outcome)
qc_cols <- c("person_id", "sample_id", "sequential_visit_number")


# Main -----------------------
# So, the final INTnorm dataset:
a <- fread(input)

# for the GNPC dataset, the prots are separate file to be merged
prot <- transpose(fread(prots), make.names = "GENE", keep.names = "sample_id") 
p_names <- names(prot)[-1]
if (!file.exists("prot_names.txt")) {
  fwrite(list(p_names), file = "prot_names.txt")
}

# Convert -1 to pseudocount (shouldn't apply to imputed samples)
prot[prot == -1] <- 0.01

a <- merge(a, prot, by = "sample_id")


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
  a <- merge(a, x, by.x = "sample_id", by.y = "SampleID")
}

# Initial QC filters
a = a[!is.na(get(outcome))]
cols = c(qc_cols, covs)

# save some mem on GNPC's tiny instances
remove(prot, x)
gc()


# Setup for baseline visit only analysis ---------------------------
# # Correction, keep only first visit where plasma available
# a = a[order(sequential_visit_number)][, head(.SD, 1), by = "person_id"]

# checking for missed duplicates or bad visit labeling
a[duplicated(person_id), ..qc_cols]

# Inverse rank normalization (INT) to z-scale
if (eigengenes != "") { # don't INT normalize WSS/eigengenes
  cols <- p_names
} else {
  cols <- c(outcome, p_names)
}

a[, (cols) := lapply(.SD, RankNorm), .SDcols=cols]

# handling incomplete cases
cols <- c("sample_id", outcome, covs)
drops <- a[!complete.cases(a[, ..cols]), sample_id]
a <- a[!sample_id %in% drops]

# generate surrogate variables ---------------------------------------
if (svs) {
  print("Generating surrogate variables...")
  
  # matrix of proteins
  cols <- c("sample_id", p_names)
  prot <- data.frame(a[, ..cols], row.names = "sample_id")
  
  # matrix of phenotypes to regress
  cols <- c("sample_id", outcome, covs)
  pheno <- data.frame(a[, ..cols], row.names = "sample_id")
  all(rownames(prot) == rownames(pheno))
  
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
  
  svobj <- sva(t(prot), mod, mod0, n.sv = 6)
  
  # extract and add to a
  svs <- data.table(svobj$sv)
  setnames(svs, gsub("V", "sv", names(svs)))
  a <- cbind(a, svs)
  
  # add to covs 
  covs <- c(covs, names(svs))
}

# GLM based feature selection -----------------------------------------
print("Testing feature importance for covariates...")
cols <- c("sample_id", outcome, covs)

# drop missing values for feature importance testing
pheno <- data.frame(a[, ..cols], row.names = "sample_id")
pheno <- pheno[complete.cases(pheno),] 

# formula
fmla <- reformulate(".", response = outcome)
# number of features
n_features <- ncol(pheno) - 1

# include all possible main effects
model = glm(fmla, data = pheno) 

# simple training
set.seed(123)
cv_model <- train(
  form = fmla, 
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
  fmla, 
  data = pheno, 
  method = "pls",
  trControl = trainControl(method = "cv", number = 5),
  preProcess = c("zv", "center", "scale"),
  tuneGrid = data.frame(ncomp = 1:min(10, n_features))
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
  fmla, 
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
    formula = fmla, 
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
  formula = fmla, 
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
  formula = fmla, 
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
# shield improper outcome
setnames(a, 
         old = c(outcome),
         new = c(orig_outcome) 
)

cols = c(qc_cols, covs, orig_outcome, p_names)
fwrite(
  a[, ..cols], 
  file = file.path(output, paste(nameset, "INTnorm.csv.gz", sep="."))
)

sink()




