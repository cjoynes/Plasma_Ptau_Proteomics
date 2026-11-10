#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 12-sample_QC_WSS_clin.R
##
## Version: 0.0.1
##
## Purpose of script: preparing WSS scores as predictors for clinical analysis
##
## Author: Robert B Butler III, Michael P Brown, Edward N Wilson III
##
## Date Created: 2026-02-19
##
## Copyright (c) 2026
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   As noted in 01a-Examining_SomaLogic_QC.R, it is really important to use the
##   data that is scale to the external reference. 
##   
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
library(matrixStats)
library(sva)
# Modeling packages
library(pls)     # for fitting PLS models
library(caret)    # for tuning models
library(ranger)   # a c++ implementation of random forest 
library(vip)       # for variable importance
library(optparse)
# library(lme4)


# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "dataProt_SS-205063.merged.QC.csv",
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
    default = "WSS_names.txt",
    help = "File path to prots file [default: %default]"
  ),
  make_option(
    c("-v", "--covs"),
    type = "character",
    default = "",
    help = "Terminal input of covariates in comma-separated format [default: %default]"
  ),
  make_option(
    c("--preds"),
    type = "character",
    default = "",
    help = "File path to predicted eigengenes csv [default: %default]"
  ),
  make_option(
    c("-m", "--pcs"),
    type = "character",
    default = "genomic/mypcair.PCs.csv",
    help = "File path to PCs csv [default: %default]"
  ),
  make_option(
    c("-c", "--outcome"),
    type = "character",
    help = "String specifying the column in data containing the outcome variable [required]"     
  ),
  make_option(
    c("-b", "--binomial"),
    action = "store_true",
    default = FALSE,
    help = "Use binomial family for logistic regression [default: %default]"
  ),
  make_option(
    c("-f", "--factors"),
    default = "0,1",
    help = "when binomial == TRUE, what are the levels (ordered) for comparison [default: %default]"
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
  "preds",
  "pcs",
  "outcome",
  "binomial",
  "factors",
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
orig_outcome <- outcome
outcome <- make.names(outcome)
factors <- strsplit(factors, ",")[[1]]

# first_prot = prots[1]                                             
# last_prot = prots[length(prots)]

qc_cols <- c(
  "PIDN", "SampleID", "ADRC_ID", "Visit", "Visit_Consensus", "ConnectivityZscore",
  "Diagnosis_group", "DS1V"
)


# Functions ------------------


# Main -----------------------
# So, the final INTnorm dataset:
a <- fread(input, colClasses = c(PIDN = "character"))

# # Map in PCs assumes your primary data table exists 
# if (pcs != "") {
  # print(sprintf("Reading PCs for inclusion from %s", pcs))
  # x <- fread(pcs, colClasses = c(PIDN = "character"))
  # a <- merge(a, x, by = "PIDN")
  # covs <- unique(c(covs, names(a)[names(a) %like% "^PC\\d{1,2}$"])) # flexible for PCs
# }

if (preds != "") {
  # Load predicted
  x <- fread(preds)
  names(x)[-1] <- paste0(names(x)[-1], ".pred")
  a <- merge(a, x, by = "SampleID")
}

# shield improper outcome
setnames(a, 
  old = c(orig_outcome), 
  new = c(outcome)
)

a = a[!is.na(get(outcome))]
# cols = c(qc_cols, covs)

# Setup for baseline visit only analysis ---------------------------
# already done in the prep for baseline_clin

# generate surrogate variables ---------------------------------------
print("Generating surrogate variables...")

# matrix of proteins
cols <- c("SampleID", prots)
prot <- data.frame(a[, ..cols], row.names = "SampleID")

# matrix of phenotypes to regress
cols <- c("SampleID", outcome, covs)
pheno <- data.frame(a[, ..cols], row.names = "SampleID")
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

# SV x module VIF diagnostic --------------------------------------------
# For each module, treat it as the exposure alongside the SVs and other
# covariates, and compute each SV's VIF in that design (1/(1-R^2) from
# regressing the SV on every other term). High VIF for an SV across most
# modules -> safe to drop that SV globally for this phenotype's -v list.
# High VIF confined to a few modules -> those modules' results should be
# treated cautiously even with the SV retained.
orig_covs <- setdiff(covs, names(svs))
sv_names <- names(svs)

vif_mat <- matrix(NA_real_, nrow = length(prots), ncol = length(sv_names),
                   dimnames = list(prots, sv_names))

for (m in prots) {
  if (!m %in% names(a)) next
  dsub <- as.data.frame(a[, c(m, orig_covs, sv_names), with = FALSE])
  dsub <- dsub[complete.cases(dsub), ]
  # drop columns that are constant within this module's complete-case
  # subset -- model.matrix() errors trying to build contrasts for a
  # single-level factor
  dsub <- dsub[, sapply(dsub, function(x) length(unique(x)) > 1), drop = FALSE]
  if (!m %in% names(dsub) || nrow(dsub) < ncol(dsub) + 5) next  # not enough rows/variation to fit reliably

  mm <- model.matrix(~ ., data = dsub)[, -1, drop = FALSE]  # dummy-codes any factor covariates
  for (sv in sv_names) {
    if (!sv %in% colnames(mm)) next
    y <- mm[, sv]
    X <- mm[, setdiff(colnames(mm), sv), drop = FALSE]
    r2 <- tryCatch(summary(lm(y ~ X))$r.squared, error = function(e) NA_real_)
    vif_mat[m, sv] <- if (!is.na(r2) && r2 < 1) 1 / (1 - r2) else NA_real_
  }
}

fwrite(as.data.table(round(vif_mat, 3), keep.rownames = "module"),
  file = file.path(output, paste(nameset, "sv_vif_matrix.csv", sep = ".")))


# GLM based feature selection -----------------------------------------
print("Testing feature importance for covariates...")
cols <- c("SampleID", outcome, covs)

# drop missing values for feature importance testing
pheno <- data.frame(a[, ..cols], row.names = "SampleID")
pheno <- pheno[complete.cases(pheno),] 

# setup for phenos
if (binomial) {
  pheno[[outcome]] <- factor(pheno[[outcome]], levels = factors, labels = c("CTRL", "CASE"))
}

# formula
fmla <- reformulate(".", response = outcome)
# number of features
n_features <- ncol(pheno) - 1

# include all possible main effects
model = glm(fmla, 
  data = pheno,
  family = if (binomial) binomial() else gaussian()
) 

# simple training
set.seed(123)
cv_model <- train(
  form = fmla, 
  data = pheno, 
  method = "glm",
  family = if (binomial) binomial() else gaussian(),
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
# Conditional trControl based on outcome type
if (binomial) {
  tc <- trainControl(
    method = "cv", 
    number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
  metric_to_use <- "ROC"
} else {
  tc <- trainControl(method = "cv", number = 5)
  metric_to_use <- "RMSE"
}

cv_model_pls <- train(
  fmla, 
  data = pheno, 
  method = "simpls",
  trControl = tc,
  preProcess = c("zv", "center", "scale"),
  tuneLength = 15,
  metric = metric_to_use
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
# p5 <- vip::vip(cv_model_pls, num_features = 55, bar = FALSE) +
  # ggtitle(paste("Partial Least Squares - Best Fit ", npcs, "PCs"))


# Random Forest feature scoring ---------------------------------------

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

# gridExtra::grid.arrange(p1, p2, nrow = 1)

# feature interaction
# layout <- "
# AACCCCCC
# BBCCCCCC
# "
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

cols = c(qc_cols, covs, orig_outcome, prots)

fwrite(
  a[, ..cols], 
  file = file.path(output, paste(nameset, "WSSnorm.csv.gz", sep="."))
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

