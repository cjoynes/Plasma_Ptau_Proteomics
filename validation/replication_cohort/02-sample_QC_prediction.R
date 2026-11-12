#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 02-sample_QC_prediction.R
##
## Version: 0.0.2
##   v0.0.2 - INT proteins before SVA and the INTnorm file (not with --clin)
##
## Purpose of script: covariate/surrogate-variable selection for the
## replication cohort, before running the Figure 6 association tests
##
## Author: Robert B Butler III
##
## Date Created: 2026-08-10
##
## Copyright (c) 2026
## Email: rrbutler@stanford.edu
##
## ---------------------------
##
## Notes:
##   Adapted from 13-sample_QC_prediction.R for a single, cross-sectional
##   external cohort (one row per subject, no repeated visits or
##   genotype data). Differences from the internal script:
##     - No --eigengenes/--preds/-m/--pcs: the WSS-score and clinical
##       merge already happened once in 01-WSSmod_prep.Rmd, so this
##       script starts directly from that single merged file. --clin is
##       still needed here (unlike those) since prots means either the
##       full protein panel (Section A) or WSS modules (Section B), and
##       the SV x feature VIF diagnostic below is only meaningful for
##       the module case.
##     - No PIDN/DS1V/repeat-visit collapsing: every row is already one
##       subject.
##     - SVA is computed directly with sva(t(prot_matrix), mod, mod0,
##       n.sv=svs) on the full sample -- the internal project_svs()
##       train-visit -> all-visits projection step is unnecessary here
##       since there is no held-out visit to project onto.
##     - New --svs (integer, default 0): 0 skips SVA entirely. Surrogate
##       variables are optional -- see README.md for guidance on when to
##       use them and how to check for overfitting.
##
##   This script is meant to be run manually, once per outcome you plan
##   to test, BEFORE editing the covs/svs values in 00-commands.sh. Look
##   at the resulting {nameset}.feature_importance.pdf (GLM/PLS/Random
##   Forest variable importance) to decide which of your candidate -v
##   covariates (and, if used, how many surrogate variables) are worth
##   including.
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
library(tools)
library(stringr)

# Nameset helper function
nameset_default <- function(input) {
  file_path_sans_ext(basename(input))
}

# read in args
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "dataProt_replication.csv.gz",
    help = "File path to the merged intermediate file from 01-WSSmod_prep.Rmd [default: %default]"
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
    default = "prot_names.txt",
    help = "File path to the feature list used for surrogate variable estimation (prot_names.txt for Section A, WSS_names.txt/WSS_names_pred.txt for Section B), and carried through as the exposure list for 03-linreg_glm.R [default: %default]"
  ),
  make_option(
    c("-v", "--covs"),
    type = "character",
    default = "",
    help = "Candidate covariates to evaluate, comma-separated (e.g. Age,Gender or Age,Gender,CountryOfOrigin) [default: %default]"
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
    type = "character",
    default = "0,1",
    help = "When binomial == TRUE, what are the levels (ordered) for comparison [default: %default]"
  ),
  make_option(
    c("--svs"),
    type = "integer",
    default = 0,
    help = "Number of surrogate variables to compute via sva(). 0 skips SVA entirely -- surrogate variables are optional, see README.md [default: %default]"
  ),
  make_option(
    c("--clin"),
    action = "store_true",
    default = FALSE,
    help = "Running clinical phenotypes (WSS modules) as prots, not raw proteins? [default: %default]"
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
  "prots",
  "covs",
  "outcome",
  "binomial",
  "factors",
  "svs",
  "clin",
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

# check outpath status
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}

# variables
covs <- strsplit(covs, ",")[[1]]
prots <- readLines(prots)
orig_outcome <- outcome
outcome <- make.names(outcome)
factors <- strsplit(factors, ",")[[1]]


# Main -----------------------
a <- fread(input)

# shield improper outcome 
setnames(a, old = orig_outcome, new = outcome)
a <- a[!is.na(get(outcome))]

# Inverse rank normalization (INT) of proteins to z-scale ----------------
# Skipped with --clin, where prots are WSS modules that are already on an INT scale.
if (!clin) {
  int_cols <- intersect(prots, names(a))
  stopifnot("no columns in --prots found in the input" = length(int_cols) > 0)
  if (length(int_cols) < length(prots)) {
    print(sprintf("%d of %d prots not found in the input and not normalized",
      length(prots) - length(int_cols), length(prots)))
  }
  a[, (int_cols) := lapply(.SD, RankNorm), .SDcols = int_cols]
}

# surrogate variables (optional) ---------------------------------------
if (svs > 0) {
  print(sprintf("Generating %d surrogate variable(s)...", svs))
  print(sprintf("Covs: %s", paste(covs, collapse = ", ")))
  print(sprintf("Prots: %s...", paste(head(prots), collapse = ", ")))
  print(sprintf("Outcome: %s", outcome))

  # matrix of features (proteins for Section A, WSS module scores for
  # Section B) used to estimate surrogate variation
  cols <- c("SampleID", prots)
  prot <- data.frame(a[, ..cols], row.names = "SampleID")

  # matrix of phenotypes to regress
  cols <- c("SampleID", outcome, covs)
  pheno <- data.frame(a[, ..cols], row.names = "SampleID")
  stopifnot(all(rownames(prot) == rownames(pheno)))

  # Remove columns with only 1 unique value
  pheno_filtered <- pheno[, sapply(pheno, function(x) length(unique(x)) > 1)]
  print(
    sprintf("These covs were dropped for sv calculation as univariate: %s",
      paste(setdiff(names(pheno), names(pheno_filtered)), collapse = ", ")
    )
  )
  pheno <- pheno_filtered

  mod <- model.matrix( # includes all covs and the outcome
    as.formula(paste("~", paste(names(pheno), collapse = " + "))),
    data = pheno
  )
  mod0 <- model.matrix( # drops the outcome
    as.formula(
      paste("~", paste(names(pheno)[names(pheno) != outcome], collapse = " + "))
    ),
    data = pheno
  )

  # cross-sectional cohort: no train/held-out-visit split needed, so the
  # SVs are computed directly on the full sample (no project_svs() step)
  svobj <- sva(t(prot), mod, mod0, n.sv = svs)
  # sva() does not preserve rownames on svobj$sv -- restore them from
  # `prot` (same sample order) before keep.rownames can pick them up
  rownames(svobj$sv) <- rownames(prot)
  svs_dt <- data.table(svobj$sv, keep.rownames = "SampleID")
  setnames(svs_dt, setdiff(names(svs_dt), "SampleID"), paste0("sv", seq_len(svs)))
  a <- merge(a, svs_dt, by = "SampleID")

  sv_names <- paste0("sv", seq_len(svs))
  covs <- c(covs, sv_names)

  # SV x feature VIF diagnostic (--clin only; prots means the full raw
  # protein panel, not modules, in the Section A path, where this
  # diagnostic isn't meaningful and would be far more expensive) --------
  if (clin) {
    # For each candidate exposure in `prots` (WSS modules here, since
    # this only runs for --clin), treat it as the exposure alongside the
    # SVs and other covariates, and compute each SV's VIF in that design
    # (1/(1-R^2) from regressing the SV on every other term). High VIF for
    # an SV across most exposures -> safe to drop that SV globally; high
    # VIF confined to a few -> those exposures' results should be treated
    # cautiously even with the SV retained. No project_svs()/repeated-visit
    # collapsing needed here (see script header) -- SVs were already fit
    # directly on this cross-sectional, one-row-per-subject sample.
    orig_covs <- setdiff(covs, sv_names)
    vif_mat <- matrix(NA_real_, nrow = length(prots), ncol = length(sv_names),
                       dimnames = list(prots, sv_names))

    for (pr in prots) {
      if (!pr %in% names(a)) next
      dsub <- as.data.frame(a[, c(pr, orig_covs, sv_names), with = FALSE])
      dsub <- dsub[complete.cases(dsub), ]
      dsub <- dsub[, sapply(dsub, function(x) length(unique(x)) > 1), drop = FALSE]
      if (nrow(dsub) < ncol(dsub) + 5) next  # not enough rows to fit reliably

      mm <- model.matrix(~ ., data = dsub)[, -1, drop = FALSE]
      for (sv in sv_names) {
        if (!sv %in% colnames(mm)) next
        y <- mm[, sv]
        X <- mm[, setdiff(colnames(mm), sv), drop = FALSE]
        r2 <- tryCatch(summary(lm(y ~ X))$r.squared, error = function(e) NA_real_)
        vif_mat[pr, sv] <- if (!is.na(r2) && r2 < 1) 1 / (1 - r2) else NA_real_
      }
    }

    fwrite(as.data.table(round(vif_mat, 3), keep.rownames = "feature"),
      file = file.path(output, paste(nameset, "sv_vif_matrix.csv", sep = ".")))
  }
}


# GLM/PLS/Random Forest feature importance -----------------------------
print("Testing feature importance for covariates...")
cols <- c("SampleID", outcome, covs)

# drop missing values for feature importance testing
pheno <- data.frame(a[, ..cols], row.names = "SampleID")
pheno <- pheno[complete.cases(pheno), ]

# Remove columns with only 1 unique value
pheno_filtered <- pheno[, sapply(pheno, function(x) length(unique(x)) > 1)]
print(
  sprintf("These covs were dropped for feature importance calculation as univariate: %s",
    paste(setdiff(names(pheno), names(pheno_filtered)), collapse = ", ")
  )
)
pheno <- pheno_filtered

# setup for phenos
if (binomial) {
  pheno[[outcome]] <- factor(pheno[[outcome]], levels = factors, labels = c("CTRL", "CASE"))
}

# formula
fmla <- reformulate(".", response = outcome)
# number of features
n_features <- ncol(pheno) - 1

# simple training
set.seed(123)
cv_model <- train(
  form = fmla,
  data = pheno,
  method = "glm",
  family = if (binomial) binomial() else gaussian(),
  trControl = trainControl(method = "cv", number = 5)
)

# Partial Least Squares feature scoring -------------------------------
set.seed(123)
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

p1 <- vip_plot(cv_model, num_features = 50) +
  ggtitle("GLM - 5-fold Cross Validation") +
  geom_hline(yintercept = min(25, floor(length(cv_model$coefnames) / 2)))
p2 <- vip_plot(cv_model_pls,num_features = 50) +
  ggtitle(paste("Partial Least Squares - Best Fit ",
                cv_model_pls$bestTune$ncomp, "PCs")) +
  geom_hline(yintercept = min(25, floor(length(cv_model_pls$coefnames) / 2)))

# Random Forest feature scoring ---------------------------------------
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
for (i in seq_len(nrow(hyper_grid))) {
  fit.rf <- ranger(
    formula = fmla,
    data = pheno,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$min.node.size[i],
    replace         = hyper_grid$replace[i],
    sample.fraction = hyper_grid$sample.fraction[i],
    verbose         = FALSE,
    respect.unordered.factors = "order",
  )
  hyper_grid$rmse[i] <- sqrt(fit.rf$prediction.error)
}

# assess top 10 models
top <- hyper_grid %>%
  arrange(rmse) %>%
  mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
  head(10)

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

p3 <- vip_plot(rf_impurity, num_features = 50) +
  ggtitle("Random Forest - Impurity") +
  geom_hline(yintercept = min(25, floor(n_features / 2)))
p4 <- vip_plot(rf_permutation, num_features = 50) +
  ggtitle("Random Forest - Permutation") +
  geom_hline(yintercept = min(25, floor(n_features / 2)))

pdf(
  file = file.path(output, paste(nameset, "feature_importance.pdf", sep = ".")),
  width = 11,
  height = 11
)
p1 + p2 + p3 + p4
dev.off()

# Write the final INTnorm file for 03-linreg_glm.R ----------------------------
# shield improper outcome back to its original name
setnames(a, old = outcome, new = orig_outcome)
cols <- c("SampleID", covs, orig_outcome, prots)

fwrite(
  a[, ..cols],
  file = file.path(output, paste(nameset, "INTnorm.csv.gz", sep = "."))
)

cat(sprintf(
  "Wrote %s with covs: %s\n",
  file.path(output, paste(nameset, "INTnorm.csv.gz", sep = ".")),
  paste(covs, collapse = ", ")
))
