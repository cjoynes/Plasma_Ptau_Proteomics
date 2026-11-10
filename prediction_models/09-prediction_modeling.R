#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: 09-prediction_modeling.R
##
## Version: 0.0.1
##
## Purpose of script: Prepare Cell enrichment data for linear regression
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
library(auditor)
library(joinet)
library(randomForestSRC)

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
    default = "",
    help = "File path to PCs csv [default: %default]"
  ),
  make_option(
    c("-c", "--outcome"),
    type = "character",
    help = "String specifying the column in data containing the outcome variable [required]"     
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
  "outcome",
  "nameset"
)

input <- "C:/Users/audweber/Box/plasma_ptau_proteomics/dataProt_SS-205063.merged.QC.csv" 
output <- "prediction_models"
prots <- "C:/Users/audweber/Box/plasma_ptau_proteomics/prot_names.txt" 
covs <- 'Age,Gender,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,APOE'   
outcome <- c("Microglia","Astrocyte","Neuron","Oligodendrocyte") 
nameset <- "dataProt_SS-205063.PC.Plasma_pTau217.APOE" 
pcs <- "C:/Users/audweber/Box/plasma_ptau_proteomics/genomic/mypcair.PCs.csv"
prediction_scores <- "C:/Users/audweber/Box/plasma_ptau_proteomics/celltype_scores/dataProt_SS-205063.Sharma.lvl1.normalized_scores.csv"

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

# # dealing with duplicated column names
# drops <- which(duplicated(names(a)[!names(a) %in% prots]))
# a[, (drops) := NULL]

# Map in PCs assumes your primary data table exists 
if (pcs != "") {
  print(sprintf("Reading PCs for inclusion from %s", pcs))
  x <- fread(pcs, colClasses = c(PIDN = "character"))
  a <- merge(a, x, by = "PIDN")
  covs <- c(covs, names(a)[names(a) %like% "^PC\\d{1,2}$"]) # flexible for PCs
}

# Add cell enrichments

if (prediction_scores != "") {
  print(sprintf("Reading Prediction scores for inclusion from %s", prediction_scores))
  x <- fread(prediction_scores, colClasses = c(SampleID = "character"))
  a <- merge(a, x, by = "SampleID")
}
# Initial QC filters
# a = a[SampleType == "Sample"]
# a = a[Study != "Control"]
# a = a[RowCheck == "PASS"]
a = a[!is.na(get(outcome))]
# # filter out missing diagnosis, PD people
# pd_folks = c(unique(a$Diagnosis_group[a$Diagnosis_group %like% "PD"]), "")
# lbd_folks = c(unique(a$Diagnosis_group[a$Diagnosis_group %like% "LBD"]), "")
# a = a[!Diagnosis_group %in% c(pd_folks, lbd_folks)]
# # # Storage days has a few NAs
# # recode Subarray to letters to avoid confusing as quantitative
# a[, Subarray := recode(Subarray, `1`="A", `2`="B", `3`="C", `4`="D", `5`="E", 
#                        `6`="F", `7`="G", `8`="H")]

cols = c(qc_cols, covs)

# Sample outlier analysis --------------------------------------------------
# check ConnectivityZscore for extreme outliers:
# a[abs(ConnectivityZscore) > 3, ..cols]
# 
# # plots connectivity histogram
# p1 = ggplot(a, aes(x=ConnectivityZscore)) + geom_histogram()
# ggsave(
#   paste(nameset, "ConnectivityZscore.pdf", sep = "."), 
#   p1, 
#   useDingbats = F, 
#   path = output
# )
# 
# # One sample had an extreme Z-score, less than -23. Participant's baseline visit
# was in a normal range, and there was no other metadata explaining the outlier.
# The rest included some samples with z in the +/- 3-4 range, but not notable

# fetch the protein columns for plotting
# cols = c("SampleID", prote) 
# b = c("SampleID", "normalized_scores", a)
# 
# # Poisson Distance Pairwise
# poisd = PoissonDistance(as.matrix(b, rownames="SampleID"))
# samplePoisDistMatrix = as.matrix(poisd$dd)
# rownames(samplePoisDistMatrix) = b[, SampleID]
# colnames(samplePoisDistMatrix) = b[, SampleID]
# 
# # set colors and draw
# colors = colorRampPalette(brewer.pal(11, "RdYlBu"))(255)
# pheatmap(samplePoisDistMatrix, col=colors, clustering_distance_rows=poisd$dd,
#          clustering_distance_cols=poisd$dd, width=110, height=85, border_color=NA,
#          filename=file.path(output, paste(nameset, "Poisson_distance.pdf", sep=".")))
# 
# # Some groupings that could be interesting, one major outlier is the same one
# from z-score connectivity, final check PCA

# PCA
# Function ---------------------------
# PCA top X genes by variance
# c_by = column to color by
# c_lab = Label for color legend
# c_type = "disc" discrete or "cont" continuous coloring
# num = number of genes, default is all
# my_pca_plot = function(vst, c_by, c_lab, c_type="disc", num=1000000) {
#   rv = rowVars(vst)
#   top_genes = order(rv, decreasing=T)[seq_len(min(num, length(rv)))]
#   pca = prcomp(t(vst[top_genes, ]))
#   percent_pc = round((pca$sdev^2 / sum(pca$sdev^2))*100, 2)
#   if (c_type == "disc") {
#     col_by = as.factor(a[[c_by]])
#   } else if (c_type == "cont") {
#     col_by = a[[c_by]]
#   } else {
#     print("Incorrect value for c_type!")
#     quit("no", 1)
#   }
# }
#   
#   # PC1 v PC2
#   ggPCA12 = ggplot(data.frame(pca$x), aes(x=PC1, y=PC2, col=col_by)) +
#     labs(x=sprintf("PC1 - %.2f%%", percent_pc[1]),
#          y=sprintf("PC2 - %.2f%%", percent_pc[2]), col=c_lab) +
#     geom_point(size=3, alpha=0.5) +
#     geom_text_repel(aes(label=rownames(data.frame(pca$x)))) +
#     ggtitle(sprintf("top %s genes", length(top_genes))) + theme_classic() + {
#       if(c_type == "disc") scale_color_brewer(type="qual", palette="Paired") 
#       else scale_color_distiller(type="div", palette="RdYlBu")
#     }
#   
#   # PC1 v PC3
#   ggPCA13 = ggplot(data.frame(pca$x), aes(x=PC1, y=PC3, col=col_by)) +
#     labs(x=sprintf("PC1 - %.2f%%", percent_pc[1]),
#          y=sprintf("PC3 - %.2f%%", percent_pc[3]), col=c_lab) +
#     geom_point(size=3, alpha=0.5) +
#     geom_text_repel(aes(label=rownames(data.frame(pca$x)))) +
#     ggtitle(sprintf("top %s genes", length(top_genes))) + theme_classic() + {
#       if(c_type == "disc") scale_color_brewer(type="qual", palette="Paired") 
#       else scale_color_distiller(type="div", palette="RdYlBu")
#     }
#   plot(ggPCA12)
#   plot(ggPCA13)
# }
# 
# # Plot PCAs
# pdf(file=file.path(output, paste(nameset, "pca_plots.pdf", sep=".")), onefile=T, paper="USr", 
#     width=11, height=8.5)
# my_pca_plot(t(as.matrix(b, rownames="SampleID")), c_by="Diagnosis_group", 
#             c_lab="Diagnosis_group", c_type="disc")
# my_pca_plot(t(as.matrix(b, rownames="SampleID")), c_by="Diagnosis_group", 
#             c_lab="Top 500 Diagnosis_group", c_type="disc", num=500)
# dev.off()

# Same outlier sample in PCs as well, same across three different metrics, remove
# drop outlier samples
drop = a[abs(ConnectivityZscore) > 10, SampleID]
a = a[!SampleID %in% drop]


# Setup for baseline visit only analysis ---------------------------
# # # drop all samples with DS1V != 0
# # a = a[DS1V == 0]
# Correction, keep only first visit where plasma available
a = a[order(DS1V)][, head(.SD, 1), by = "PIDN"]

# checking for missed duplicates or bad visit labeling
a[duplicated(PIDN), ..qc_cols]
# 
# # Inverse rank normalization (INT) to z-scale
# cols = c(outcome, prots)
# a[, (cols) := lapply(.SD, RankNorm), .SDcols=cols]

# handling incomplete cases
cols <- c("SampleID", outcome, covs)
drops <- a[!complete.cases(a[, ..cols]), SampleID]
a <- a[!SampleID %in% drops]

# # generate surrogate variables ---------------------------------------
# # matrix of proteins
# cols <- c("SampleID", prots)
# prot <- t(as.matrix(a[, ..cols], rownames = "SampleID"))
# 
# # matrix of phenotypes to regress
# cols <- c("SampleID", outcome, covs)
# pheno <- data.frame(a[, ..cols], row.names = "SampleID")
# all(colnames(prot) == rownames(pheno))
# 
# # defining models with static set of 6 SVs
# mod <- model.matrix( # this will include all covs and the outcome
#   as.formula(paste("~", paste(names(pheno), collapse=" + "))), 
#   data = pheno
# ) 
# mod0 <- model.matrix( # this will drop the outcome
#   as.formula(
#     paste("~", paste(names(pheno)[names(pheno) != outcome], collapse = " + "))
#   ), 
#   data = pheno
# )
# svobj <- sva(prot, mod, mod0, n.sv = 6)
# 
# # extract and add to a
# svs <- data.table(svobj$sv)
# setnames(svs, gsub("V", "sv", names(svs)))
# a <- cbind(a, svs)
# 
# # add to covs 
# covs <- c(covs, names(svs))
# cols <- c("SampleID", outcome, covs)

# GLM based feature selection -----------------------------------------
# examined potential covariates (use pheno data.frame above)
# drop missing values for feature importance testing
pheno <- data.frame(a[, ..cols], row.names = "SampleID")
pheno <- pheno[complete.cases(pheno),] 
pheno <- pheno[1:12]

# include all possible main effects
# model = glm(as.formula(paste(outcome, "~", ".")), data = pheno) 
# 
# # simple training
# set.seed(123)
# cv_model <- train(
#   as.formula(paste(outcome, "~", ".")), 
#   data = pheno, 
#   method = "glm",
#   trControl = trainControl(method = "cv", number = 50)
# )

# Create a tuning grid that evaluates all three models
# Ridge: alpha = 0
# Elastic Net: alpha = 0.5
# Lasso: alpha = 1

y_cols <- c("Microglia","Astrocyte","Neuron","Oligodendrocyte")
feature_sets <- list(
  full_biomarkers   = setdiff(names(pheno2), c(y_cols, "SampleID")),  # everything except outcomes/ID
  selected_markers  = c("PlasmaAB142P","PlasmaAB140P","PlasmaABRatio","PlasmapTau217",
                        "PlasmapTau217_AB42Ratio","PlasmaGFAP","PlasmaNfL")                 
)
tuning_grid <- expand.grid(
  alpha = c(0, 0.1, 0.2, 0.3, 0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1),  
  lambda = seq(0.001, 0.1, length = 20)  # Test 10 lambda values
)



# Run cross-validation with all three methods
cv_model <- train(
  as.formula(paste(outcome, "~", ".")), 
  data = pheno,
  method = "glmnet",
  trControl = trainControl(
    method = "cv", 
    number = 50,
    savePredictions = "final",
    verboseIter = TRUE
  ),
  tuneGrid = tuning_grid,  # Use custom grid
  metric = "RMSE"  # Or "Accuracy" for classification
)

print(cv_model)
print(cv_model$bestTune) 


coefs <- coef(cv_model2$finalModel, s = cv_model2$bestTune$lambda)
# Visualize the results across methods- only works when grid is set to 0,0.5,1
ggplot(cv_model$results, aes(x = lambda, y = RMSE, color = factor(alpha))) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("blue", "purple", "red"),
                     labels = c("Ridge (α=0)", "Elastic Net (α=0.5)", "Lasso (α=1)"),
                     name = "Method") +
  labs(title = "Performance Comparison of Regularization Methods",
       x = "Lambda (regularization strength)",
       y = "RMSE") +
  theme_minimal()

# Save the best model
saveRDS(cv_model, file.path(output, paste0(args_list$nameset, "_best_glmnet_model.rds")))

# Create separate models for each method 
# Ridge
ridge_model <- train(
  as.formula(paste(outcome, "~", ".")), 
  data = pheno,
  method = "glmnet",
  trControl = trainControl(method = "cv", number = 50),
  tuneGrid = expand.grid(alpha = 0, lambda = seq(0.001, 0.1, length = 10))
)

# Elastic Net
elastic_model <- train(
  as.formula(paste(outcome, "~", ".")), 
  data = pheno,
  method = "glmnet",
  trControl = trainControl(method = "cv", number = 50),
  tuneGrid = expand.grid(alpha = 0.5, lambda = seq(0.001, 0.1, length = 10))
)

# Lasso
lasso_model <- train(
  as.formula(paste(outcome, "~", ".")), 
  data = pheno,
  method = "glmnet",
  trControl = trainControl(method = "cv", number = 50),
  tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 0.1, length = 10))
)

# Compare the three models
models <- list(Ridge = ridge_model, ElasticNet = elastic_model, Lasso = lasso_model)
resamples <- resamples(models)
summary(resamples)

# Visualize comparison
bwplot(resamples, metric = "RMSE")

## Run model again including features (our features will be the plasma biomarkers)
outcome <- "Microglia"
feature <- c(
  "PlasmaPTau181", "PlasmaAB142P","PlasmaAB140P",  
  "PlasmaABRatio", "PlasmapTau217","PlasmapTau217_AB42Ratio",
  "PlasmaGFAP","PlasmaNfL" 
)
col2 <- c("SampleID", outcome, covs, feature)
pheno2 <- data.frame(a[, ..col2], row.names = "SampleID")
pheno2 <- pheno2[complete.cases(pheno2),] 
pheno2 <- pheno2[,c(2)]


cv_model2 <- train(
  as.formula(paste(outcome, "~", ".")), 
  data = pheno2,
  method = "glmnet",
  trControl = trainControl(
    method = "cv", 
    number = 50,
    savePredictions = "final",
    verboseIter = TRUE
  ),
  tuneGrid = tuning_grid,  # Use custom grid
  metric = "RMSE"  # Or "Accuracy" for classification
)

#Best lamda and alpha 
print(cv_model2)
print(cv_model2$bestTune)


# coefficients
cm_micro <- coef(cv_model2$finalModel, s = cv_model2$bestTune$lambda)

micro_model <- cv_model2
# Microglia: alpha     lambda
#             0.1 0.04789474

#Astroctyes alpha     lambda
#            0.7 0.02705263
#Neurons :    alpha     lambda
#                1 0.01663158
#Olgos   alpha lambda
#          0.6    0.1
# Residual error analysis 
cv_exp <- DALEX::explain(cv_model2, label = "cv", data = pheno2, y = pheno2$Astrocytes)
ridge_exp <- DALEX::explain(ridge_model, label = "ridge", data = pheno2, y = pheno2$Microglia)
lasso_exp <- DALEX::explain(lasso_model, label = "lasso", data = pheno2, y = pheno2$Microglia)
elastic_exp <- DALEX::explain(elastic_model, label = "elastic", data = pheno2, y = pheno2$Microglia)

#Model residuals 
cv_mr <- model_residual(cv_exp)
ridge_mr <- model_residual(ridge_exp)
lasso_mr <- model_residual(lasso_exp)
elastic_mr <- model_residual(elastic_exp)

#plots 
plot(cv_mr, type = "prediction", abline = TRUE)
plot(cv_mr, cv_mr, ridge_mr, lasso_mr, elastic_mr, type = "rec")
plot(cv_mr, type = "rroc")

# Visualize the results across methods
ggplot(cv_model2$results, aes(x = lambda, y = RMSE, color = factor(alpha))) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("blue", "purple", "red"),
                     labels = c("Ridge (α=0)", "Elastic Net (α=0.5)", "Lasso (α=1)"),
                     name = "Method") +
  labs(title = "Performance Comparison of Regularization Methods",
       x = "Lambda (regularization strength)",
       y = "RMSE") +
  theme_minimal()

# Save the best model
saveRDS(cv_model, file.path(output, paste0(args_list$nameset, "_best_glmnet_model.rds")))

# Create separate models for each method 
# Ridge
ridge_model2 <- train(
  as.formula(paste(outcome, "~", ".")), 
  data = pheno2,
  method = "glmnet",
  trControl = trainControl(method = "cv", number = 50),
  tuneGrid = expand.grid(alpha = 0, lambda = seq(0.001, 0.1, length = 10))
)

# Elastic Net
elastic_model2 <- train(
  as.formula(paste(outcome, "~", ".")), 
  data = pheno2,
  method = "glmnet",
  trControl = trainControl(method = "cv", number = 50),
  tuneGrid = expand.grid(alpha = 0.5, lambda = seq(0.001, 0.1, length = 10))
)

# Lasso
lasso_model2 <- train(
  as.formula(paste(outcome, "~", ".")), 
  data = pheno2,
  method = "glmnet",
  trControl = trainControl(method = "cv", number = 50),
  tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 0.1, length = 10))
)

# Compare the three models
models2 <- list(Ridge = ridge_model2, ElasticNet = elastic_model2, Lasso = lasso_model2)
resamples2 <- resamples(models2)
summary(resamples2)

# Visualize comparison
bwplot(resamples2, metric = "RMSE")

###### Comparing different variates against each other 
### NOTE: Index is important because it allows us to compare across conditions
######    by using the same folds for each run 

y_cols <- c("Microglia","Astrocyte","Neuron","Oligodendrocyte")  # outcomes

# Variants: by default treat entries as "DROP these columns".
# Special case: "Only_biomarkers" means KEEP ONLY these predictors.
variants <- list(
  ALL             = character(0),               # keep all predictors
  NO_APOE         = c("APOE"),
  NO_APOE_Gender  = c("APOE","Gender"),
  NO_PCS          = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"),
  Only_biomarkers = c("PlasmaPTau181","PlasmaAB142P","PlasmaAB140P","PlasmaABRatio",
                      "PlasmapTau217","PlasmapTau217_AB42Ratio","PlasmaGFAP","PlasmaNfL")
)

# Pool of predictors = everything except the outcome columns
predictor_pool <- setdiff(names(pheno2), y_cols)

# glmnet tuning grid
tuning_grid <- expand.grid(
  alpha  = seq(0, 1, 0.1),
  lambda = 10^seq(-4, 0, length.out = 50)
)

# Output folder
out_dir <- "caret_runs"
dir.create(out_dir, showWarnings = FALSE)

# =========================
# 1) FIX CV FOLDS (reproducible)
# =========================
set.seed(42)
K <- 50
folds_train <- createFolds(seq_len(nrow(pheno2)), k = K, returnTrain = TRUE)
index    <- folds_train
indexOut <- lapply(folds_train, function(tr) setdiff(seq_len(nrow(pheno2)), tr))

# =========================
# 2) LOOP: outcome × variant
# =========================
rows <- list()
coef_txt <- file.path(out_dir, "coefficients_printouts.txt")
if (file.exists(coef_txt)) file.remove(coef_txt)  # fresh file each run

for (outcome in y_cols) {
  for (tag in names(variants)) {
    
    # Build predictor set for this variant
    if (tag == "Only_biomarkers") {
      # KEEP ONLY these biomarker columns
      keep <- intersect(variants[[tag]], predictor_pool)
      if (length(keep) == 0)
        warning(sprintf("[%s|%s] None of the requested 'Only_biomarkers' found in data.", outcome, tag))
      x_cols <- keep
    } else {
      # Default: DROP these columns
      drops <- variants[[tag]]
      missing_drops <- setdiff(drops, names(pheno2))
      if (length(missing_drops))
        warning(sprintf("[%s|%s] Drop list contains missing columns: %s",
                        outcome, tag, paste(missing_drops, collapse=", ")))
      x_cols <- setdiff(predictor_pool, drops)
    }
    
    # Data for caret
    df <- pheno2[, c(outcome, x_cols), drop = FALSE]
    
    # Train
    cv_model <- train(
      as.formula(paste(outcome, "~ .")),
      data = pheno2,
      method = "glmnet",
      trControl = trainControl(method = "cv", number = K,
                               index = index, indexOut = indexOut,
                               savePredictions = "final", verboseIter = TRUE),
      tuneGrid = tuning_grid,
      metric = "RMSE"
    )
    
    # Best tune + OOF metrics
    bt  <- cv_model$bestTune
    oof <- subset(cv_model$pred, alpha == bt$alpha & lambda == bt$lambda)
    rmse_val <- RMSE(oof$pred, oof$obs)
    r2_val   <- R2(oof$pred,   oof$obs)
    
    rows[[length(rows) + 1]] <- data.frame(
      outcome = outcome,
      variant = tag,
      best_alpha = bt$alpha,
      best_lambda = bt$lambda,
      rmse_oof = rmse_val,
      r2_oof   = r2_val,
      n_features = length(x_cols),
      stringsAsFactors = FALSE
    )
    
    # === COEFFICIENTS (exactly like coef()) ===
    co <- coef(cv_model$finalModel, s = bt$lambda)  # sparse matrix (".": zero)
    cat(sprintf("\n===== %s | %s =====\n", outcome, tag))
    print(co)  # console
    # also append the same printout to a text file
    cat(sprintf("\n===== %s | %s =====\n", outcome, tag), file = coef_txt, append = TRUE)
    capture.output(print(co), file = coef_txt, append = TRUE)
    
    # OPTIONAL: also save non-zero coefficients to CSV (tidy)
    cz <- as.matrix(co)
    coef_df <- data.frame(feature = rownames(cz), coef = as.numeric(cz), row.names = NULL)
    coef_df_nz <- subset(coef_df, coef != 0)  # drop zeros; keep intercept if you want
    coef_df_nz$outcome <- outcome
    coef_df_nz$variant <- tag
    fwrite(coef_df_nz,
           file = file.path(out_dir, "coefficients_nonzero.csv"),
           append = file.exists(file.path(out_dir, "coefficients_nonzero.csv")))
  }
}

# =========================
# 3) SAVE METRICS SUMMARY
# =========================
results <- rbindlist(rows)
fwrite(results, file.path(out_dir, "glmnet_outcomes_by_variant.csv"))

# Quick macro comparison (lower is better)
macro <- results[, .(macro_rmse = mean(rmse_oof), macro_r2 = mean(r2_oof)), by = variant][order(macro_rmse)]
fwrite(macro, file.path(out_dir, "macro_by_variant.csv"))

cat("\nSaved:\n",
    "- metrics: ", file.path(out_dir, "glmnet_outcomes_by_variant.csv"), "\n",
    "- macro:   ", file.path(out_dir, "macro_by_variant.csv"), "\n",
    "- coefs:   ", file.path(out_dir, "coefficients_printouts.txt"),
    " (raw coef() printouts)\n",
    "- coefs:   ", file.path(out_dir, "coefficients_nonzero.csv"),
    " (tidy non-zero)\n", sep = "")

## RROC plots 

elastic_exp <- DALEX::explain(elastic_model2, label = "elastic", data = pheno2, y = pheno2$Microglia)
lasso_exp <- DALEX::explain(lasso_model2, label = "lasso", data = pheno2, y = pheno2$Microglia)
ridge_exp <- DALEX::explain(ridge_model2, label = "ridge", data = pheno2, y = pheno2$Microglia)
rf_exp <- DALEX::explain(rf_model, label = "rf", data = pheno2, y = pheno2$Microglia)


elastic_mr <- model_residual(elastic_exp)
lasso_mr <- model_residual(lasso_exp)
ridge_mr <- model_residual(ridge_exp)
rf_mr <- model_residual(rf_exp)

plot(lasso_mr, ridge_mr, type = "prediction", abline = TRUE)
plot(cv_mr, rf_mr, type = "prediction", abline = TRUE)
plot(cv_mr, ridge_mr, lasso_mr, elastic_mr, type = "prediction", abline = TRUE)
plot(cv_mr, cv_mr, ridge_mr, lasso_mr, elastic_mr, type = "rec")
plot(cv_mr, cv_mr, ridge_mr, lasso_mr, elastic_mr, type = "rroc")

plot(cv_mr, rf_mr, type = "rec")
plot(cv_mr, rf_mr, type = "rroc")

##Test predictions

y_hat <- predict(cv_model2, newdata = pheno2)
resid <- pheno2[[outcome]] - y_hat

est <- cv_model2$bestTune
oof <- cv_model2$pred %>%
  filter(alpha == best$alpha, lambda == best$lambda) %>%
  arrange(rowIndex) %>%
  rename(y_hat = pred, y = obs)
rmse_oof <- RMSE(oof$y_hat, oof$y)
mae_oof  <- MAE(oof$y_hat,  oof$y)
r2_oof   <- R2(oof$y_hat,   oof$y)


pheno2_with_oof <- pheno2 %>%
  mutate(rowIndex = row_number()) %>%
  inner_join(oof %>% select(rowIndex, y_hat), by = "rowIndex") %>%
  mutate(resid = .data[[outcome]] - y_hat)
ggplot(pheno2_with_oof, aes(x = y_hat, y = .data[[outcome]])) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "predicted", y = "observed", title = "Observed vs Predicted")

# summary(cv_model$finalModel)


##Multivariate elastic net regression ----------------------------------
feature <- c(
  "PlasmaPTau181", "PlasmaAB142P","PlasmaAB140P",  
  "PlasmaABRatio", "PlasmapTau217","PlasmapTau217_AB42Ratio",
  "PlasmaGFAP","PlasmaNfL" 
)
col2 <- c("SampleID", outcome, covs, feature)
pheno2 <- data.frame(a[, ..col2], row.names = "SampleID")
pheno2 <- pheno2[complete.cases(pheno2),] 
pheno2 <- pheno2[,c(1:15,29:31)]

#All variables need to be numeric 
pheno2$Gender <- as.numeric(pheno2$Gender == "M")
pheno2$APOE <- ifelse(grepl("4", pheno2$APOE), 1, 0) 

#cv.joinet will find the best MSE 
cv_joinet_model <- cv.joinet(
  Y = as.matrix(pheno2[, c(1:4)]),
  X = as.matrix(pheno2[, c(5:23)]),
  family = "gaussian",
  type.measure = "mse",
  alpha.base = 0.1,
  alpha.meta = 0.2
)

##Testing multiple alphas 

# Inputs
idx <- c(5:13)
Y <- as.matrix(pheno2[, 1:4, drop = FALSE])
X <- as.matrix(pheno2[,idx, drop = FALSE])
if (is.null(colnames(X))) colnames(X) <- paste0("feat", seq_len(ncol(X)))
 
n <- nrow(Y)
alphas <- seq(0, 1, 0.1)
 
# 1) FIX CV SPLITS so every grid point uses the SAME folds
   set.seed(42)
   foldid.ext <- sample(rep(1:5, length.out = n))
   foldid.int <- sample(rep(1:10, length.out = n))

# 2) Run grid; compute macro-average MSE (rowMeans), keep mapping
grid <- CJ(alpha.base = alphas, alpha.meta = alphas)
res_list <- vector("list", nrow(grid))
 
  for (k in seq_len(nrow(grid))) {
         ab <- grid$alpha.base[k]; am <- grid$alpha.meta[k]
         m <- cv.joinet(
               Y, X, family = "gaussian", type.measure = "mse",
               alpha.base = ab, alpha.meta = am,
               nfolds.ext = 5, nfolds.int = 10,
               foldid.ext = foldid.ext, foldid.int = foldid.int
           )
#m: rows = c("none","base","meta"), cols = outcomes
           res_list[[k]] <- data.table(
                 alpha.base = ab,
                 alpha.meta = am,
                 method = rownames(m),
                 mse = rowMeans(m)   # macro-average across outcomes
             )
}

sums_wide <- dcast(rbindlist(res_list), alpha.base + alpha.meta ~ method, value.var = "mse")
sums_wide[, `:=`(
       diff = base - meta,           # >0 means META better than BASE (lower MSE)
       base.rmse = sqrt(base),
       meta.rmse = sqrt(meta),
       none.rmse = sqrt(none)
  )]
 

sums_marked <- copy(sums_wide)
sums_marked[, is_best := FALSE]
sums_marked[which.min(meta), is_best := TRUE]
setorder(sums_marked, -is_best, meta)   # best row on top, then by meta

# single CSV containing everything, with best flagged
one_csv <- file.path(output, "alpha_grid_with_best_nobiomakers.csv")
fwrite(sums_marked, one_csv)

  
#Run the model 
  joinet_model <- joinet(
    Y = as.matrix(pheno2[, c(1:4)]),
    X = as.matrix(pheno2[, c(5:23)]),
    family = "gaussian",
    nfolds = 50,
    foldid = NULL,
    type.measure = "mse",
    alpha.base = 1,
    alpha.meta = 1,
    weight = NULL,
    sign = NULL, 
    compare = TRUE
  )
  

  
# Random Forest feature scoring ---------------------------------------
# number of features
  n_features <- ncol(pheno2) - 1
  
  # train a default random forest model
  rf1 <- ranger(
    as.formula(paste(outcome, "~", ".")), 
    data = pheno2,
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
    num.trees = 1000,
    rmse = NA                                               
  )
  
  
  # execute full cartesian grid search
  for(i in seq_len(nrow(hyper_grid))) {
    # fit model for ith hyperparameter combination
    fit.rf <- ranger(
      formula = as.formula(paste(outcome, "~", ".")), 
      data = pheno2,
      num.trees        = hyper_grid$num.trees[i],
      mtry            = hyper_grid$mtry[i],
      min.node.size   = hyper_grid$min.node.size[i],
      replace         = hyper_grid$replace[i],
      sample.fraction = hyper_grid$sample.fraction[i],
      verbose         = FALSE,
      respect.unordered.factors = 'order',
      seed= 42,
    )
    # export OOB error 
    hyper_grid$rmse[i] <- sqrt(fit.rf$prediction.error)
  }
  
  # assess top 10 models
  top = hyper_grid %>%
    arrange(rmse) %>%
    mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
    head(10)
  # For Neuron:
# mtry min.node.size replace sample.fraction      rmse perc_gain
# 22             5   FALSE            0.63 0.9101213 1.2430475
## At seed 123:
# mtry min.node.size replace sample.fraction num.trees      rmse perc_gain
# 1     8            10   FALSE             0.5      1000 0.2090959 1.0440298
# 2     8             5    TRUE             0.5      1000 0.2091018 1.0412515
# #at seed 42: 
# mtry min.node.size replace sample.fraction num.trees      rmse perc_gain
# 1     8             5    TRUE            0.50      1000 0.2086859 1.2380649
# 2     8            10    TRUE            0.50      1000 0.2086906 1.2358362

  # For Microglia:
  # mtry min.node.size replace sample.fraction num.trees      rmse perc_gain
  # 1     6             5    TRUE            0.50      1000 0.2067569  2.150978
  # 2     6            10    TRUE            0.63      1000 0.2069148  2.076281
  # 3     6            10    TRUE            0.50      1000 0.2069324  2.067929
  # # For Astrocytes:  
  # mtry min.node.size replace sample.fraction num.trees      rmse perc_gain
  # 1     6            10    TRUE            0.50      1000 0.2214970  1.314299
  # 2     6             5    TRUE            0.50      1000 0.2217411  1.205556
  # # For Olgo:
  # mtry min.node.size replace sample.fraction num.trees      rmse perc_gain
  # 1     2             1    TRUE            0.50      1000 0.3513820 -56.55463
  # 2     2             3    TRUE            0.50      1000 0.3516030 -56.65312
  
  
  # re-run model with impurity-based variable importance
  rf_impurity <- ranger(
    formula = as.formula(paste(outcome, "~", ".")), 
    data = pheno2,
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
    data = pheno2,
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
    ggtitle("Random Forest - Impurity") + geom_hline(yintercept = 20)
  p4 <- vip_plot(rf_permutation) +
    ggtitle("Random Forest - Permutation") + geom_hline(yintercept = 20)
  
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
  

out_dir <- file.path("prediction_model", "vip plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

outcome <- "Oligodendrocyte"  # change to cell type your interested in 
col2 <- c("SampleID", outcome, covs, feature)
pheno2 <- data.frame(a[, ..col2], row.names = "SampleID")
pheno2 <- pheno2[complete.cases(pheno2),] 
pheno2 <- pheno2[,c(1:15, 24:31)]

  
covar_variants <- list(
    ALL            = character(0),                 # use all predictors
    NO_APOE        = c("APOE"),
    NO_APOE_SEX    = c("APOE","Gender"),
    NO_PCS         = c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8"),
    BIOMARKERS_ONLY= c("PlasmaPTau181","PlasmaAB142P","PlasmaAB140P","PlasmaABRatio",
                       "PlasmapTau217","PlasmapTau217_AB42Ratio","PlasmaGFAP","PlasmaNfL")
)
  
  # tuned hyperparams
  hp <- list(num.trees=1500, mtry=2, min.node.size=1, sample.fraction=0.50, replace=TRUE)
  
  predictor_pool <- setdiff(names(pheno2), outcome)
  
  # ---- helper: build df safely for each variant ----
  .build_df <- function(ph, outcome, drops_or_keep, keep_only = FALSE) {
    pool <- setdiff(names(ph), outcome)
    if (keep_only) {
      keep <- intersect(drops_or_keep, pool)
      if (!length(keep)) stop("Keep-only set left 0 predictors; check names.")
      ph[, c(outcome, keep), drop = FALSE]
    } else {
      drops <- intersect(drops_or_keep, names(ph))
      xcols <- setdiff(pool, drops)
      if (!length(xcols)) stop("Drop set removed all predictors; check names.")
      ph[, c(outcome, xcols), drop = FALSE]
    }
  }
  
  run_one <- function(tag, drops_or_keep, keep_only=FALSE) {
    df <- .build_df(pheno2, outcome, drops_or_keep, keep_only)
    set.seed(42)
    fit <- ranger(
      as.formula(paste(outcome, "~ .")), data=df,
      num.trees=hp$num.trees, mtry=hp$mtry, min.node.size=hp$min.node.size,
      replace=hp$replace, sample.fraction=hp$sample.fraction,
      importance="permutation", scale.permutation.importance=TRUE,
      respect.unordered.factors="partition", seed=42
    )
    rmse_oob <- sqrt(fit$prediction.error)
    vi <- sort(fit$variable.importance, decreasing=TRUE)
    list(tag=tag, p=ncol(df)-1, rmse_oob=rmse_oob, vi=vi, fit=fit)
  }
  
  # ---- run variants ----
  results <- list(
    run_one("ALL",             covar_variants$ALL,             keep_only=FALSE),
    run_one("NO_APOE",         covar_variants$NO_APOE,         keep_only=FALSE),
    run_one("NO_APOE_SEX",     covar_variants$NO_APOE_SEX,     keep_only=FALSE),
    run_one("NO_PCS",          covar_variants$NO_PCS,          keep_only=FALSE),
    run_one("BIOMARKERS_ONLY", covar_variants$BIOMARKERS_ONLY, keep_only=TRUE)
  )
  
  # ---- summary table (unchanged) ----
  tbl <- rbindlist(lapply(results, function(x)
    data.table(variant=x$tag, p=x$p, rmse_oob=x$rmse_oob)))
  tbl[order(rmse_oob)]
  
  fwrite(tbl, file.path(out_dir, sprintf("vip_rmse_summary_%s.csv", outcome)))
  
  
  # ---- VIP plots for each variant ----
  for (res in results) {
    p <- vip_plot(res$fit) +
      ggtitle(sprintf("VIP: %s | %s", outcome, res$tag)) +
      geom_hline(yintercept = 20)
    print(p)  # show in the viewer
    ggsave(file.path("prediction_model", "vip plots", sprintf("vip_%s_%s.png", outcome, res$tag)),
           p, width=7, height=5, dpi=300)
  }
  

## Multivariate Random Forest #################
  outcome <- c("Microglia","Astrocyte","Neuron","Oligodendrocyte")
  
#Fit for all covariates
  fit_all <- rfsrc(Multivar(Microglia, Astrocyte, Neuron, Oligodendrocyte) ~ 
          APOE + Gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
          Age + PlasmaAB142P + PlasmaAB140P + PlasmaABRatio + PlasmapTau217 +          
          PlasmapTau217_AB42Ratio + PlasmaGFAP +PlasmaNfL, 
        data=pheno2, importance=TRUE)
  fit_no_APOE <- rfsrc(Multivar(Microglia, Astrocyte, Neuron, Oligodendrocyte) ~ 
                     Gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + 
                     Age + PlasmaAB142P + PlasmaAB140P + PlasmaABRatio + PlasmapTau217 +          
                     PlasmapTau217_AB42Ratio + PlasmaGFAP +PlasmaNfL, 
                   data=pheno2, importance=TRUE)
  fit_no_PC <- rfsrc(Multivar(Microglia, Astrocyte, Neuron, Oligodendrocyte) ~ 
                       APOE + Gender + Age + PlasmaAB142P + PlasmaAB140P + PlasmaABRatio + PlasmapTau217 +          
                       PlasmapTau217_AB42Ratio + PlasmaGFAP +PlasmaNfL, 
                     data=pheno2, importance=TRUE)
  fit_no_covars <- rfsrc(Multivar(Microglia, Astrocyte, Neuron, Oligodendrocyte) ~ 
                       PlasmaAB142P + PlasmaAB140P + PlasmaABRatio + PlasmapTau217 +          
                       PlasmapTau217_AB42Ratio + PlasmaGFAP +PlasmaNfL, 
                     data=pheno2, importance=TRUE)
  
## Plot individual cell type outputs 

pdf_path <- file.path(output, "rfsrc_celltype_plots_no_covars.pdf")
pdf(file = pdf_path, width = 10, height = 6)
  plot(fit_no_covars, m.target = "Astrocyte")
  plot(fit_no_covars, m.target = "Microglia")
  plot(fit_no_covars, m.target = "Neuron")
  plot(fit_no_covars, m.target = "Oligodendrocyte")
 
  dev.off()
  
  
