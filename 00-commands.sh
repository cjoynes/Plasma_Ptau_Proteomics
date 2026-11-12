#!/usr/bin/env bash

###################################################################
#Script Name  : 00-commands.sh 
#Description  : Main association analysis wrapper script
#Usage        : bash 00-commands.sh
#Author       : Robert R Butler III
#Date Created : 2024-11-13
#Email        : rrbutler@stanford.edu
#Copyright (c) 2024
###################################################################

# Notes       :
# must have run scripts in the files_to_join and genomic folders first
# files_to_join - cleaning and merging covs, prots, biomarkers together
# genomic - using genetic data make PCs and GRM
# 
# must update values below to match your covariates After running 01-sample_QC.R,
# assess covariate importance for inclusion
#
# requires GNU parallel or further modification for batch cluster submission
# run with R v4.4.1


###################### Baseline cohort no genetics ############################
# biomarker names to be tested
names=(
  "PlasmaPTau181"
  "PlasmaAB140P" 
  "PlasmaAB142P" 
  "PlasmaNfL" 
  "PlasmaGFAP" 
  "PlasmapTau217"
  "PlasmaABRatio" 
  "PlasmapTau217_AB42Ratio" 
)
# enrichment modules
modules=("TMT-AD" "TMT-Resilience" "TMT-PTXtrial" "TMTSoma-CSF" "Wan_etal" "Soma-plasma")
# apoe or not
apoe=("" ".APOE")


# baseline ---------------------------------------------------------------------
# covariates to be included, biomarker ratios do not include HybControlNormScale
covs=(
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage'
)

# 01-sample_QC.R  ------------------------------------------
# Rscript 01-sample_QC.R \
  # -i "dataProt_SS-205063.merged.QC.csv" \
  # -o "baseline_PCs_relatedness/PlasmapTau217" \
  # -p "prot_names.txt" \
  # -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'   \
  # -c "PlasmapTau217" \
  # -n "dataProt_SS-205063.PC.PlasmapTau217" \
  # -m "genomic/mypcair.PCs.csv"

# baseline
parallel -j 8 Rscript 01-sample_QC.R \
  -i "dataProt_SS-205063.merged.QC.csv" \
  -o "baseline/{1}{3}"   \
  -p "prot_names.txt" \
  -v {2}{=3 s/\./\,/ =}  \
  -c {1} \
  -n "dataProt_SS-205063.{1}{3}" \
  ::: "${names[@]}" \
  :::+ "${covs[@]}"  \
  ::: "${apoe[@]}"

# covariates to be included, biomarker ratios do not include HybControlNormScale
covs=(
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6'
)

# PlateIds for SomaLogic run as covariates in cases where it showed up in the top 
#   half of feature importance for 01-sample_QC.R 
recodes=(
  'PlateId,P0028452,P0028439,P0028437,P0028449,P0028465,P0028565,P0028466'
  'PlateId,P0028451,P0028447,P0028439,P0028463,P0028462,P0028438,P0028430,P0028437'
  'PlateId,P0028462,P0028449,P0028439,P0028452,P0028565,P0028451,P0028430,P0028447,P0028437'
  'PlateId,P0028430,P0028437,P0028438'
  'PlateId,P0028449,P0028565,P0028439,P0028462,P0028452,P0028428,P0028463'
  'PlateId,P0028451,P0028438,P0028439,P0028449,P0028565,P0028452,P0028430'
  'PlateId,P0028448,P0028438,P0028437,P0028451,P0028565,P0028449,P0028447,P0028428,P0028452'
  'PlateId,P0028451,P0028449,P0028565,P0028439,P0028452'
)

# 02-linreg_glm.R ------------------------------------------
parallel -j 8 Rscript 02-linreg_glm.R \
  -i "baseline/{1}{4}/dataProt_SS-205063.{1}{4}.INTnorm.csv.gz"   \
  -o "baseline/{1}{4}" \
  -p "prot_names.txt" \
  -v {3}{=4 s/\./\,/ =}  \
  -d {2}  \
  -c "{1}"   \
  -n "dataProt_SS-205063.{1}{4}" \
  ::: "${names[@]}" \
  :::+ "${recodes[@]}"  \
  :::+ "${covs[@]}"  \
  ::: "${apoe[@]}"

# 03-Protein_QC.R ------------------------------------------
# Rscript 03-Protein_QC.R \
#   -i "baseline/PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.linreg_glm.csv.gz"  \
#   -o "baseline/PlasmaAB140P" \
#   -r "2024-08-30_Somalogic_reference_Key.txt" \
#   -n "dataProt_SS-205063.PlasmaAB140P"

parallel -j 1 Rscript 03-Protein_QC.R \
  -i "baseline/{1}{2}/dataProt_SS-205063.{1}{2}.linreg_glm.csv.gz" \
  -o "baseline/{1}{2}" \
  -r "2024-08-30_Somalogic_reference_Key.txt" \
  -n "dataProt_SS-205063.{1}{2}" \
  ::: "${names[@]}" \
  ::: "${apoe[@]}"



###################### Baseline cohort with PCs ############################
# same covs to generate svs and add in PCs  
# covariates to be included, biomarker ratios do not include HybControlNormScale
covs=(
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage'
)

# 01-sample_QC.R  ------------------------------------------
# Rscript 01-sample_QC.R \
  # -i "dataProt_SS-205063.merged.QC.csv" \
  # -o "baseline_PCs/Plasma_pTau217.APOE" \
  # -p "prot_names.txt" \
  # -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,APOE'   \
  # -c "Plasma_pTau217" \
  # -n "dataProt_SS-205063.PC.Plasma_pTau217.APOE" \
  # -m "genomic/mypcair.PCs.csv"

# baseline_PCs and baseline_PCs_relatedness use the same covariates
parallel -j 8 Rscript 01-sample_QC.R \
  -i "dataProt_SS-205063.merged.QC.csv" \
  -o "baseline_PCs/{1}{3}"   \
  -p "prot_names.txt" \
  -v {2}{=3 s/\./\,/ =}  \
  -c {1} \
  -n "dataProt_SS-205063.PC.{1}{3}" \
  -m "genomic/mypcair.PCs.csv" \
  ::: "${names[@]}" \
  :::+ "${covs[@]}"  \
  ::: "${apoe[@]}"

# 02-linreg_glm.R ------------------------------------------
# covariates to be included, biomarker ratios do not include HybControlNormScale
covs=(
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
)
# PlateIds for SomaLogic run as covariates in cases where it showed up in the top 
#   half of feature importance for 01-sample_QC.R 
recodes=(
  'PlateId,P0028462,P0028452,P0028449,P0028428,P0028565,P0028463,P0028437'
  'PlateId,P0028451,P0028463,P0028462,P0028438,P0028430,P0028437'
  'PlateId,P0028438,P0028439,P0028452,P0028449,P0028565,P0028451,P0028447,P0028430,P0028437'
  'PlateId,P0028430'
  'PlateId,P0028451,P0028449,P0028565,P0028452,P0028428,P0028439'
  'PlateId,P0028451,P0028438,P0028449,P0028565,P0028430,P0028439,P0028452'
  'PlateId,P0028448,P0028451,P0028438,P0028565,P0028449,P0028452,P0028428,P0028447'
  'PlateId,P0028451,P0028449,P0028565,P0028452'
)

# baseline_PCs still uses glm
parallel -j 8 Rscript 02-linreg_glm.R \
  -i "baseline_PCs/{1}{4}/dataProt_SS-205063.PC.{1}{4}.INTnorm.csv.gz"   \
  -o "baseline_PCs/{1}{4}" \
  -p "prot_names.txt" \
  -v {3}{=4 s/\./\,/ =}  \
  -d {2}  \
  -c "{1}"   \
  -n "dataProt_SS-205063.PC.{1}{4}" \
  ::: "${names[@]}" \
  :::+ "${recodes[@]}"  \
  :::+ "${covs[@]}"  \
  ::: "${apoe[@]}"

# 03-Protein_QC.R ------------------------------------------
# Rscript 03-Protein_QC.R \
#   -i "baseline/PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.linreg_glm.csv.gz"  \
#   -o "baseline/PlasmaAB140P" \
#   -r "2024-08-30_Somalogic_reference_Key.txt" \
#   -n "dataProt_SS-205063.PlasmaAB140P"

parallel -j 1 Rscript 03-Protein_QC.R \
  -i "baseline_PCs/{1}{2}/dataProt_SS-205063.PC.{1}{2}.linreg_glm.csv.gz" \
  -o "baseline_PCs/{1}{2}" \
  -r "2024-08-30_Somalogic_reference_Key.txt" \
  -n "dataProt_SS-205063.PC.{1}{2}" \
  ::: "${names[@]}" \
  ::: "${apoe[@]}"



###################### Baseline cohort PCs + GRM ############################
# same covs, svs and PCs, but we are adding a relatedness matrix and running lme4qtl

# 02-linreg_lme4qtl.R ------------------------------------------
# Rscript 02-linreg_lme4qtl.R \
#   -i "baseline_PCs_relatedness/Plasma_pTau217.APOE/dataProt_SS-205063.PC.Plasma_pTau217.APOE.INTnorm.csv.gz"   \
#   -o "baseline_PCs_relatedness/Plasma_pTau217.APOE" \
#   -p "prot_names.txt" \
#   -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,PlateId,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,APOE'   \
#   -d 'PlateId,P0028428,P0028439,P0028462,P0028437,P0028430,P0028448,P0028452,P0028447'  \
#   -c "Plasma_pTau217" \
#   -n "dataProt_SS-205063.PC.related.Plasma_pTau217.APOE" \
#   -r "genomic/mypcrelate.rds" \

# in our experience, increasing the threads at the script level is a negligible gain versus running single-threaded in parallel
parallel -j 8 Rscript 02-linreg_lme4qtl.R \
  -i "baseline_PCs/{1}{4}/dataProt_SS-205063.PC.{1}{4}.INTnorm.csv.gz" \
  -o "baseline_PCs_relatedness/{1}{4}"   \
  -p "prot_names.txt" \
  -v {3}{=4 s/\./\,/ =}  \
  -d {2}  \
  -t 1 \
  -c {1} \
  -n "dataProt_SS-205063.PC.related.{1}{4}" \
  -r "genomic/mypcrelate.rds" \
  ::: "${names[@]}" \
  :::+ "${recodes[@]}"  \
  :::+ "${covs[@]}"  \
  ::: "${apoe[@]}"


# 03-Protein_QC.R ------------------------------------------
# Rscript 03-Protein_QC.R \
#   -i "baseline_PCs_relatedness/PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.linreg_glm.csv.gz"  \
#   -o "baseline_PCs_relatedness/PlasmaAB140P" \
#   -r "2024-08-30_Somalogic_reference_Key.txt" \
#   -n "dataProt_SS-205063.PlasmaAB140P"

parallel -j 8 Rscript 03-Protein_QC.R \
  -i "baseline_PCs_relatedness/{1}{2}/dataProt_SS-205063.PC.related.{1}{2}.lme4qtl.csv.gz" \
  -o "baseline_PCs_relatedness/{1}{2}" \
  -r "2024-08-30_Somalogic_reference_Key.txt" \
  -n "dataProt_SS-205063.PC.related.{1}{2}" \
  ::: "${names[@]}" \
  ::: "${apoe[@]}"


# Looks like relatedness works well!

###################### Pathways, Celltypes, Enrichments ########################
# 04-TREAT-AD_correlation.R --------------------------------
# Rscript 04-TREAT-AD_correlation.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PlasmaAB140P" \
#   -c "PlasmaAB140P"

parallel -j 8 Rscript 04-TREAT-AD_correlation.R   \
  -i "baseline_PCs_relatedness/{1}{3}/dataProt_SS-205063.PC.related.{1}{3}.csv"  \
  -o "baseline_PCs_relatedness/{1}{3}" \
  -n "dataProt_SS-205063.PC.related.{1}{3}" \
  -c "{1}"  \
  -m "{2}" \
  -r "biodomains" \
  ::: "${names[@]}" \
  ::: "TREAT-AD" "TMTSoma-CSF" \
  ::: "${apoe[@]}"


# 05-Combined-Enrich.R -------------------------------------
# Rscript 05-Combined-Enrich.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PlasmaAB140P" \
#   -m "TMT-AD"

parallel -j 8 Rscript 05-Combined-Enrich.R  \
  -i "baseline_PCs_relatedness/{1}{3}/dataProt_SS-205063.PC.related.{1}{3}.csv"  \
  -o "baseline_PCs_relatedness/{1}{3}" \
  -n "dataProt_SS-205063.PC.related.{1}{3}" \
  -m "{2}" \
  -r "biodomains" \
  ::: "${names[@]}" \
  ::: "${modules[@]}" \
  ::: "${apoe[@]}"

# 06-TREAT-AD_GSEA.R ---------------------------------------
# be sure to have a functioning synapse client
# Rscript 06-TREAT-AD_GSEA.R \
  # -i "baseline_PCs_relatedness/PlasmaAB142P/dataProt_SS-205063.PC.related.PlasmaAB142P.csv"  \
  # -o "baseline_PCs_relatedness/PlasmaAB142P" \
  # -n "dataProt_SS-205063.PC.related.PlasmaAB142P" 

# conda activate synpase
parallel -j 8 Rscript 06-TREAT-AD_GSEA.R   \
  -i "baseline_PCs_relatedness/{1}{2}/dataProt_SS-205063.PC.related.{1}{2}.csv"  \
  -o "baseline_PCs_relatedness/{1}{2}" \
  -n "dataProt_SS-205063.PC.related.{1}{2}"  \
  ::: "${names[@]}" \
  ::: "${apoe[@]}"


# 07-Sharma_etal_correlation.R -----------------------------
# Rscript 07-Sharma_etal_correlation.R \
#   -i "PlasmaAB140P/dataProt_SS-205063.PlasmaAB140P.csv"  \
#   -o "PlasmaAB140P" \
#   -n "dataProt_SS-205063.PlasmaAB140P" \
#   -c "PlasmaAB140P" 

parallel -j 8 Rscript 07-Sharma_etal_correlation.R   \
  -i "baseline_PCs_relatedness/{1}{2}/dataProt_SS-205063.PC.related.{1}{2}.csv"  \
  -o "baseline_PCs_relatedness/{1}{2}" \
  -n "dataProt_SS-205063.PC.related.{1}{2}" \
  -c "{1}" \
  -r "cell_type_prot" \
  ::: "${names[@]}" \
  ::: "${apoe[@]}"


###################### STRING Modules ######################################

# 08-String_module_construction.R ----------------------
Rscript 08-String_module_construction.R \
  -i "PlasmaGFAP,PlasmaNfL,PlasmaPTau181,PlasmapTau217,PlasmaAB140P,PlasmaAB142P,PlasmaABRatio,PlasmapTau217_AB42Ratio" \
  -n "dataProt_SS-205063.STRING" \
  -o "module_scores"

Rscript 08a-String_module_testing.R


# 09-WSS_construction.R ----------------------------
Rscript 09-WSS_construction.R \
  -i "module_scores/STRING_merged_scaled_scores.csv" \
  -n "dataProt_SS-205063.STRING" \
  --INT


# 10 - Prediction Model construction ---------------------
Rscript 10-all_prediction_models.R

# 11-Eigengene_clustering.Rmd ----------------------------
# Run interactively



###################### WSS modules clin analysis ###############################
# 12-sample_QC_WSS_clin.R ------------------------------
# Now just run WSS x clin analyses on well conserved WSS modules

# Binary variables -----------------------
names=(
  "HYPERCHO"
  "HYPERTEN"
  "DEP2YRS"
  "ANXIETY"
  "DIABETES"
  "Diagnosis_bin"
  "last_diag_bin"
)

# QC and surrogates
parallel -j 7 Rscript 12-sample_QC_WSS_clin.R \
  -i "dataProt_SS-205063.baseline_clin.csv" \
  -o "WSS_clin/{}"   \
  -v "Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,Storage,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8" \
  -c {} \
  -n "dataProt_SS-205063.baseline_clin_binary.{}" \
  --binomial \
  ::: "${names[@]}" 

# 02 - lme4qtl ------
# covariates to be included, biomarker ratios do not include HybControlNormScale
covs=(
  'Age,Gender,Storage_days,Subarray,Fed.Fasted,HybControlNormScale,Storage,sv3,sv4,sv5,PC1,PC7'
  'Age,Gender,Subarray,Fed.Fasted,Storage,sv4,sv5,PC1,PC3,PC4,PC6,PC8'
  'Age,Gender,Storage_days,Subarray,Fed.Fasted,HybControlNormScale,Storage,sv3,sv4,PC6,PC7'
  'Age,Gender,Subarray,Fed.Fasted,HybControlNormScale,Storage,sv3,sv6,PC1,PC2,PC4,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Subarray,Fed.Fasted,Storage,sv3,sv5,sv6,PC1,PC7,PC8'
  'Age,Gender,Storage_days,Subarray,Storage,sv2,sv3,sv4,sv6,PC5,PC6,PC7,PC8'
  'Age,Gender,Storage_days,Subarray,HybControlNormScale,Storage,sv3,PC1,PC5,PC6,PC7,PC8'
)
# Subarray for SomaLogic run as covariates in cases where it showed up in the top 
#   half of feature importance for 01-sample_QC.R 
recodes=(
  'Subarray,B,D'
  'Subarray,B,E|Storage,Negative.80C'
  'Subarray,B,D,H'
  'Subarray,D,E,H|Storage,Negative.20C.LowVolume'
  'Storage,Negative.80C'
  'Subarray,D,F|Storage,Negative.80C'
  'Subarray,D,F,H|Storage,Negative.80C'
)

parallel -j 7 Rscript 02-linreg_lme4qtl.R \
  -i "WSS_clin/{1}/dataProt_SS-205063.baseline_clin_binary.{1}.WSSnorm.csv.gz" \
  -o "WSS_clin/{1}"   \
  --prots "WSS_names.txt" \
  -v {3} \
  -d {2} \
  -s 'Age,Storage_days,Fed.Fasted' \
  -c "{1}" \
  -t 2 \
  -n "dataProt_SS-205063.baseline_clin_binary.{1}" \
  --binomial \
  -f "0,1" \
  ::: "${names[@]}" \
  :::+ "${recodes[@]}"  \
  :::+ "${covs[@]}"  


# Continuous variables -----------------------
names=(
  "MoCA"
  "CDRSUM"
  "last_MoCA"
  # "last_CDRSUM"
  # "MoCA_velocity"
  # "MoCA_v_floored"
)

# QC and surrogates
parallel -j 4 Rscript 12-sample_QC_WSS_clin.R \
  -i "dataProt_SS-205063.baseline_clin.csv" \
  -o "WSS_clin/{}"   \
  -v "Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,Storage,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8" \
  -c {} \
  -n "dataProt_SS-205063.baseline_clin_continuous.{}" \
  ::: "${names[@]}" 

# 02 - lme4qtl ------
# covariates to be included, biomarker ratios do not include HybControlNormScale
covs=(
  'Age,Gender,Storage_days,Subarray,sv2,sv3,sv4,sv6,PC1,PC2,PC5,PC6,PC7'
  'Age,Gender,Storage_days,Subarray,Fed.Fasted,Storage,sv3,sv4,PC1,PC4,PC5,PC8'
  'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,sv3,sv4,sv5,sv6,PC1,PC2,PC4,PC6,PC7,PC8'
  # 'Age,Gender,Storage_days,Subarray,Fed.Fasted,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  # 'Age,Gender,Storage_days,Subarray,Fed.Fasted,HybControlNormScale,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
  # 'Age,Gender,Storage_days,Subarray,Fed.Fasted,HybControlNormScale,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'
)
# Subarray for SomaLogic run as covariates in cases where it showed up in the top 
#   half of feature importance for 01-sample_QC.R 
recodes=(
  'Subarray,D'
  'Subarray,D,G|Storage,Negative.80C'
  'Subarray,D'
  # 'Subarray,D,G'
  # 'Subarray,B,C'
  # 'Subarray,B,C,H'
)

parallel -j 3 Rscript 02-linreg_lme4qtl.R \
  -i "WSS_clin/{1}/dataProt_SS-205063.baseline_clin_continuous.{1}.WSSnorm.csv.gz" \
  -o "WSS_clin/{1}"   \
  --prots "WSS_names.txt" \
  -v {3} \
  -d {2} \
  -s 'Age,Storage_days,Fed.Fasted' \
  -c "{1}" \
  -t 2 \
  -n "dataProt_SS-205063.baseline_clin_continuous.{1}" \
  ::: "${names[@]}" \
  :::+ "${recodes[@]}"  \
  :::+ "${covs[@]}"  


###################### WSS modules prot analysis ###############################
# should just run with the default association scripts
names=(
  "Merged.M1" 
  "Merged.M11" 
  "Merged.M7" 
  "Merged.M9" 
  "Merged.M2" 
  "PlasmaPTau181.M8" 
  "Merged.M6" 
  "PlasmaAB140P.M17" 
  "Merged.M12" 
  "Merged.M4" 
  "Merged.M5" 
  "PlasmaAB140P.M20" 
  "Merged.M8" 
  "PlasmaPTau181.M35" 
  "Merged.M3" 
  "Merged.M26" 
  "Merged.M13" 
  "Merged.M10" 
  "PlasmaPTau181.M18" 
  "Merged.M14" 
  "Merged.M19" 
  "Merged.M16" 
  "PlasmaPTau181.M26" 
  "Merged.M29" 
  "Merged.M24" 
  "PlasmaPTau181.M13" 
  "PlasmaPTau181.M47" 
  "Merged.M25" 
  "PlasmaAB140P.M35" 
  "PlasmaPTau181.M46" 
  "Merged.M20" 
  "Merged.M32" 
  "Merged.M22" 
  "Merged.M15" 
  "Merged.M21" 
  "Merged.M18" 
  "PlasmaPTau181.M5" 
  "PlasmaAB140P.M15" 
  "Merged.M27" 
  "PlasmaAB140P.M40" 
  "Merged.M17" 
  "PlasmaNfL.M3" 
  "PlasmaAB140P.M18" 
  "PlasmaPTau181.M30" 
  "PlasmaGFAP.M12"
)

# 01-sample_QC.R  ------------------------------------------
parallel -j 8 Rscript 01-sample_QC.R \
  -i "dataProt_SS-205063.merged.QC.csv" \
  -o "WSS_prot/{}"   \
  -p "prot_names.txt" \
  -e "module_scores/dataProt_SS-205063.STRING.INT.baseline.WSS.csv" \
  -v "Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,Storage"  \
  -c {} \
  -n "dataProt_SS-205063.baseline_WSS.{}" \
  --skip-qc \
  ::: "${names[@]}" 


# 02-linreg_lme4qtl.R ------------------------------------------
# Rscript 02-linreg_lme4qtl.R \
#   -i "WSS_prot/Merged.M1/dataProt_SS-205063.baseline_WSS.Merged.M1.INTnorm.csv.gz"   \
#   -o "WSS_prot/Merged.M1" \
#   -p "prot_names.txt" \
#   -v 'Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,APOE'   \
#   -d 'Subarray,E,F,H'  \
#   -c "Merged.M1" \
#   -n "dataProt_SS-205063.baseline_WSS.Merged.M1" \
#   -r "genomic/mypcrelate.rds" \

# refactored to make use of multi-threaded in future.lapply
parallel -j 1 Rscript 02-linreg_lme4qtl.R \
  -i "WSS_prot/{}/dataProt_SS-205063.baseline_WSS.{}.INTnorm.csv.gz" \
  -o "WSS_prot/{}"   \
  -v "Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8"  \
  -t 18 \
  -c {} \
  -n "dataProt_SS-205063.baseline_WSS.{}" \
  ::: "${names[@]}" 


# 03-Protein_QC.R ------------------------------------------
# Rscript 03-Protein_QC.R \
  # -i "WSS_prot/Merged.M3/dataProt_SS-205063.baseline_WSS.Merged.M3.lme4qtl.csv.gz"  \
  # -o "WSS_prot/Merged.M3" \
  # -r "2024-08-30_Somalogic_reference_Key.txt" \
  # -n "dataProt_SS-205063.baseline_WSS.Merged.M3"
  # -p "cohens_d"

parallel -j 1 Rscript 03-Protein_QC.R \
  -i "WSS_prot/{}/dataProt_SS-205063.baseline_WSS.{}.lme4qtl.csv.gz" \
  -o "WSS_prot/{}" \
  -r "2024-08-30_Somalogic_reference_Key.txt" \
  -n "dataProt_SS-205063.baseline_WSS.{}" \
  ::: "${names[@]}" 

# 06-TREAT-AD_GSEA.R ---------------------------------------
# be sure to have a functioning synapse client
# Rscript 06-TREAT-AD_GSEA.R \
  # -i "WSS_prot/Merged.M3/dataProt_SS-205063.baseline_WSS.Merged.M3.csv"  \
  # -o "WSS_prot/Merged.M3" \
  # -n "dataProt_SS-205063.baseline_WSS.Merged.M3" \
  # -p "cohens_d"

# conda activate synpase
parallel -j 1 Rscript 06-TREAT-AD_GSEA.R   \
  -i "WSS_prot/{}/dataProt_SS-205063.baseline_WSS.{}.csv"  \
  -o "WSS_prot/{}" \
  -n "dataProt_SS-205063.baseline_WSS.{}"  \
  ::: "${names[@]}" 





###################### Follow-up cohort prot predictions #######################
# 13-sample_QC_prediction.R --------------------------------
# module names to be tested (all R2 > 0.25)
names=(
  # "Merged.M13"
  # "Merged.M8"
  # "Merged.M27"
  # "Merged.M5"
  # "Merged.M2"
  # "Merged.M29"
  # "Merged.M1"
  # "Merged.M32"
  # "Merged.M15"
  # "Merged.M18"
  # "Merged.M21"
  # "Merged.M14"
  "Merged.M4"
)

# prediction or not
pred=("" ".pred")
# pred=("")


# Rscript 11-sample_QC_eigengene.R \
  # -i "dataProt_SS-205063.merged.QC.csv" \
  # -e "module_scores/dataProt_SS-205063.STRING.INT.eigengenes.csv" \
  # --preds "prediction_models/STRING_WSS/joinet/INT/Biomarkers_Age_Sex.csv" \
  # -o "WSS_predict/Merged.M13" \
  # -p "prot_names.txt" \
  # -v "Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,Storage"   \
  # -c "Merged.M13" \
  # -n "dataProt_SS-205063.follow-up.Merged.M13" \
  # -m "genomic/mypcair.PCs.csv"

# All models use the same covariates here
parallel -j 2 Rscript 13-sample_QC_prediction.R \
  -o "WSS_predict/{}"   \
  -v "Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,Storage" \
  -c {} \
  -n "dataProt_SS-205063.follow-up.{}" \
  ::: "${names[@]}" 


# 02-linreg_lme4qtl.R ------------------------------------------
# Rscript 02-linreg_lme4qtl.R \
  # -i "dichotomous_eigengenes/MEMerged-M11/dataProt_SS-205063.quartile.MEMerged-M11.INTnorm.csv.gz"   \
  # -o "dichotomous_eigengenes/MEMerged-M11" \
  # -p "prot_names.txt" \
  # -v 'Age,Gender,Storage_days,Study,Fed.Fasted,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'   \
  # -s 'Age,Storage_days,Fed.Fasted' \
  # -c "quartile" \
  # -n "dataProt_SS-205063.quartile.MEMerged-M11" \
  # -r "genomic/mypcrelate.rds" \
  # --binomial

# multi-threaded, with future.lapply solving 
parallel -j 1 Rscript 02-linreg_lme4qtl.R \
  -i "WSS_predict/{1}/dataProt_SS-205063.follow-up.{1}.INTnorm.csv.gz" \
  -o "WSS_predict/{1}{2}"   \
  -v 'Age,Gender,Storage_days,Study,Fed.Fasted,HybControlNormScale,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8'   \
  -s 'Age,Storage_days,Fed.Fasted' \
  -c "{1}{2}" \
  -t 12 \
  -n "dataProt_SS-205063.follow-up.{1}{2}" \
  ::: "${names[@]}" \
  ::: "${pred[@]}" 


# 03-Protein_QC.R ------------------------------------------
# Rscript 03-Protein_QC.R \
  # -i "WSS_predict/Merged.M3.pred/dataProt_SS-205063.follow-up.Merged.M3.pred.lme4qtl.csv.gz"  \
  # -o "WSS_predict/Merged.M3.pred" \
  # -r "2024-08-30_Somalogic_reference_Key.txt" \
  # -n "dataProt_SS-205063.follow-up.Merged.M3.pred"
  # -p "cohens_d"



parallel -j 1 Rscript 03-Protein_QC.R \
  -i "WSS_predict/{1}{2}/dataProt_SS-205063.follow-up.{1}{2}.lme4qtl.csv.gz" \
  -o "WSS_predict/{1}{2}" \
  -r "2024-08-30_Somalogic_reference_Key.txt" \
  -n "dataProt_SS-205063.follow-up.{1}{2}" \
  ::: "${names[@]}" \
  ::: "${pred[@]}"

# 06-TREAT-AD_GSEA.R ---------------------------------------
# be sure to have a functioning synapse client
# Rscript 06-TREAT-AD_GSEA.R \
  # -i "WSS_predict/Merged.M3/dataProt_SS-205063.follow-up.Merged.M3.csv"  \
  # -o "WSS_predict/Merged.M3" \
  # -n "dataProt_SS-205063.follow-up.Merged.M3" \
  # -p "cohens_d"

# conda activate synpase
parallel -j 1 Rscript 06-TREAT-AD_GSEA.R   \
  -i "WSS_predict/{1}{2}/dataProt_SS-205063.follow-up.{1}{2}.csv"  \
  -o "WSS_predict/{1}{2}" \
  -n "dataProt_SS-205063.follow-up.{1}{2}"  \
  ::: "${names[@]}" \
  ::: "${pred[@]}"



###################### Follow-up cohort clin predictions #######################
# Continuous variables -----------------------
names=(
  "MoCA"
  # "last_MoCA"
)

# QC and surrogates
parallel -j 2 Rscript 13-sample_QC_prediction.R \
  -i "dataProt_SS-205063.follow-up_clin.csv" \
  -p "WSS_names.txt" \
  --clin \
  -o "WSS_clin_actual/{}" \
  -v "Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,Storage,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8" \
  -c {} \
  -n "dataProt_SS-205063.follow-up_clin_continuous.{}" \
  ::: "${names[@]}" 

# 02 - lme4qtl ------
# covariates to be included, biomarker ratios do not include HybControlNormScale
covs=(
  'Age,Gender,Storage_days,Subarray,HybControlNormScale,sv2,sv4,sv6,PC1,PC2,PC3,PC5,PC6'
  # 'Age,Gender,Storage_days,Subarray,HybControlNormScale,Storage,sv1,sv2,sv3,sv4,sv5,sv6,PC1,PC2,PC3,PC4,PC5,PC6,PC7'
)
# Subarray for SomaLogic run as covariates in cases where it showed up in the top 
#   half of feature importance for 01-sample_QC.R 
recodes=(
  'Subarray,B,E,H'
  # 'Subarray,B,E,H'
)

parallel -j 1 Rscript 02-linreg_lme4qtl.R \
  -i "WSS_clin_actual/{1}/dataProt_SS-205063.follow-up_clin_continuous.{1}.INTnorm.csv.gz" \
  -o "WSS_clin_actual/{1}"   \
  --prots "WSS_names.txt" \
  -v {3} \
  -d {2} \
  -s 'Age,Storage_days' \
  -c "{1}" \
  -t 18 \
  -n "dataProt_SS-205063.follow-up_clin_continuous.{1}" \
  ::: "${names[@]}" \
  :::+ "${recodes[@]}"  \
  :::+ "${covs[@]}"  


# QC and surrogates
parallel -j 4 Rscript 13-sample_QC_prediction.R \
  -i "dataProt_SS-205063.follow-up_clin.csv" \
  -p "WSS_names_pred.txt" \
  --clin \
  -o "WSS_clin_predict/{}" \
  -v "Age,Gender,Storage_days,Study,Subarray,Fed.Fasted,HybControlNormScale,Storage,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8" \
  -c {} \
  -n "dataProt_SS-205063.follow-up_clin_continuous.{}" \
  ::: "${names[@]}" 

# covariates to be included, biomarker ratios do not include HybControlNormScale
covs=(
  'Age,Gender,Storage_days,Subarray,HybControlNormScale,sv1,sv6,PC1,PC2,PC6,PC8'
  # 'Age,Gender,Fed.Fasted,Storage,sv1,sv2,sv3,sv4,PC1,PC2,PC5,PC6,PC8'
)
# Subarray for SomaLogic run as covariates in cases where it showed up in the top 
#   half of feature importance for 01-sample_QC.R 
recodes=(
  'Subarray,B,E'
  # ''
)
parallel -j 1 Rscript 02-linreg_lme4qtl.R \
  -i "WSS_clin_predict/{1}/dataProt_SS-205063.follow-up_clin_continuous.{1}.INTnorm.csv.gz" \
  -o "WSS_clin_predict/{1}"   \
  --prots "WSS_names_pred.txt" \
  -v {3} \
  -d {2} \
  -s 'Age,Storage_days' \
  -c "{1}" \
  -t 18 \
  -n "dataProt_SS-205063.follow-up_clin_continuous.{1}" \
  ::: "${names[@]}" \
  :::+ "${recodes[@]}"  \
  :::+ "${covs[@]}"  


# Binary variables -----------------------
names=(
  "Diagnosis_bin"
)

# QC and surrogates
parallel -j 1 Rscript 13-sample_QC_prediction.R \
  -i "dataProt_SS-205063.follow-up_clin.csv" \
  -p "WSS_names.txt" \
  --clin \
  -o "WSS_clin_actual/{}" \
  -v "Age,Gender,Storage_days,Subarray,Fed.Fasted,HybControlNormScale,Storage,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8" \
  -c {} \
  -n "dataProt_SS-205063.follow-up_clin_binary.{}" \
  --binomial \
  ::: "${names[@]}" 

# 02 - lme4qtl ------
# covariates to be included, biomarker ratios do not include HybControlNormScale
covs=(
  'Age,Gender,Storage_days,Subarray,Fed.Fasted,HybControlNormScale,sv2,sv3,sv6,PC1,PC3,PC5,PC6,PC7,PC8'
)
# Subarray for SomaLogic run as covariates in cases where it showed up in the top 
#   half of feature importance for 01-sample_QC.R 
recodes=(
  'Subarray,B'
)

parallel -j 1 Rscript 02-linreg_lme4qtl.R \
  -i "WSS_clin_actual/{1}/dataProt_SS-205063.follow-up_clin_binary.{1}.INTnorm.csv.gz" \
  -o "WSS_clin_actual/{1}"   \
  --prots "WSS_names.txt" \
  -v {3} \
  -d {2} \
  -s 'Age,Storage_days' \
  -c "{1}" \
  -t 18 \
  -n "dataProt_SS-205063.follow-up_clin_binary.{1}" \
  --binomial \
  -f "0,1" \
  ::: "${names[@]}" \
  :::+ "${recodes[@]}"  \
  :::+ "${covs[@]}"  


# QC and surrogates
parallel -j 1 Rscript 13-sample_QC_prediction.R \
  -i "dataProt_SS-205063.follow-up_clin.csv" \
  -p "WSS_names_pred.txt" \
  --clin \
  -o "WSS_clin_predict/{}" \
  -v "Age,Gender,Storage_days,Subarray,Fed.Fasted,HybControlNormScale,Storage,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8" \
  -c {} \
  -n "dataProt_SS-205063.follow-up_clin_binary.{}" \
  --binomial \
  ::: "${names[@]}" 

# covariates to be included, biomarker ratios do not include HybControlNormScale
covs=(
  'Age,Gender,Storage_days,Fed.Fasted,Storage,sv1,PC1,PC3,PC5'
)
# Subarray for SomaLogic run as covariates in cases where it showed up in the top 
#   half of feature importance for 01-sample_QC.R 
recodes=(
  'Storage,Negative.80C'
)
parallel -j 1 Rscript 02-linreg_lme4qtl.R \
  -i "WSS_clin_predict/{1}/dataProt_SS-205063.follow-up_clin_binary.{1}.INTnorm.csv.gz" \
  -o "WSS_clin_predict/{1}"   \
  --prots "WSS_names_pred.txt" \
  -v {3} \
  -d {2} \
  -s 'Age,Storage_days' \
  -c "{1}" \
  -t 18 \
  -n "dataProt_SS-205063.follow-up_clin_binary.{1}" \
  --binomial \
  -f "0,1" \
  --min-sep 1 \
  ::: "${names[@]}" \
  :::+ "${recodes[@]}"  \
  :::+ "${covs[@]}"  



###################### External cohort prot predictions #######################

# First steps 01, 02 run locally

# 13-sample_QC_prediction.R --------------------------------
# module names to be tested (all R2 > 0.25)
names=(
  "Merged.M13"
  "Merged.M8"
  "Merged.M27"
  "Merged.M5"
  "Merged.M2"
  "Merged.M29"
  "Merged.M1"
  "Merged.M32"
  "Merged.M15"
  "Merged.M18"
  "Merged.M21"
  "Merged.M14"
  "Merged.M4"
)

# prediction or not
pred=("" ".pred")
# pred=("")


# 03-Protein_QC.R ------------------------------------------
# Rscript 03-Protein_QC.R \
  # -i "WSS_predict/Merged.M3.pred/dataProt_SS-205063.follow-up.Merged.M3.pred.lme4qtl.csv.gz"  \
  # -o "WSS_predict/Merged.M3.pred" \
  # -r "2024-08-30_Somalogic_reference_Key.txt" \
  # -n "dataProt_SS-205063.follow-up.Merged.M3.pred"
  # -p "cohens_d"

ln -s "/mnt/box/Robert Butler's Externally Shareable Files/plasma_proteomics_replication_cohort" ~/rep
rep_dir="~/rep"


parallel -j 1 Rscript 03-Protein_QC.R \
  -i "${rep_dir}/WSS_predict/{1}{2}/dataProt_replication.{1}{2}.linreg_glm.csv.gz" \
  -o "${rep_dir}/WSS_predict/{1}{2}" \
  -r "2024-08-30_Somalogic_reference_Key.txt" \
  -n "dataProt_replication.{1}{2}" \
  ::: "${names[@]}" \
  ::: "${pred[@]}"

# 06-TREAT-AD_GSEA.R ---------------------------------------
# be sure to have a functioning synapse client
# Rscript 06-TREAT-AD_GSEA.R \
  # -i "WSS_predict/Merged.M3/dataProt_SS-205063.follow-up.Merged.M3.csv"  \
  # -o "WSS_predict/Merged.M3" \
  # -n "dataProt_SS-205063.follow-up.Merged.M3" \
  # -p "cohens_d"

# conda activate synpase
parallel -j 1 Rscript 06-TREAT-AD_GSEA.R   \
  -i "${rep_dir}/WSS_predict/{1}{2}/dataProt_replication.{1}{2}.csv"  \
  -o "${rep_dir}/WSS_predict/{1}{2}" \
  -n "dataProt_replication.{1}{2}"  \
  ::: "${names[@]}" \
  ::: "${pred[@]}"





