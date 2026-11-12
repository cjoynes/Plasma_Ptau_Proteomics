#!/usr/bin/env bash

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

# plotting dists
Rscript 01-Prediction_power.R  \
  -i "Merged.M13,Merged.M8,Merged.M27,Merged.M5,Merged.M2,Merged.M29,Merged.M1,Merged.M32,Merged.M15,Merged.M18,Merged.M21,Merged.M14,Merged.M4" \
  -c "../../prediction_models/STRING_WSS_v2/model_info/joinet/joinet_model_correlations_STRING_WSS_INT.csv"

# biodomain enrichments
Rscript 02-Combined-TREAT-AD_GSEA.R  \
  -i "Merged.M13,Merged.M8,Merged.M27,Merged.M5,Merged.M2,Merged.M29,Merged.M1,Merged.M32,Merged.M15,Merged.M18,Merged.M21,Merged.M14,Merged.M4" \
  --inpath "../../WSS_predict" \
  -n "dataProt_SS-205063.follow-up" \
  -s "top10.pred"




# External cohort
# ln -s "/mnt/box/Robert Butler's Externally Shareable Files/plasma_proteomics_replication_cohort" ~/rep
Rscript 02-Combined-TREAT-AD_GSEA.R  \
  -i "Merged.M13,Merged.M8,Merged.M27,Merged.M5,Merged.M2,Merged.M29,Merged.M1,Merged.M32,Merged.M15,Merged.M18,Merged.M21,Merged.M14,Merged.M4" \
  --inpath "~/rep/WSS_predict" \
  -n "dataProt_replication" \
  -s "top13.pred"
