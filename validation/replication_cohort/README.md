# Replication cohort pipeline (Figure 6: WSS module prediction)

This is a self-contained script set for an independent collaborator
cohort to replicate Figure 6: does a WSS module score predicted from 8
core plasma biomarkers (+ Age/Gender) reproduce the same
cognitive/diagnosis associations as the module score measured directly
from the full proteomic panel?

## What you need to provide

1. **Proteomics matrix** (CSV) — one row per sample, first column
   `SampleID`. By default, remaining columns are expected to be
   **SomaLogic probe IDs**, either probe names (e.g. `CRYBB2.10000.28.3`)
   or SeqIds (e.g. `seq.10000.28`), matching the bundled
   `soma_probe_to_symbol_key.csv`. Supply protein levels on their usual
   (normalized RFU) scale, not log-transformed or scaled; `01-WSSmod_prep.Rmd`
   rank-normalizes them before scoring. If your proteomics is from a
   different, already-symbol-named platform, set
   `map_probes_to_symbols: FALSE` in `01-WSSmod_prep.Rmd`'s params.
2. **Biomarker data** (CSV) — one row per sample: `SampleID`, `Age`,
   `Gender` (1 = male, 0 = not male), and the 8 raw core plasma
   biomarkers (`PlasmaPTau181`, `PlasmaAB142P`, `PlasmaAB140P`,
   `PlasmaABRatio`, `PlasmapTau217`, `PlasmapTau217_AB42Ratio`,
   `PlasmaGFAP`, `PlasmaNfL`). Note which of these were measured on a
   different platform than the reference (Lumipulse for pTau217, GFAP and
   NfL); see `self_normalize_biomarkers` in `01-WSSmod_prep.Rmd`.
3. **Clinical data** (CSV) — one row per sample: `SampleID`, `Age`,
   `Gender`, `CognitiveScore` (MoCA or equivalent), `Diagnosis_bin`
   (0 = control, 1 = case), plus any other demographic/technical
   variables you recorded — see "Choosing covariates" below.

## Bundled reference files (already in this folder)

These are platform/curation constants, not something you generate:

- `2024-08-30_Somalogic_reference_Key.txt` — SomaLogic probe → gene
  annotation key (`Key_2` matches your proteomics matrix's probe IDs).
- `prot_names.txt` — the full SomaScan panel's probe IDs.
- `soma_probe_to_symbol_key.csv` — a derived `(ix, seq, Symbol)` lookup
  used by `01-WSSmod_prep.Rmd` to build a gene-symbol matrix for
  `WSSmod::calculate_WSS()` from your proteomics, whether it is named by
  probe (`ix`, e.g. `CRYBB2.10000.28.3`) or by SeqId (`seq`, e.g.
  `seq.10000.28`). Built by `00-make_probe_key.R` from the probe list
  used to construct the WSS modules, so the module scores reproduce the
  original ones: one probe per gene symbol; multi-gene "complex" probes
  (e.g. `IL12B.IL23A.10365.132.3`) are duplicated across each constituent
  gene's column rather than averaged.
- `WSS_names.txt` / `WSS_names_pred.txt` — the subset of 45 WSS 
  modules preserved in the GNPC (actual / `.pred` predicted) 
  tested in the "clin predictions" section below.

Your cohort is assumed to be **cross-sectional** (one row/visit per
subject). This pipeline does not model repeated measures.

## Setup

```r
# install.packages("pak")
pak::pak("Longo-Lab/WSSmod")
```

You'll also need `data.table`, `optparse`, `tools`, `stringr`, `tictoc`,
`plyr`, `biomaRt`, `ggplot2`, `ggpubr`, `ggrepel`, `gprofiler2`,
`patchwork`, `qqman`, `Cairo`, `RNOmni`, `dplyr`, `RColorBrewer`,
`matrixStats`, `sva`, `pls`, `caret`, `ranger`, `vip`.

## Run order

1. **`01-WSSmod_prep.Rmd`** — edit the `params` block for your file
   paths, then knit/run. Computes actual WSS scores
   (`WSSmod::calculate_WSS`, after mapping your probe-level proteomics to
   gene symbols and rank-normalizing each protein; checked on the internal
   follow-up cohort, where it reproduces the original module scores at
   median r = 0.99) and predicted WSS scores
   (`WSSmod::normalize_wss_biomarkers` + `WSSmod::predict_WSS`), merges
   everything, and writes the **raw, un-QC'd** intermediate file
   `dataProt_replication.csv.gz`. Biomarkers are projected onto the bundled
   reference distribution, except those listed in `self_normalize_biomarkers`
   (default: pTau217, GFAP and the pTau217/AB42 ratio), which are
   rank-normalized within your cohort because Simoa/Quanterix values differ
   in scale from the Lumipulse reference (about 2.5x for pTau217 and GFAP;
   NfL is on the same scale, so it is projected). Check the printed mean/SD
   table: a projected biomarker with a mean near +/-1 or more has a platform
   offset and belongs in `self_normalize_biomarkers`. On the internal cohort
   with Simoa values substituted, this removed the bias from projecting
   (predicted WSS bias +0.057 to -0.002).
2. **`02-sample_QC_prediction.R`** — run this **first** for each module/clin
   you plan to test (each of the 12 modules for the "prot
   predictions" section; `CognitiveScore` and `Diagnosis_bin`, actual and
   predicted, for the "clin predictions" section). It evaluates candidate
   covariates via GLM/PLS/Random Forest feature importance
   (`*.feature_importance.pdf`) and, optionally, computes surrogate
   variables (specified with `--svs 6`, zero by default). 
   **Inspect its output before moving on** — see "Choosing covariates" 
   and "Surrogate variables" below.
3. **`00-commands.sh`** — edit the `$COVS` / `covs=(...)` placeholders at
   the top of each section to reflect what you decided in step 2, then
   `bash 00-commands.sh`. Runs the association analyses:
   - Section A: each WSS module (actual + predicted) vs. all proteins
     (`03-linreg_glm.R` → `04-Protein_QC.R`, which maps probe IDs back to
     genes via `-r 2024-08-30_Somalogic_reference_Key.txt`)
   - Section B: Cognitive Score and Diagnosis vs. each WSS module
     (actual + predicted) — the `covs` array is linked to
     `WSS_names.txt`/`WSS_names_pred.txt`, since actual vs. predicted
     scores may warrant different covariates/SVs.

## Choosing covariates

`Age,Gender` in the examples above are just a starting point, **not a
fixed recommendation** — different cohorts record different things. For
instance, a different real validation cohort (GNPC) included
`CountryOfOrigin` as a covariate that mattered there; your cohort may
have its own relevant demographic or technical variables (site, batch,
assay date, fasting status, etc.), or none at all recorded.

Use `02-sample_QC_prediction.R`'s feature-importance output
(`*.feature_importance.pdf`, showing GLM cross-validation, PLS, and
Random Forest impurity/permutation importance) to see which of your
candidate `-v` covariates are actually worth controlling for, rather
than assuming `Age,Gender` is sufficient or necessary.

## Surrogate variables (optional)

Surrogate variables (SVs, via `sva::sva()`) are **optional** —
`02-sample_QC_prediction.R`'s `--svs` defaults to `0` (skip). They're
most useful when your cohort has few or no recorded technical
covariates and you still want to account for unmeasured batch-like
structure.

If you use them, they are typically added iteratively rather than all at once:

1. Calulate `--svs 6` for a good sized cohort (500-1000). Can have sva estimate
   for a different size, but it usually grossly overpredicts number of svs
2. If specifying `--clin` for module to clinical phenotype comparisions, the 
   script will also generate a variance inflation factor (VIF) matrix to examine
   collinearity between a surrogate variable and each module being tested. If 
   these values are > 10, especially for multiple modules in a column, drop that 
   sv. Balance this with feature-importance scores (i.e. if an sv is important 
   and has borderline collinearity with a single module, you can keep it).
3. Fit 1 or 2 in `03-linreg_glm.R` and check resulting association `beta`/`se` 
   distributions (e.g. via `04-Protein_QC.R`'s `*_distributions.pdf`, or a quick
   histogram of your own) — inflated betas or implausibly small SEs are
   a sign of overfitting.
4. Add more SVs one at a time and refit. If the distributions spread beyond 
   beta +/- 10ish, or strong association pvalues flip wildly, back off.
5. For very large cohorts (i.e. GNPC), svs tend to have very small eigenvalues, 
   and contribute negligibly to a large balanced cohort.

## What's different from the internal "follow-up cohort" pipeline

This folder deliberately drops internal-cohort-only machinery that
doesn't apply to an external, cross-sectional collaborator dataset:

- No genomic PC merge, and no automatic repeated-visit projection of
  surrogate variables (`project_svs()`/train-vs-all split) — surrogate
  variables here are computed directly on the full (single-visit)
  sample if you choose to use them.
- Plain `glm` (`03-linreg_glm.R`) instead of the mixed-model
  `02-linreg_lme4qtl.R`, since a cross-sectional cohort has no repeated
  measures to model.
- No `06-TREAT-AD_GSEA.R` BioDomain enrichment step.
- No fixed internal batch covariates (`Storage_days`, `Subarray`,
  `Fed.Fasted`, `HybControlNormScale`) — covariates are chosen per-cohort.
