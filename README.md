# Predicting anti-TNF therapy failure in ulcerative colitis

R code for the study:

> Derakhshan Nazari MH, Shahrokh S, Ghanbari-Maman L, Maleknia S, Ghorbaninejad M,
> Meyfour A. **Prediction of anti-TNF therapy failure in ulcerative colitis patients
> by ensemble machine learning: A prospective study.** *Heliyon* 9 (2023) e21154.
> https://doi.org/10.1016/j.heliyon.2023.e21154

Roughly 40 % of ulcerative colitis (UC) patients do not respond to anti-TNF
therapy. Using baseline mucosal biopsies from 79 anti-TNF-naive patients across
four cohorts, this study characterises the transcriptional differences between
responders and non-responders and builds a classifier that flags likely
non-responders *before* treatment begins. Two microarray series were integrated
into a discovery cohort, differential expression and functional/network analyses
identified the altered genes and pathways, and repeated LASSO feature selection
followed by a two-layer stack ensemble classifier produced a five-mRNA panel
(IL13RA2, HCAR3, CSF3, INHBA, MMP1) that predicts response with ~95 % average
accuracy. The panel was validated on an independent cohort and experimentally by
RT-qPCR in a prospective adalimumab-treated patient cohort.

## Data

| Dataset | Platform | Role |
| --- | --- | --- |
| GSE12251, GSE16879 | Affymetrix U133 Plus 2.0 (GPL570) | Integrated to form the discovery cohort (35 UC samples) |
| GSE73661 | Affymetrix Human Gene 1.0 ST (GPL6244) | Independent test cohort |
| RT-qPCR cohort | — | 22 prospective patients treated with adalimumab (CinnoRa), Taleghani hospital |

Microarray data are public and available from GEO. The scripts expect the series
matrices and the GPL570 annotation file as local text files.

## Scripts

`00_full_pipeline_as_submitted.R` is the consolidated script that accompanied the
manuscript, covering preprocessing through the ensemble classifier in one file.
The numbered scripts are the per-stage working versions of the same analyses,
plus the validation steps the consolidated file does not include.

| Script | Description |
| --- | --- |
| `00_full_pipeline_as_submitted.R` | Consolidated pipeline submitted with the paper: preprocessing and batch correction, differential expression, correlation heatmap, pathway plot, Pearson correlation, network construction, LASSO feature selection, and the ensemble classifier. |
| `01_preprocessing_batch_correction.R` | Merges the two discovery series with the platform annotation, applies quantile normalisation, and removes batch effects with ComBat. |
| `02_differential_expression.R` | Differential expression between non-responders and responders with limma. |
| `03_correlation_heatmap.R` | Correlation heatmap of the differentially expressed genes. |
| `04_gene_correlation.R` | Pairwise Pearson correlation used to define the co-expression edges of the network. |
| `05_network_construction.R` | Builds the node and edge attribute tables for the co-expression/protein-interaction network analysed in Cytoscape. |
| `06_pathway_bubble_plot.R` | Bubble plot of the enriched pathways. |
| `07_lasso_feature_selection.R` | Repeated (50×) 10-fold cross-validated LASSO feature selection, ranking genes by average absolute coefficient. |
| `08_ensemble_machine_learning.R` | Two-layer stack ensemble (LDA, random forest, neural network, KNN, combined by a GLM meta-learner) with repeated cross-validation and ROC/AUC evaluation. |
| `09_independent_cohort_validation.R` | Applies the trained model to the independent test cohort (GSE73661). |
| `10_rt_qpcr_analysis.R` | Processing and analysis of the RT-qPCR measurements from the prospective cohort. |
| `11_rt_qpcr_machine_learning.R` | Trains and evaluates the classifier on the RT-qPCR panel data. |
| `12_roc_curve_analysis.R` | ROC curves and AUC for the final panel. |

## Requirements

R 4.0 or later with, among others: `limma`, `sva`, `glmnet`, `caret`,
`caretEnsemble`, `randomForest`, `e1071`, `MASS`, `pROC`, `igraph`, `pheatmap`,
`ggplot2`, `dplyr` and `data.table`. Each script loads what it needs at the top.

Several steps of the published analysis were carried out in external tools rather
than in R, and are therefore not reproduced here:

* **CIBERSORTx** (LM22 signature) for immune cell deconvolution
* **Cytoscape** 3.8.2 with **STRING**, **CentiScape 2.2** and **CytoHubba** for
  network construction and centrality analysis
* **Gene Ontology Resource** and **MSigDB** for enrichment analysis
* **SPSS** for the immune cell fraction statistics

## A note on file paths

Every path in these scripts has been replaced with a `path/to/…` placeholder.
Point them at your own directories before running; most scripts also begin with a
`setwd()` call that needs the same treatment.

## Citation

If you use this code, please cite the paper above.
