Dissertation "Assessment of modeling strategies for drug response prediction in cell lines and xenografts"
================

Roman Kurilov, 2018
--------
***


Summary
--------

Despite significant progress in cancer research, effective cancer treatment is still a challenge. Cancer treatment approaches are shifting from standard cytotoxic chemotherapy regimens towards a precision oncology paradigm, where a choice of treatment is personalized, i.e. based on a tumor’s molecular features. In order to match tumor molecular features with therapeutics we need to identify biomarkers of response and build predictive models. Recent growth of large-scale pharmacogenomics resources which combine drug sensitivity and multi-omics information on a large number of samples provides necessary data for biomarker identification and drug response modelling. However, although many efforts of using this information for drug response prediction have been made, our ability to accurately predict drug response using genetic data remains limited.

In this work we used pharmacogenomics data from the largest publicly available studies in order to systematically assess various aspects of the drug response model-building process with the ultimate goal of improving prediction accuracy. We applied several machine learning methods (regularized regression, support vector machines, random forest) for predicting response to a number of drugs. We found that while accuracy of response prediction varies across drugs (in most of the cases R2 values vary between 0.1 and 0.3), different machine learning algorithms applied for the the same drug have similar prediction performance. Experiments with a range of different training sets for the same drug showed that predictive power of a model depends on the type of molecular data, the selected drug response metric, and the size of the training set. It depends less on number of features selected for modelling and on class imbalance in training set. We also implemented and tested two methods for improving consistency for pharmacogenomics data coming from different datasets.

We tested our ability to correctly predict response in xenografts and patients using models trained on cell lines. Only in a fraction of the tested cases we managed to get reasonably accurate predictions, particularly in case of response to erlotinib in the NSCLC xenograft cohort, and in cases of responses to erlotinib and docetaxel in the NSCLC and BRCA patient cohorts respectively.

This work also includes two applied pharmacogenomics analyses. The first is an analysis of a drug-sensitivity screen performed on a panel of Burkitt cell lines. This combines unsupervised data exploration with supervised modelling. The second is an analysis of drug-sensitivity data for the DKFZ-608 compound and the generation of the corresponding response prediction model.

In summary, we applied machine learning techniques to available high-throughput pharmacogenomics data to study the determinants of accurate drug response prediction. Our results can help to draft guidelines for building accurate models for personalized drug response prediction and therefore contribute to advancing of precision oncology.



Reproducibility of the Analyses
--------------------------------------------
In all scripts, the user needs to change the variable **path** to the repository root folder.  
Prior to running scripts from the chapters **data.R** script from the root folder should be executed in order produce the necessary data files for the Chapters 2,3,7.

Code for Chapters 4 and 5 is available via this repository: https://github.com/RomaHD/DrugRespPrediction/ 

Chapter 2
-------------------------------

|     | name | thesis's section                            |
|-----|------|------------------------------------------|
|1. | drug_resp_consistency.R  | "Drug response consistency" |
|2. | biomarkers_consistency.R  | "Biomarkers' consistency" |

Chapter 3
-------------------------------

|     | name | thesis's section                             |
|-----|------|------------------------------------------|
|1. | multi_task.R  | "Multi-task glmnet models" |
|2. | aggregated.R  | "Modelling on aggregated data" |
|3. | feature_interactions.R  | "Modelling with feature interactions" |
|4. | weights.R  | "Modelling with weights" |
|5. | imbalance_and_consistency.R   | "How class imbalance and cross-set inconsistency affect prediction accuracy"|



Chapter 6
-------------------------------
Data: patient expression data files “bortezomib.patient.RData”, “erlotinib_data.RData”, 
and “pp.RData” should be downloaded from [here](http://compbio.cs.toronto.edu/cp2p/) and extracted in the folder /data

Scripts should be executed in the following order:

|     | name | description                                   |
|-----|------|-----------------------------------------------|
|1. |data_prep.R                  | preprocess all data necessary for the modelling |
|2. |main_analysis.R                       | performs the main anlysis |

Chapter 7
-------------------------------

|     | name | thesis's section/description                             |
|-----|------|------------------------------------------|
|1. | bl_unsupervised_analysis.R  | "Burkitt lymphoma drug sensitivity screen analysis" -- usupervised part |
|2. | bl_modelling.R  | "Burkitt lymphoma drug sensitivity screen analysis" -- modelling |
|3. | dkfz608_modelling.R  | "Drug response prediction model for DKFZ-608 compound" |
|4. | shiny_data.R | "Shiny application for complex drug response visualization" -- data preparation |
|5. | app.R   | "Shiny application for complex drug response visualization" -- app.R file of the application|

SessionInfo
------------------------

Session environment information 

#sessionInfo()

R version 3.5.0 (2018-04-23)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] PharmacoGx_1.8.3 e1071_1.7-0      doMC_1.3.5       iterators_1.0.9  foreach_1.4.4    caret_6.0-80     ggplot2_2.2.1   
[8] lattice_0.20-35 

loaded via a namespace (and not attached):
  [1] fgsea_1.4.1         colorspace_1.3-2    class_7.3-14        rprojroot_1.3-2     lsa_0.73.1          pls_2.6-0          
  [7] DRR_0.0.3           SnowballC_0.5.1     prodlim_2018.04.18  lubridate_1.7.4     codetools_0.2-15    splines_3.5.0      
 [13] mnormt_1.5-5        robustbase_0.93-0   knitr_1.20          RcppRoll_0.3.0      magicaxis_2.0.3     broom_0.4.4        
 [19] ddalpha_1.3.3       cluster_2.0.7-1     kernlab_0.9-26      sfsmisc_1.1-2       mapproj_1.2.6       compiler_3.5.0     
 [25] backports_1.1.2     assertthat_0.2.0    Matrix_1.2-14       lazyeval_0.2.1      limma_3.34.9        htmltools_0.3.6    
 [31] tools_3.5.0         bindrcpp_0.2.2      igraph_1.2.1        gtable_0.2.0        glue_1.2.0          RANN_2.5.1         
 [37] reshape2_1.4.3      dplyr_0.7.5         maps_3.3.0          fastmatch_1.1-0     Rcpp_0.12.17        slam_0.1-43        
 [43] Biobase_2.38.0      gdata_2.18.0        nlme_3.1-137        psych_1.8.4         timeDate_3043.102   gower_0.1.2        
 [49] stringr_1.3.1       gtools_3.5.0        DEoptimR_1.0-8      MASS_7.3-50         scales_0.5.0        ipred_0.9-6        
 [55] relations_0.6-8     RColorBrewer_1.1-2  sets_1.0-18         yaml_2.1.19         gridExtra_2.3       downloader_0.4     
 [61] rpart_4.1-13        stringi_1.2.3       NISTunits_1.0.1     plotrix_3.7-2       randomForest_4.6-14 caTools_1.17.1     
 [67] BiocGenerics_0.24.0 BiocParallel_1.12.0 lava_1.6.1          geometry_0.3-6      rlang_0.2.1         pkgconfig_2.0.1    
 [73] bitops_1.0-6        evaluate_0.10.1     pracma_2.1.4        purrr_0.2.5         bindr_0.1.1         recipes_0.1.3      
 [79] labeling_0.3        CVST_0.2-2          tidyselect_0.2.4    plyr_1.8.4          magrittr_1.5        R6_2.2.2           
 [85] gplots_3.0.1        dimRed_0.1.0        sm_2.2-5.5          pillar_1.2.3        foreign_0.8-70      withr_2.1.2        
 [91] survival_2.42-3     abind_1.4-5         nnet_7.3-12         tibble_1.4.2        KernSmooth_2.23-15  rmarkdown_1.10     
 [97] grid_3.5.0          data.table_1.11.4   marray_1.56.0       piano_1.18.1        ModelMetrics_1.1.0  digest_0.6.15      
[103] tidyr_0.8.1         stats4_3.5.0        munsell_0.5.0       celestial_1.4.1     magic_1.5-8         tcltk_3.5.0   