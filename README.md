# White_Clover_WSC_Outlier_Detection_GWAS

This repository contains the scripts and workflows used to analyse phenotypic data, perform population structure assessment, 
outlier detection and genome-wide association studies for Pearson et al. (2022) Outlier analyses and genome-wide association study identify 
*glgC* and *ERD6-like 4* as candidate genes for foliar water-soluble carbohydrate accumulation in *Trifolium repens*.


## Estimated phenotype means

Assessment of the WSC phenotype and average leaf area using linear models accounting for treatment and spatial variation.

[SSS_EMM.R](https://github.com/SofiePearson/White_Clover_WSC_Outlier_Detection_GWAS/blob/main/CODE/SSS_EMM.R)

[LA_EMM.R](https://github.com/SofiePearson/White_Clover_WSC_Outlier_Detection_GWAS/blob/main/CODE/LA_EMM.R)


## Correlation and regression analyses

Correlation and regression analyses of the phenotypic data.

[Correlation_SSS_LA.R](https://github.com/SofiePearson/White_Clover_WSC_Outlier_Detection_GWAS/blob/main/CODE/Correlation_SSS_LA.R)

[Regression_SSS_LA.R](https://github.com/SofiePearson/White_Clover_WSC_Outlier_Detection_GWAS/blob/main/CODE/Regression_SSS_LA.R)


## Population structure

PCA assessment of genetic data and controls, Pairwise F<sub>ST</sub> analysis, DAPC determination of genetic clusters and AMOVA.

[Population_genetic_structure.R](https://github.com/SofiePearson/White_Clover_WSC_Outlier_Detection_GWAS/blob/main/CODE/Population_genetic_structure.R)


## Outlier detection

PCAdapt, Bayescan and KGD-F<sub>ST</sub> analyses.

[Outlier_detection.R](https://github.com/SofiePearson/White_Clover_WSC_Outlier_Detection_GWAS/blob/main/CODE/Outlier_detection.R)


## GWAS

Genome-wide association study for a subset of 605 samples with genotype and phenotype data.

[GWAS.R](https://github.com/SofiePearson/White_Clover_WSC_Outlier_Detection_GWAS/blob/main/CODE/GWAS.R)
