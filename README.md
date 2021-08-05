# PCA_SPOP_RNA_classifier
A novel gene expression signature, classifier (SCaPT), and decision tree to predict the SPOP mutant subclass from RNA gene expression data.

## Publication
Impact of the SPOP Mutant Subtype on the Interpretation of Clinical Parameters in Prostate Cancer

https://ascopubs.org/doi/full/10.1200/PO.18.00036

## Abstract
Purpose
Molecular characterization of prostate cancer, including The Cancer Genome Atlas, has revealed distinct subtypes with underlying genomic alterations. One of these core subtypes, SPOP (speckle-type POZ protein) mutant prostate cancer, has previously only been identifiable via DNA sequencing, which has made the impact on prognosis and routinely used risk stratification parameters unclear.

Methods
We have developed a novel gene expression signature, classifier (Subclass Predictor Based on Transcriptional Data), and decision tree to predict the SPOP mutant subclass from RNA gene expression data and classify common prostate cancer molecular subtypes. We then validated and further interrogated the association of prostate cancer molecular subtypes with pathologic and clinical outcomes in retrospective and prospective cohorts of 8,158 patients.

Results
The subclass predictor based on transcriptional data model showed high sensitivity and specificity in multiple cohorts across both RNA sequencing and microarray gene expression platforms. We predicted approximately 8% to 9% of cases to be SPOP mutant from both retrospective and prospective cohorts. We found that the SPOP mutant subclass was associated with lower frequency of positive margins, extraprostatic extension, and seminal vesicle invasion at prostatectomy; however, SPOP mutant cancers were associated with higher pretreatment serum prostate-specific antigen (PSA). The association between SPOP mutant status and higher PSA level was validated in three independent cohorts. Despite high pretreatment PSA, the SPOP mutant subtype was associated with a favorable prognosis with improved metastasis-free survival, particularly in patients with high-risk preoperative PSA levels.

Conclusion
Using a novel gene expression model and a decision tree algorithm to define prostate cancer molecular subclasses, we found that the SPOP mutant subclass is associated with higher preoperative PSA, less adverse pathologic features, and favorable prognosis. These findings suggest a paradigm in which the interpretation of common risk stratification parameters, particularly PSA, may be influenced by the underlying molecular subtype of prostate cancer.


## Install
The current R version is 3.6.1. 

Required SVM package:

install.packages("https://cran.r-project.org/src/contrib/Archive/e1071/e1071_1.7-2.tar.gz", repos=NULL)

## Prerequisite

### 1. SPOP mutant signature
Signature is downloaded from Supplementary Table 1a (https://ascopubs.org/doi/suppl/10.1200/PO.18.00036/)

### 2. SPOP signature normalization based on TCGA FPKM expression data
sig212 <- read.table("del2017/SPOP-RNA-classifier/TCGA_333_SPOP_sig_212genes.txt", sep="\t", header=T, check.names=F)

### 3. SPOP mutant status is derived from TCGA PCA study (PMID: 26544944)
Table S1: https://www.cell.com/fulltext/S0092-8674(15)01339-2#supplementaryMaterial
