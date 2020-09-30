# SPOP-RNA-classifier
A novel gene expression signature, classifier (SCaPT), and decision tree to predict the SPOP mutant subclass from RNA gene expression data.

## Publication
"Impact of the SPOP Mutant Subtype on the Interpretation of Clinical Parameters in Prostate Cancer"

https://ascopubs.org/doi/full/10.1200/PO.18.00036

## Install
The current R version is 3.6.1. 

Required SVM package:

install.packages("https://cran.r-project.org/src/contrib/Archive/e1071/e1071_1.7-2.tar.gz", repos=NULL)

## Prerequisite

### 1. SPOP mutant signature
Signature is downloaded from Supplementary Table 1a (https://ascopubs.org/doi/suppl/10.1200/PO.18.00036/)

### 2. SPOP mutant signature input normalization based on TCGA FPKM expression data
sig212 <- read.table("del2017/SPOP-RNA-classifier/TCGA_333_SPOP_sig_212genes.txt", sep="\t", header=T, check.names=F)

