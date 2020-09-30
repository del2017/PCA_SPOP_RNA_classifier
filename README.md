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

### 2. SPOP mutant signature input from TCGA FPKM expression data
sig212 <- read.table("del2017/SPOP-RNA-classifier/TCGA_333_SPOP_sig_212genes.txt", sep="\t", header=T, check.names=F)

sig212[1:5, 1:5]
    Gene TCGA-2A-A8VL-01 TCGA-2A-A8VO-01 TCGA-2A-A8VT-01 TCGA-2A-A8VV-01
1  ABCC1       2203.9135       1360.4809       2473.8971       2194.1405
2 ABHD11       1522.1421       1391.7181       1814.3382       1170.9613
3 ABLIM3        234.2945        820.1434         92.2794        251.7943
4  ACAD9        886.2925        913.1505        816.5368        918.2822
5  ADAM7          1.0299          0.0000          5.5147         43.2992











