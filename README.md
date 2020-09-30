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
sig212_zscore <- t(scale(t(sig212[,2:(ncol(sig212)-3)]), center=T, scale=T)[,1:nrow(sig212)])
rownames(sig212_zscore) <- sig212$Gene
sigSPOP_zscore <- sig212_zscore

### 3. 
#>           Id   Hugo DHcR_Normal PDR_Normal GEXP_Normal    Reptime  PDR_Tumor
#> 1 SRR2069925 ABCB10 -0.12649901 -0.3666188   1.0297728 -0.6863157 0.06439566
#> 2 SRR2069925  ABCD3 -0.12649901 -0.4862479   1.6453271 -0.3921280 0.02900665
#> 3 SRR2069925   ABL2 -0.12649901 -0.3398278   2.1078829 -0.5868720 0.02704556
#> 4 SRR2069925  ACADM -0.12649901 -0.5184135   1.4140661  0.5608746 0.04907679
#> 5 SRR2069925  ACAP3  0.07178416 -0.1074514   0.2976902 -0.7028897 0.06073043
#> 6 SRR2069925  ACBD3 -0.12649901 -0.2668689   2.9805750 -0.9473556 0.06904483
#>   Depth_Tumor CpGs_Tumor DHcR_Tumor
#> 1    44.91818        110      0.000
#> 2    38.69118         68      0.000
#> 3    36.35714         28      0.000
#> 4    32.36364         33      0.000
#> 5    40.36500        200      0.025
#> 6    36.91892         37      0.000










