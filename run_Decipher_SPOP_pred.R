library(survival)
library(e1071)
library(RColorBrewer)
library(gdata)
library(survminer)

#Subclass color
sur_col_SPOPo<-brewer.pal(9,"Set1")[5]
sur_col_wt<-brewer.pal(9,"Set1")[3]


#Decipher data input
#dat_gdx_expm<-readRDS("/Users/deli/Desktop/Chris/TCGA-PCA/GenomeDX/gdx_20180530_gdxExon_expression-gene.rds")
dat_gdx_expm[1:3, 1:3]
#Decipher data info saved in WCM cluster athana:
datinfo<-read.xls("/home/del2017/lab_del2017/Decipher/gndx-2016-04-26_clinical_expression_data.xlsx", 1) 
datinfo2<-read.xls("/home/del2017/lab_del2017/Decipher/gndx-2016-04-26_clinical_data_key.xlsx", 1)

#ETS exp
dat_gdx_expm_sub<-t(merge(dat_gdx_expm, data.frame(c("ERG", "ETV1", "ETV4", "ETV5", "FLI1", "SPINK1")), by.x="gene", by.y="c..ERG....ETV1....ETV4....ETV5....FLI1....SPINK1.."))
dat_gdx_expm_sub<-data.frame(dat_gdx_expm_sub[2:nrow(dat_gdx_expm_sub),])
colnames(dat_gdx_expm_sub)<-c("ERG","ETV1","ETV4","ETV5","FLI1","SPINK1")
#Molecular subtype ERG, ETS, SPINK1 and TripleNeg
datinfom<-merge(datinfo, dat_gdx_expm_sub, by.x="celfile_name", by.y="row.names", sort=F)
#Exclude 116 duplicates cases from JHMI cohort
datinfom<-unique(datinfom[,c(1,3:(ncol(datinfom)-9), (ncol(datinfom)-7):ncol(datinfom))])
datinfom$ERG_overexp<-ifelse(as.numeric(as.character(datinfom$ERG))>0.6, "Pos", "Neg")
datinfom$ETS_overexp<-ifelse(as.numeric(as.character(datinfom$ETV1))>0.41| as.numeric(as.character(datinfom$ETV4))>0.32|as.numeric(as.character(datinfom$ETV5))>0.48|as.numeric(as.character(datinfom$FLI1))>0.53, "Pos", "Neg")
datinfom$SPINK1_overexp<-ifelse(as.numeric(as.character(datinfom$SPINK1))>1.03, "Pos", "Neg")
datinfom$TripleNeg<-ifelse(as.numeric(as.character(datinfom$ERG))<=0.6 & as.numeric(as.character(datinfom$ETV1))<=0.41& as.numeric(as.character(datinfom$ETV4))<=0.32&as.numeric(as.character(datinfom$ETV5))<=0.48&as.numeric(as.character(datinfom$FLI1))<=0.53
                           &as.numeric(as.character(datinfom$SPINK1))<=1.03,"Pos", "Neg")

#TCGA PCA 333 sample information
dinfo_TCGA <- read.csv("/Users/deli/Dropbox/Deli_LabMeeting/SPOP_signature/SPOP_prediction/TCGA_PRAD_333_dinfo.csv", check.names = F)
dinfo_TCGA_SPOP<-dinfo_TCGA[, c("SAMPLE_ID", "Sample.ID", "SPOP_mut")]

#SPOP sig 212 genes
sig212<-read.table("/Users/deli/Dropbox/Deli_LabMeeting/SPOP_signature/SPOP_prediction/TCGA_333_SPOP_sig_212genes.txt", sep="\t", header=T, check.names=F)
sig212_zscore<-t(scale(t(sig212[,2:(ncol(sig212)-3)]), center=T, scale=T)[,1:nrow(sig212)])
rownames(sig212_zscore)<-sig212$Gene
sigSPOP_zscore <- sig212_zscore

#tesing data from GenomeDx after Z-score scale from each cohort 
dinfo_gdx<-datinfom
#SPOP
dat_gdx_SPOP<-merge(dat_gdx_expm, data.frame(rownames(sigSPOP_zscore)), by.x="gene", by.y="rownames.sigSPOP_zscore.")
dat_gdx_SPOP<-data.frame(dat_gdx_SPOP$gene, data.frame(apply(dat_gdx_SPOP[,2:ncol(dat_gdx_SPOP)], 2, function(x) as.numeric(as.character(x))), check.rows = F), check.rows = F)
colnames(dat_gdx_SPOP)<-colnames(dat_gdx_expm)

cost_value_SPOP<-0.03

#Identify the subtypes from each cohort
svm.pred.info_clas_cohort<-NULL
for (i in 1: length(summary(datinfom$study_name)))
{
  datinfom_i<-datinfom[datinfom$'study_name'==names(summary(datinfom$study_name)[i]),]
  
  #SVM model to predict SPOP
  dat_gdx_SPOP_i<-cbind(dat_gdx_SPOP$gene, dat_gdx_SPOP[, as.matrix(datinfom[datinfom$'study_name'==names(summary(datinfom$study_name)[i]),1])])
  dat_gdx_SPOP_i_zscore<-t(scale(t(dat_gdx_SPOP_i[,2:ncol(dat_gdx_SPOP_i)]), center=T, scale=T)[,1:nrow(dat_gdx_SPOP_i)])
  rownames(dat_gdx_SPOP_i_zscore)<-dat_gdx_SPOP_i[,1]
  dat_gdx_SPOP_i_zscore<-dat_gdx_SPOP_i_zscore[complete.cases(dat_gdx_SPOP_i_zscore),]
  datm<-merge(sigSPOP_zscore, dat_gdx_SPOP_i_zscore, by="row.names")
  #traning data from TCGA 
  traindat<-t(datm[,2:(ncol(sigSPOP_zscore)+1)])
  colnames(traindat)<-datm[,1]
  traindat<-merge(traindat, dinfo_TCGA_SPOP[,c(1,ncol(dinfo_TCGA_SPOP))] , by.x="row.names", by.y="SAMPLE_ID")
  traindatm<-traindat[,2:ncol(traindat)]
  rownames(traindatm)<-rownames(traindat)
  #testing data from GenomeDx 
  testdat<-t(datm[,(ncol(sigSPOP_zscore)+2):ncol(datm)])
  colnames(testdat)<-datm[,1]
  testdatm<-testdat
  # svm predict linear cost=1 type="C-classification"
  svm.model_cost_SPOP <- svm(SPOP_mut ~ ., data = traindatm,  kernel="linear",  type="C-classification", cost=cost_value_SPOP)
  svm.pred_clas_SPOP <- predict(svm.model_cost_SPOP, testdatm)
  
  #Combine prediction
  svm.pred.info_clas_i<-merge(datinfom_i, testdatm[,1], by.x="celfile_name", by.y="row.names", sort=F)
  svm.pred.info_clas_i$pred_SPOP<-as.matrix(svm.pred_clas_SPOP)
  svm.pred.info_clas_cohort<-rbind(svm.pred.info_clas_cohort, svm.pred.info_clas_i)
}


svm.pred.info_clas_cohort$ERG_overexp_status<-ifelse(svm.pred.info_clas_cohort$ERG_overexp=="Pos", 1, 0)
svm.pred.info_clas_cohort$ETS_overexp_status<-ifelse(svm.pred.info_clas_cohort$ETS_overexp=="Pos", 1, 0)

#Call the ERG/ETS first, then SPOP
svm.pred.info_clas_cohort$ERG_status<-ifelse(svm.pred.info_clas_cohort$ERG_overexp=="Pos" & svm.pred.info_clas_cohort$ETS_overexp=="Neg", "Pos", "Neg")
svm.pred.info_clas_cohort$ETS_status<-ifelse(svm.pred.info_clas_cohort$ETS_overexp=="Pos" & svm.pred.info_clas_cohort$ERG_overexp=="Neg", "Pos", "Neg")
svm.pred.info_clas_cohort$SPOPmut<-ifelse(svm.pred.info_clas_cohort$pred_SPOP==1 & (svm.pred.info_clas_cohort$ERG_status=="Neg" & svm.pred.info_clas_cohort$ETS_status=="Neg"), "Pos", "Neg")

nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ERG_status=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$ETS_status=="Pos",])
nrow(svm.pred.info_clas_cohort[svm.pred.info_clas_cohort$SPOPmut=="Pos",])


svm.pred.info_clas_cohort_sur<-svm.pred.info_clas_cohort
svm.pred.info_clas_cohort_sur$Subtype<-ifelse(svm.pred.info_clas_cohort_sur$SPOPmut=="Pos", "SPOP_mut", "SPOP_wt")
svm.pred.info_clas_cohort_sur_fit<-survfit( Surv(met_time, met) ~ Subtype, data = svm.pred.info_clas_cohort_sur)
ggsurvplot(svm.pred.info_clas_cohort_sur_fit,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c(sur_col_SPOPo, sur_col_wt),
           break.time.by = 50, 
           xlim = c(0, 250),
           xlab="MET time (month)")


