args=commandArgs(trailingOnly=TRUE)

TRAIT = args[1]
RDA_FILE = args[2]
QC_FILE = args[3]
FAM_FILE = args[4]
PCA_FILE = args[5]
OUT_DIR = args[6] 

library(dplyr)
library(data.table)
library(mice)

#-------------------------------------------------
# IMPORT RDA_FILE
#-------------------------------------------------

load(RDA_FILE)

# Interesting columns

InterestingColumns=list(
  ## Main columns
  IID="eid",
  WMH="f25781_2_0",
  BrainSize="f25009_2_0",
  
  ## Covariates
  AgeBaseline="f21022_0_0",
  Sex="f31_0_0",
  YearBirth="f34_0_0",
  MonthBirth="f52_0_0",
  DateAssessment="f53_2_0",
  UKBcentre="f54_2_0",
  MeanrfMRI="f25741_2_0",
  MeantfMRI="f25742_2_0",
  
  ## Correlation
  SysBloodPress_auto="f4080_2_",
  DiasBloodPress_auto="f4079_2_"
)

Classes=list(
  ## Main columns
  IID="character",
  WMH="numeric",
  BrainSize="numeric",
  
  ## Covariates
  AgeBaseline="numeric",
  Sex="factor",
  YearBirth="numeric",
  MonthBirth="numeric",
  DateAssessment="character",
  UKBcenter="factor",
  MeanrfMRI="numeric",
  MeantfMRI="numeric",
  
  ## Correlation
  SysBloodPress_auto=rep("numeric",2),
  DiasBloodPress_auto=rep("numeric",2)
)

ColumnsIndices=lapply(InterestingColumns,function(x) grep(pattern=x, colnames(res_ukb)))

# Extract interesting columns
res_ukb_sub=res_ukb %>%
  select(unlist(ColumnsIndices) )
colnames(res_ukb_sub)=names(unlist(ColumnsIndices))

# Recode with the correct class
RecodeVariables=function(data, classes){
  for(i in 1:length(classes)){
    if(classes[i]=="numeric"){
      data[,i]=as.numeric(as.character(data[,i]))
    }else if(classes[i]=="character"){
      data[,i]=as.character(data[,i])
    }else if(classes[i]=="factor"){
      data[,i]=as.factor(data[,i])
    }
  }
  return(data)
}

res_ukb_sub=RecodeVariables(res_ukb_sub, unlist(Classes))

#-------------------------------------------------
# IMPORT PCA_FILE
#-------------------------------------------------
# Extract samples from PCA
pca=fread(PCA_FILE, header=TRUE)
colnames(pca)=c("FID",
                "IID", 
                paste("PC", 1:10 ,sep="_"))
pca = pca %>%
  select(IID,FID, eval(paste("PC", 1:10, sep="_")))
pca$IID=as.character(pca$IID)
for(i in grep( colnames(pca), pattern="PC")){
  pca[[i]]=round(pca[[i]], digits=4)
}

data2=res_ukb_sub %>%
  right_join(pca, by= "IID")


#-------------------------------------------------
# IMPORT QC_FILE
#-------------------------------------------------
QC_data=fread(QC_FILE, header=FALSE)
FAM_data=fread(FAM_FILE, header=FALSE)

FAM_Columns=c(
  "FID",
  "IID",
  "PID",
  "MID",
  "Sex",
  "Batch"
)

QC_Columns=c(
  paste("genotyping.array", 1:2, sep="_"),
  paste("Batch", 1:2, sep="_"),
  "Plate.Name",
  "Well",
  "Cluster.CR",
  "dQC",
  "Internal.Pico..ng.uL.",
  "Submitted.Gender",
  "Inferred.Gender",
  "X.intensity",
  "Y.intensity",
  "Submitted.Plate.Name",
  "Submitted.Well",
  "sample.qc.missing.rate",
  "heterozygosity",
  "heterozygosity.pc.corrected",
  "het.missing.outliers",
  "putative.sex.chromosome.aneuploidy",
  "in.kinship.table",
  "excluded.from.kinship.inference",
  "excess.relatives",
  "in.white.British.ancestry.subset",
  "used.in.pca.calculation",
  paste("PC", 1:40, sep=""),
  "in.Phasing.Input.chr1_22",
  "in.Phasing.Input.chrX",
  "in.Phasing.Input.chrXY"
)

colnames(FAM_data)=FAM_Columns
colnames(QC_data)=QC_Columns
QC_data=QC_data %>%
  bind_cols(FAM_data)

QC_data=QC_data %>%
  select(IID, Batch_1)
QC_data[["IID"]]=as.character(QC_data[["IID"]])
QC_data[["Batch_1"]]=as.factor(QC_data[["Batch_1"]])

output=data2 %>%
  left_join(QC_data, by= "IID")

output=output%>%
  select(FID,IID,everything())

#-------------------------------------------------
# TRANSFORMATION OF VARIABLES
#-------------------------------------------------
# Phenotype : Normalize by BrainSize
output=output %>%
  mutate(logWMHnorm=log(WMH)*mean(BrainSize)/BrainSize)

# Sex
output$SexCode=factor(output$Sex, levels=c("0", "1"))
levels(output$SexCode)=c("Female", "Male")

# Age at MRI : Compute from DateAssessment, MonthBirth, YearBirth
output$DateAssessment=as.Date(output$DateAssessment, format="%Y-%m-%d")
output$Birthday=unlist(lapply(1:nrow(output), function(x) paste(output$YearBirth[x],output$MonthBirth[x],1, sep="-")))
output$Birthday=as.Date(output$Birthday, format="%Y-%m-%d")

output=output%>%
  mutate(	AgeMRI=as.numeric(difftime(DateAssessment,Birthday, unit="days"))/365, 
          AgeMRI2=AgeMRI^2,
          AgeMRI3=AgeMRI^3
  )

# UKB centre
output$UKBcentre=factor(output$UKBcentre, levels=c("11025", "11027"))
levels(output$UKBcentre)=c("Cheadle_imaging", "Newcastle_imaging")
output$UKBcentre_name=output$UKBcentre
levels(output$UKBcentre)=c("0", "1")

# Batch
output$Batch_name = output$Batch_1
levels(output$Batch_1)=c("0", "1")

# MeanrfMRI, MeantfMRI: Imputation of missing values
subimputation=output %>%
  select(Sex,
         UKBcentre,
         MeanrfMRI, 
         MeantfMRI,
         eval(paste("PC_",1:10,sep="")), 
         AgeMRI, 
         Batch_1)
imputeddata=mice(subimputation, method="pmm")

output[is.na(output$MeanrfMRI),"MeanrfMRI"]=imputeddata$imp$MeanrfMRI[,1]
output[is.na(output$MeantfMRI),"MeantfMRI"]=imputeddata$imp$MeantfMRI[,1]

# Interactions
output = output %>%
  mutate(	AgeMRI_Sex=AgeMRI*as.numeric(as.character(Sex)),
          AgeMRI2_Sex=AgeMRI2*as.numeric(as.character(Sex)), 
          AgeMRI3_Sex=AgeMRI3*as.numeric(as.character(Sex))
  )

#-------------------------------------------------
# OUPTPUT
#-------------------------------------------------
write.table(output, paste(OUT_DIR, "/", TRAIT, "_covar.txt", sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")



