args=commandArgs(trailingOnly=TRUE)

rootname=args[1]
trait=args[2]
relatedness=args[3]
QC=args[4]
fam=args[5]
OUT_DIR=args[6]

library(data.table)
library(ukbtools)
library(dplyr)

#-----------------------------------------------------------------------
# IMPORT
#-----------------------------------------------------------------------
sample_IDs=read.table(paste(OUT_DIR,"/samples/",rootname, "_", trait, ".samples.txt", sep=""), header=F)
sample_IDs=sample_IDs[,1]
rel_data=fread(relatedness, header=TRUE)
QC_data=fread(QC, header=FALSE)
FAM_data=fread(fam, header=FALSE)


#-----------------------------------------------------------------------
# INDIVIDUALS WITH WMH AND BRAIN SIZE 
#-----------------------------------------------------------------------
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

QC_data=cbind.data.frame(FAM_data, QC_data)

QC_sub=QC_data %>%
  filter(IID %in% sample_IDs)

#-----------------------------------------------------------------------
# QC
#-----------------------------------------------------------------------
# Relatedness
RemovedRelated=ukb_gen_samples_to_remove(data=rel_data, ukb_with_data=QC_sub[["IID"]], cutoff=0.0884)
nbRelated=length(RemovedRelated)

# Gender mismatch
QC_sub[["Submitted.Gender"]]=as.factor(QC_sub[["Submitted.Gender"]])
QC_sub[["Inferred.Gender"]]=as.factor(QC_sub[["Inferred.Gender"]])

print(summary(QC_sub[["Inferred.Gender"]]))
print(summary(QC_sub[["Submitted.Gender"]]))

RemovedGender=QC_sub %>%
  filter(Submitted.Gender != Inferred.Gender) %>%
  select(IID)
nbGender=nrow(RemovedGender)

# Outliers in heterozygosity and missingness
QC_sub[["het.missing.outliers"]]=as.factor(QC_sub[["het.missing.outliers"]])
print(summary(QC_sub[["het.missing.outliers"]]))
RemovedOutliers=QC_sub %>%
  filter(het.missing.outliers == 1) %>%
  select(IID)
nbOutliers=nrow(RemovedOutliers)

# Missingness
RemovedMissingness=QC_sub %>%
  filter(sample.qc.missing.rate > 0.05) %>%
  select(IID)
nbMissingness=nrow(RemovedMissingness)

# White British ancestry
QC_sub[["in.white.British.ancestry.subset"]]=as.factor(QC_sub[["in.white.British.ancestry.subset"]])

RemovedBritAncestry=QC_sub %>%
  filter(!IID %in% unique(c(RemovedRelated,
                            RemovedGender[["IID"]],
                            RemovedOutliers[["IID"]],
                            RemovedMissingness[["IID"]]))) %>%
  filter(in.white.British.ancestry.subset == 0) %>%
  select(IID)
nbBritAncestry=nrow(RemovedBritAncestry)

# European ancestry - Joni's Code
##Set seed
set.seed(1204688)

##K means clustering on each PC
K_MEAN<-4

PCs=QC_data %>%
  select(IID, PC1, PC2)
PC1_K<-kmeans(PCs$PC1, K_MEAN)
PC2_K<-kmeans(PCs$PC2, K_MEAN)
PC1_2_K<-kmeans(PCs[,c("PC1","PC2")], K_MEAN)

##Add clusters to PC dataframe
PCs$PC1.Cluster<-PC1_K$cluster
PCs$PC2.Cluster<-PC2_K$cluster
#PCs$Clusters<-as.factor(paste(PC1_K$cluster,PC2_K$cluster,sep="."))

##WWE group is the majority
MAX_PC1<-ifelse(match(max(table(PCs$PC1.Cluster, PCs$PC2.Cluster)), table(PCs$PC1.Cluster, PCs$PC2.Cluster)) %% K_MEAN == 0, K_MEAN, match(max(table(PCs$PC1.Cluster, PCs$PC2.Cluster)), table(PCs$PC1.Cluster, PCs$PC2.Cluster)) %% K_MEAN)
MAX_PC2<-ceiling(match(max(table(PCs$PC1.Cluster, PCs$PC2.Cluster)), table(PCs$PC1.Cluster, PCs$PC2.Cluster))/K_MEAN)

##Make list of WWE IDs
WWE<-as.data.frame(PCs[PCs$PC1.Cluster == MAX_PC1 & PCs$PC2.Cluster == MAX_PC2,1])
names(WWE)<-"ID"

RemovedEurAncestry=PCs %>%
  filter(IID %in% QC_sub[["IID"]]) %>%
  filter(!IID %in% unique(c(RemovedRelated,
                            RemovedGender[["IID"]],
                            RemovedOutliers[["IID"]],
                            RemovedMissingness[["IID"]]))) %>%
  filter(PC1.Cluster != MAX_PC1 | PC2.Cluster != MAX_PC2)%>%
  select(IID)					
nbEurAncestry=nrow(RemovedEurAncestry)

#-----------------------------------------------------------------------
# LIST OF REMOVED INDIVIDUALS
#-----------------------------------------------------------------------
RemovedAncestry=RemovedEurAncestry
nbAncestry=nrow(RemovedAncestry)

# All removed individuals
RemovedIndividuals=unique(c(
  RemovedRelated,
  RemovedGender[["IID"]],
  RemovedOutliers[["IID"]],
  RemovedMissingness[["IID"]],
  RemovedAncestry[["IID"]]
))
nbTotal=length(RemovedIndividuals)

ExtractedIndividuals=QC_sub %>%
  filter(! IID %in% RemovedIndividuals) %>%
  select(IID)
ExtractedIndividuals=ExtractedIndividuals[["IID"]]


#-----------------------------------------------------------------------
# PLOT OF PCs
#-----------------------------------------------------------------------
couleurs=rep("all",nrow(QC_data))
couleurs[which(QC_data[["IID"]] %in% ExtractedIndividuals)]="WMH_1"
couleurs=as.factor(couleurs)

library(ggplot2)

RemovedIndividuals_noancestry=unique(c(
  RemovedRelated,
  RemovedGender[["IID"]],
  RemovedOutliers[["IID"]],
  RemovedMissingness[["IID"]]
))

png(paste(OUT_DIR, "/pop_structure/plot_PC_",trait,"_keptindividuals.png",sep=""), width=20, height=20, units="cm", res=300)
ggplot(data=QC_data)  +
  geom_point(aes(x=PC1, y=PC2, col=couleurs), size=0.5)+
  geom_point(data=filter(QC_data, IID %in% ExtractedIndividuals),aes(x=PC1, y=PC2, col="WMH_1"), size=0.5)+
  scale_color_manual(values=c(adjustcolor(c("grey"), alpha=0.2),  "tomato"), breaks=c("all",  "WMH_1"))
dev.off()

png(paste(OUT_DIR, "/pop_structure/plot_PC_",trait,"_keptindividuals2.png",sep=""), width=20, height=20, units="cm", res=300)
ggplot(data=QC_data)  +
  geom_point(aes(x=PC1, y=PC2, col=couleurs), size=0.5)+
  geom_point(data=filter(QC_data, IID %in% ExtractedIndividuals),aes(x=PC1, y=PC2, col="WMH_1"), size=0.5)+
  geom_point(data=filter(QC_data, IID %in% filter(RemovedAncestry,! IID %in% RemovedIndividuals_noancestry)[["IID"]]),aes(x=PC1, y=PC2, col="WMH_0"), size=0.5)+
  scale_color_manual(values=c(adjustcolor(c("grey"), alpha=0.2), "black", "tomato"),  breaks=c("all", "WMH_0" ,"WMH_1"))
dev.off()


#-----------------------------------------------------------------------
# SUMMARY OF ALL REMOVED INDIVIDUALS
#-----------------------------------------------------------------------
tab=matrix(
  c(
    "Initial number of samples", length(sample_IDs),
    "Number of samples in QC file", nrow(QC_sub),
    "Kinship >= 0.0884",nbRelated,
    "Gender mismatch", nbGender,
    "Outliers (het. missing.)", nbOutliers,
    "Missing rate > 0.05", nbMissingness,
    "No European ancestry", nbAncestry,
    "Total of removed samples", nbTotal,
    "Remaining samples", nrow(QC_sub) - nbTotal
  ),
  ncol=2,
  byrow=TRUE
)


#-----------------------------------------------------------------------
# TXT OUTPUT
#-----------------------------------------------------------------------
write.table(tab, paste(OUT_DIR, "/samples/summarytab_removedsamples_",trait,".txt",sep=""), quote=F, row.names=F, col.names=F, sep="\t")
write.table(ExtractedIndividuals, paste(OUT_DIR, "/samples/extractedsamples_",trait,".txt",sep=""), col.names=F, quote=F, row.names=F)
write.table(RemovedIndividuals, paste(OUT_DIR, "/samples/removedsamples_",trait,".txt",sep=""), col.names=F, quote=F, row.names=F)
