library(data.table)
library(dplyr)

# In the meta-analysis, only a Z-score is computed.
# Some softwares to do post-GWAS analyses need BETA and SE estimates.
# BETA and SE are estimated from Z-score, A1_FREQ and TotalSampleSize according to https://www.nature.com/articles/ng.3538#methods .

#------------------------
# read files
#------------------------
ukbiobank=fread("~/results_GWAS_imaging_WMH_15.10.2018/GWAS_logWMHnorm_ALLchr.logWMHnorm_info.glm.linear.chrpos", header=TRUE) # contains Z-score and TotalSampleSize
metaanalysis=fread("~/results_GWAS_imaging_WMH_15.10.2018/metaGWAS_logWMHnorm_chrpos.tbl", header=TRUE) # contains A1_freq

#------------------------
# Data processing of the GWAS in UK Biobank
#------------------------
# doubling rows in ukbiobank GWAS to get the combinations A1, AX and AX, A1 to merge with metaGWAS
ukbiobank2=bind_rows(ukbiobank%>%select(ID,A1,AX,A1_FREQ),
                     ukbiobank%>%mutate(Allele1=AX, AlleleX=A1)%>%
                       select(ID,Allele1,AlleleX,A1_FREQ)%>%
                       rename(A1=Allele1, AX=AlleleX))%>%
  mutate(ID_A1_AX=paste(ID,A1,AX, sep="_"))
tabukbiobank2=count(ukbiobank2, ID_A1_AX)
tabukbiobank2_variants=filter(tabukbiobank2, n>1)

ukbiobank3=ukbiobank%>%
  mutate(ID_A1_AX=paste(ID,A1,AX, sep="_"), ID_AX_A1=paste(ID,AX,A1, sep="_"))%>%
  filter((!ID_A1_AX%in%tabukbiobank2_variants$ID_A1_AX) & (!ID_AX_A1%in%tabukbiobank2_variants$ID_A1_AX))

#------------------------
# Data processing of the meta-analysis
#------------------------
# setting allele1 and allele2 to upper cases
metaanalysis=metaanalysis%>%
  mutate(Allele1up=toupper(Allele1),
         Allele2up=toupper(Allele2))

#------------------------
# Merging the 2 datasets
#------------------------
merged1=inner_join(ukbiobank3, metaanalysis, by=c("ID"="MarkerName", "A1"="Allele1up", "AX"="Allele2up"))
merged2=inner_join(ukbiobank3, metaanalysis, by=c("ID"="MarkerName", "A1"="Allele2up", "AX"="Allele1up"))
merged=bind_rows(merged1, merged2)

#------------------------
# Computing BETA and SE
#------------------------
merged_modified=merged%>%
  mutate(SEmeta=1/sqrt(2*A1_FREQ*(1-A1_FREQ)*(TotalSampleSize + Zscore^2)))%>%
  mutate(BETAmeta=Zscore*SEmeta)%>%
  mutate(BETAmeta_A1=sign(BETA)*abs(BETAmeta),
         Zscoremeta_A1=sign(BETA)*abs(Zscore))

#------------------------
# Writing the dataset
#------------------------
merged_final=merged_modified%>%
  select(ID,
         A1,
         AX,
         A1_FREQ,
         Weight,
         `P-value`,
         TotalSampleSize,
         CHROM.y,
         POS.y,
         Zscoremeta_A1,
         SEmeta,
         BETAmeta_A1)%>%
  rename(CHROM=CHROM.y,
         POS=POS.y
  )
write.table(merged_final, "~/results/results_GWAS_imaging_WMH_15.10.2018/metaGWAS_logWMHnorm_BETA_SE_A1_FREQ.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
