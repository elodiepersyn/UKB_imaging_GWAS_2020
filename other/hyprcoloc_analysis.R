library(hyprcoloc)
library(data.table)
library(dplyr)

# Summary statistics
# run th other/add_BETA_SE_metaGWAS.R to get a file for the meta-analysis with BETA and SE estimates.
WMH_sumstats=fread("~/results/results_GWAS_imaging_WMH_15.10.2018/metaGWAS_logWMHnorm_BETA_SE_A1_FREQ.txt", header=TRUE)
FA_sumstats=fread("~/results/results_GWAS_imaging_FA_15.10.2018/GWAS_PC1_trait_ALLchr.PC1_trait_info.glm.linear", header=TRUE)
MD_sumstats=fread("~/results/results_GWAS_imaging_MD_15.10.2018/GWAS_PC1_trait_ALLchr.PC1_trait_info.glm.linear", header=TRUE)
AS_sumstats=fread("~/data/MEGASTROKE_data/MEGASTROKE.1.AS.TRANS.out", header=TRUE)
SVS_sumstats=fread("~/data/MEGASTROKE_data/MEGASTROKE.5.SVS.TRANS.out", header=TRUE)
AIS_sumstats=fread("~/data/MEGASTROKE_data/MEGASTROKE.2.AIS.TRANS.out", header=TRUE)
LAS_sumstats=fread("~/data/MEGASTROKE_data/MEGASTROKE.3.LAS.TRANS.out", header=TRUE)
CES_sumstats=fread("~/data/MEGASTROKE_data/MEGASTROKE.4.CES.TRANS.out", header=TRUE)

# Sentinels
# run the other/sentinel_proxy_files.R to get these files
sentinelsWMH=read.table("~/results/results_GWAS_imaging_WMH_15.10.2018/sentinels.txt", header=TRUE)
sentinelsFA=read.table("~/results/results_GWAS_imaging_FA_15.10.2018/sentinels.txt", header=TRUE)
sentinelsMD=read.table("~/results/results_GWAS_imaging_MD_15.10.2018/sentinels.txt", header=TRUE)
sentinels=bind_rows(sentinelsWMH,sentinelsFA, sentinelsMD)

results=list()
for(snp in 1:nrow(sentinels)){
  # For SNP 1
  sentinel=sentinels$rsID[snp]
  chr=sentinels$CHR[snp]
  position=sentinels$START[snp]
  
  subsetWMH=WMH_sumstats%>%
    filter(CHROM==chr & (POS>=(position-5e5) & POS<=(position+5e5)))%>%
    dplyr::select(ID, A1, AX, BETAmeta_A1, SEmeta, POS, CHROM)%>%
    dplyr::rename(A1_WMH=A1, AX_WMH=AX, BETA_WMH=BETAmeta_A1, SE_WMH=SEmeta)
  
  subsetFA=FA_sumstats%>%
    filter(ID %in% subsetWMH$ID)%>%
    dplyr::select(ID, A1, AX, BETA, SE)%>%
    dplyr::rename(A1_FA=A1, AX_FA=AX, BETA_FA=BETA, SE_FA=SE)
  
  subsetMD=MD_sumstats%>%
    filter(ID %in% subsetWMH$ID)%>%
    dplyr::select(ID, A1, AX, BETA, SE)%>%
    dplyr::rename(A1_MD=A1, AX_MD=AX, BETA_MD=BETA, SE_MD=SE)
  
  subsetAS=AS_sumstats%>%
    filter(MarkerName %in% subsetWMH$ID)%>%
    dplyr::select(MarkerName, Allele1, Allele2, Effect, StdErr)%>%
    mutate(Allele1=toupper(Allele1), Allele2=toupper(Allele2))%>%
    dplyr::rename(ID=MarkerName, A1_AS=Allele1, AX_AS=Allele2, BETA_AS=Effect, SE_AS=StdErr)
  
  subsetAIS=AIS_sumstats%>%
    filter(MarkerName %in% subsetWMH$ID)%>%
    dplyr::select(MarkerName, Allele1, Allele2, Effect, StdErr)%>%
    mutate(Allele1=toupper(Allele1), Allele2=toupper(Allele2))%>%
    dplyr::rename(ID=MarkerName, A1_AIS=Allele1, AX_AIS=Allele2, BETA_AIS=Effect, SE_AIS=StdErr)
  
  subsetSVS=SVS_sumstats%>%
    filter(MarkerName %in% subsetWMH$ID)%>%
    dplyr::select(MarkerName, Allele1, Allele2, Effect, StdErr)%>%
    mutate(Allele1=toupper(Allele1), Allele2=toupper(Allele2))%>%
    dplyr::rename(ID=MarkerName, A1_SVS=Allele1, AX_SVS=Allele2, BETA_SVS=Effect, SE_SVS=StdErr)
  
  subsetCES=CES_sumstats%>%
    filter(MarkerName %in% subsetWMH$ID)%>%
    dplyr::select(MarkerName, Allele1, Allele2, Effect, StdErr)%>%
    mutate(Allele1=toupper(Allele1), Allele2=toupper(Allele2))%>%
    dplyr::rename(ID=MarkerName, A1_CES=Allele1, AX_CES=Allele2, BETA_CES=Effect, SE_CES=StdErr)
  
  subsetLAS=LAS_sumstats%>%
    filter(MarkerName %in% subsetWMH$ID)%>%
    dplyr::select(MarkerName, Allele1, Allele2, Effect, StdErr)%>%
    mutate(Allele1=toupper(Allele1), Allele2=toupper(Allele2))%>%
    dplyr::rename(ID=MarkerName, A1_LAS=Allele1, AX_LAS=Allele2, BETA_LAS=Effect, SE_LAS=StdErr)
  
  mergeddatasets=subsetWMH%>%
    inner_join(subsetFA, by=c("ID"))%>%
    inner_join(subsetMD,  by=c("ID"))%>%
    inner_join(subsetAS, by=c("ID"))%>%
    inner_join(subsetAIS, by=c("ID"))%>%
    inner_join(subsetLAS, by=c("ID"))%>%
    inner_join(subsetCES, by=c("ID"))%>%
    inner_join(subsetSVS, by=c("ID"))
  
  mergeddatasets2=mergeddatasets
  for(i in 1:nrow(mergeddatasets2)){
    if(mergeddatasets2[i,"A1_WMH"]==mergeddatasets2[i,"AX_FA"]){
      mergeddatasets2[i,"A1_FA"]=mergeddatasets2[i,"A1_WMH"]
      mergeddatasets2[i,"AX_FA"]=mergeddatasets2[i,"AX_WMH"]
      mergeddatasets2[i,"BETA_FA"]=-mergeddatasets2[i,"BETA_FA"]      
    }
    if(mergeddatasets2[i,"A1_WMH"]==mergeddatasets2[i,"AX_MD"]){
      mergeddatasets2[i,"A1_MD"]=mergeddatasets2[i,"A1_WMH"]
      mergeddatasets2[i,"AX_MD"]=mergeddatasets2[i,"AX_WMH"]
      mergeddatasets2[i,"BETA_MD"]=-mergeddatasets2[i,"BETA_MD"]      
    }
    if(mergeddatasets2[i,"A1_WMH"]==mergeddatasets2[i,"AX_AS"]){
      mergeddatasets2[i,"A1_AS"]=mergeddatasets2[i,"A1_WMH"]
      mergeddatasets2[i,"AX_AS"]=mergeddatasets2[i,"AX_WMH"]
      mergeddatasets2[i,"BETA_AS"]=-mergeddatasets2[i,"BETA_AS"]      
    }
    if(mergeddatasets2[i,"A1_WMH"]==mergeddatasets2[i,"AX_AIS"]){
      mergeddatasets2[i,"A1_AIS"]=mergeddatasets2[i,"A1_WMH"]
      mergeddatasets2[i,"AX_AIS"]=mergeddatasets2[i,"AX_WMH"]
      mergeddatasets2[i,"BETA_AIS"]=-mergeddatasets2[i,"BETA_AIS"]      
    }
    if(mergeddatasets2[i,"A1_WMH"]==mergeddatasets2[i,"AX_LAS"]){
      mergeddatasets2[i,"A1_LAS"]=mergeddatasets2[i,"A1_WMH"]
      mergeddatasets2[i,"AX_LAS"]=mergeddatasets2[i,"AX_WMH"]
      mergeddatasets2[i,"BETA_LAS"]=-mergeddatasets2[i,"BETA_LAS"]      
    }
    if(mergeddatasets2[i,"A1_WMH"]==mergeddatasets2[i,"AX_CES"]){
      mergeddatasets2[i,"A1_CES"]=mergeddatasets2[i,"A1_WMH"]
      mergeddatasets2[i,"AX_CES"]=mergeddatasets2[i,"AX_WMH"]
      mergeddatasets2[i,"BETA_CES"]=-mergeddatasets2[i,"BETA_CES"]      
    }
    if(mergeddatasets2[i,"A1_WMH"]==mergeddatasets2[i,"AX_CES"]){
      mergeddatasets2[i,"A1_CES"]=mergeddatasets2[i,"A1_WMH"]
      mergeddatasets2[i,"AX_CES"]=mergeddatasets2[i,"AX_WMH"]
      mergeddatasets2[i,"BETA_CES"]=-mergeddatasets2[i,"BETA_CES"]      
    }
  }
  
  mergeddatasets2=mergeddatasets2%>%
    filter(A1_WMH==A1_FA & A1_WMH==A1_MD & A1_WMH==A1_AS & A1_WMH==A1_AIS & A1_WMH==A1_LAS & A1_WMH==A1_CES & A1_WMH==A1_SVS)%>%
    dplyr::group_by(ID)%>%
    dplyr::mutate(duplicate=sum(!is.na(ID))) %>%
    filter(duplicate==1)
  
  # Regression coefficients and standard errors from ten GWAS studies 
  betas <- mergeddatasets2%>%
    dplyr::select(ID, BETA_WMH, BETA_FA, BETA_MD, BETA_AS, BETA_AIS, BETA_LAS, BETA_CES, BETA_SVS)%>%
    tibble::column_to_rownames("ID")%>%
    dplyr::rename(WMH=BETA_WMH, FA=BETA_FA, MD=BETA_MD, AS=BETA_AS, AIS=BETA_AIS, LAS=BETA_LAS, CES=BETA_CES, SVS=BETA_SVS)
  betas=as.matrix(betas)
  
  ses <-  mergeddatasets2%>%
    dplyr::select(ID, SE_WMH, SE_FA, SE_MD, SE_AS, SE_AIS, SE_LAS, SE_CES, SE_SVS)%>%
    tibble::column_to_rownames("ID")%>%
    dplyr::rename(WMH=SE_WMH, FA=SE_FA, MD=SE_MD, AS=SE_AS, AIS=SE_AIS, LAS=SE_LAS, CES=SE_CES, SVS=SE_SVS)
  ses=as.matrix(ses)
  
  # Trait names and SNP IDs
  traits <- colnames(betas)
  rsid <- rownames(betas)
  
  # Colocalization analysis
  results[[sentinel]]=hyprcoloc(betas[,1:3], ses[,1:3], trait.names=traits[1:3], snp.id=rsid)
}

results2=lapply(results, function(x) as.data.frame(x$results))

library (plyr)
df <- ldply (results2, data.frame)
df=df%>%
  filter(traits!="None")

data2$P_WMH=as.numeric(gsub(data2$P_WMH, replacement = "e", pattern="×10"))
data2$P_FA=as.numeric(gsub(data2$P_FA, replacement = "e", pattern="×10"))
data2$P_MD=as.numeric(gsub(data2$P_MD, replacement = "e", pattern="×10"))
data2=data2%>%
  filter(!is.na(Chr))%>%
  left_join(df, by=c("Lead.SNP"=".id"))
write.csv(data2, "~/tables/results_tophits_20_06_2019_hyprcoloc.csv")

data2%>%
  filter(!is.na(traits) & (grepl("FA", traits) | grepl("WMH", traits) | grepl("MD", traits)))%>%
  dplyr::select(Lead.SNP)%>%
  distinct()

