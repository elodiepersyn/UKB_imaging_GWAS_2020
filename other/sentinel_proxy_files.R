library(data.table)
library(dplyr)

# FILES FROM THE GWAS_FA_MD and GWAS_WMH SNAKEMAKE PIPELINES
# for FA
SUMSTATS_FILE="~/results/results_GWAS_imaging_FA_15.10.2018/GWAS_PC1_trait_ALLchr.PC1_trait_info.glm.linear"
CLUMPED_FILE = "~/results/results_GWAS_imaging_FA_15.10.2018/GWAS_PC1_trait_ALLchr.PC1_trait_info.glm.linear.clumped"
files_ld=dir("~/results/results_GWAS_imaging_FA_15.10.2018/Regional_Plots_PC1_trait/", pattern=".ld")
files_ld=paste("~/results/results_GWAS_imaging_FA_15.10.2018/Regional_Plots_PC1_trait/", files_ld, sep="")

# for MD
# SUMSTATS_FILE="~/results/results_GWAS_imaging_MD_15.10.2018/GWAS_PC1_trait_ALLchr.PC1_trait_info.glm.linear"
# CLUMPED_FILE = "~/results/results_GWAS_imaging_MD_15.10.2018/GWAS_PC1_trait_ALLchr.PC1_trait_info.glm.linear.clumped"
# files_ld=dir("~/results/results_GWAS_imaging_MD_15.10.2018/Regional_Plots_PC1_trait/", pattern=".ld")
# files_ld=paste("~/results/results_GWAS_imaging_MD_15.10.2018/Regional_Plots_PC1_trait/", files_ld, sep="")

# for WMH
# SUMSTATS_FILE="~/results/results_GWAS_imaging_WMH_15.10.2018/GWAS_logWMHnorm_ALLchr.logWMHnorm_info.glm.linear"
# CLUMPED_FILE = "~/results/results_GWAS_imaging_WMH_15.10.2018/metaGWAS_logWMHnorm_1.tbl.clumped"
# files_ld=dir("~/results/results_GWAS_imaging_WMH_15.10.2018/Regional_Plots/", pattern=".ld")
# files_ld=paste("~/results/results_GWAS_imaging_WMH_15.10.2018/Regional_Plots/4", files_ld, sep="")

#-------------------------------------------
# reading files
#-------------------------------------------
sumstats=fread(SUMSTATS_FILE, verbose = FALSE, showProgress = FALSE)
clumped=read.table(CLUMPED_FILE, header=TRUE)

ldinfo_tot=NULL
for(file in files_ld){
  ldinfo=fread(file, header=TRUE)
  ldinfo_tot=rbind(ldinfo_tot,ldinfo)
}

#-------------------------------------------
# getting a column CHROM_POS in sumstats file
#-------------------------------------------
sumstats=sumstats%>%
  mutate(CHROM_POS=paste(CHROM,POS, sep="_"))

#-------------------------------------------
# getting the most significant snp per 250 kb window
#-------------------------------------------
# order by position
clumped=clumped%>%
  arrange(CHR, BP)

# get the distance between two consecutive snps
distance_vectors=list()

for (chr in unique(clumped$CHR)){
  sub=clumped%>%filter(CHR==chr)
  ncolsub=nrow(sub)
  distance_vectors[[as.character(chr)]]=abs(sub$BP[-ncolsub]-sub$BP[-1])
}

# grouping by distance

group_by_distance=function(x){
  group=1
  groups=1
  if(length(x)!=0){
    for(i in 1:length(x)){
      if(x[i]>250000){
        group=group+1
      }
      groups=c(groups,group)
    }
  }
  return(groups)
}

groups_by_chr=lapply(distance_vectors, function(x) group_by_distance(x))
clumped$group=unlist(groups_by_chr)

clumped=clumped%>%
  group_by(CHR, group)%>%
  mutate(min=min(P))%>%
  filter(P==min)%>%
  ungroup()

#-------------------------------------------
# merging clumped file with sumstats file to get rsids
#-------------------------------------------
assoc_clumped=inner_join(clumped, sumstats, by=c("SNP"="CHROM_POS"))
assoc_clumped$P=as.numeric(formatC(assoc_clumped$P.x,format="e",digits=2))
assoc_clumped=assoc_clumped%>%
  arrange(CHROM,POS)%>%
  dplyr::select(CHROM, POS, ID, A1, AX, A1_FREQ, BETA, SE, P, SP2)%>%
  mutate(BETA=round(BETA,3), SE=round(SE,3), A1_FREQ=round(A1_FREQ,2))

#-------------------------------------------
# getting the sentinel.txt file
#-------------------------------------------
## tab separated .txt filr containing rsids, chromosomes, and GRCh37 coordinates of sentinel variants

sentinels=assoc_clumped%>%
  dplyr::select(ID,CHROM,POS)%>%
  dplyr::rename(rsID=ID, CHR=CHROM, START=POS)%>%
  mutate(END=START)%>%
  filter(grepl("rs",rsID))
write.table(sentinels,"~/results/results_GWAS_imaging_FA_15.10.2018/PROGEM/sentinels.txt", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)

#-------------------------------------------
# getting the proxies.txt file
#-------------------------------------------
#tab separated .txt filr containing rsids, chromosomes, and GRCh37 coordinates of proxy variants

# starting with getting the clumps
proxies=assoc_clumped%>%
  dplyr::select(ID, CHROM, POS, SP2)%>%
  mutate(CHROM_POS=paste(CHROM,POS, sep="_"))%>%
  as.data.frame()

# getting the proxies IDs (CHROM_POS) from the clumped file
proxieslist=apply(proxies, 1, function(x) gsub(strsplit(x[4], split=",")[[1]], pattern="\\(1\\)", replacement = ""))

# getting the ld between proxies and the sentinel snp from the ld files
proxiesfile=NULL
for(i in 1:length(proxieslist)){
  sub=proxieslist[[i]]
  numrows=length(sub)
  sentinel_proxy=cbind.data.frame(LEAD_rsID=rep(proxies[i,1],numrows), CHROM_POS=rep(proxies[i,5],numrows), PROXY_chrpos=sub)
  submerge=inner_join(sentinel_proxy,ldinfo_tot, by=c("CHROM_POS"="SNP_A", "PROXY_chrpos"="SNP_B"))
  submerge=submerge%>%
    dplyr::select(LEAD_rsID,PROXY_chrpos,CHR_B,BP_B,R2)
  proxiesfile=rbind(proxiesfile,submerge)
}

# getting proxies rsids from the sumstats file
proxiesfile2=inner_join(proxiesfile,
                        sumstats%>%dplyr::select(CHROM_POS,ID)%>%distinct(),
                        by=c("PROXY_chrpos"="CHROM_POS")
)

# final data processig to get the proxies.txt file
proxiesfile3=proxiesfile2%>%
  dplyr::rename(PROXY_rsID=ID,
                PROXY_CHR=CHR_B,
                PROXY_START=BP_B,
                r2=R2)%>%
  mutate(PROXY_END=PROXY_START)%>%
  dplyr::select(PROXY_rsID,PROXY_CHR,PROXY_START,PROXY_END,LEAD_rsID,r2)%>%
  filter(grepl("rs",PROXY_rsID))

write.table(proxiesfile3,"~/results/results_GWAS_imaging_FA_15.10.2018/PROGEM/proxies.txt", sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
