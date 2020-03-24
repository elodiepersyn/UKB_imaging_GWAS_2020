library(phenoscanner)
library(dplyr)
library(circlize)

#-----------------------
# READ FILES
#-----------------------
# run the other/sentinel_proxy_files.R to get these files
sentinelsWMH=read.table("~/results/results_GWAS_imaging_WMH_15.10.2018/sentinels.txt", header=TRUE)
proxiesWMH=read.table("~/results/results_GWAS_imaging_WMH_15.10.2018/proxies.txt", header=TRUE)
proxiesWMH=proxiesWMH%>%
  filter(r2>=0.8)

sentinelsFA=read.table("~/results/results_GWAS_imaging_FA_15.10.2018/sentinels.txt", header=TRUE)
proxiesFA=read.table("~/results/results_GWAS_imaging_FA_15.10.2018/proxies.txt", header=TRUE)
proxiesFA=proxiesFA%>%
  filter(r2>=0.8)

sentinelsMD=read.table("~/results/results_GWAS_imaging_MD_15.10.2018/sentinels.txt", header=TRUE)
proxiesMD=read.table("~/results/results_GWAS_imaging_MD_15.10.2018/proxies.txt", header=TRUE)
proxiesMD=proxiesMD%>%
  filter(r2>=0.8)

sentinels=bind_rows(sentinelsWMH,sentinelsFA, sentinelsMD)
proxies=bind_rows(proxiesWMH, proxiesFA, proxiesMD)

#-----------------------
# PHENOSCANNER ANNOTATION
#-----------------------
phenoscannerGWAS=list()
for(rsID in unique(sentinels$rsID)){
  sub=proxies[which(proxies$LEAD_rsID==rsID),]
  variants_ids=as.character(c(rsID,as.character(sub$PROXY_rsID)))
  l=length(variants_ids)
  
  results=list()
  count=1
  nsnps=0
  while(nsnps<l){
    variants=variants_ids[((count-1)*100+1):min((count*100),l)]
    res <- phenoscanner(snpquery=variants, pvalue = 5e-8)
    results[[count]]=res
    count=count+1
    nsnps=nsnps+length(variants)
  }
  results2=results[[1]]
  if(count>2){
    for(i in 2:(count-1)){
      results2$snps=bind_rows(results2$snps, results[[i]]$snps)
      results2$results=bind_rows(results2$results, results[[i]]$results)
    }
  }
  phenoscannerGWAS[[rsID]]=results2
}

phenoscannerGWAS_bind=NULL
for(i in 1:length(names(phenoscannerGWAS))){
  if(nrow(phenoscannerGWAS[[i]]$results)!=0){
    phenoscannerGWAS_bind=bind_rows(phenoscannerGWAS_bind, inner_join(phenoscannerGWAS[[i]]$results, phenoscannerGWAS[[i]]$snps)%>%mutate(sentinelsnp=names(phenoscannerGWAS)[[i]]))
  }else{
    phenoscannerGWAS_bind=bind_rows(phenoscannerGWAS_bind, phenoscannerGWAS[[i]]$snps%>%mutate(sentinelsnp=names(phenoscannerGWAS)[[i]]))
  }
}

#-----------------------
# DATA PROCESSING
#-----------------------

# Getting a table with all annotations
pheno=phenoscannerGWAS_bind%>%
  mutate(chr=as.numeric(as.character(chr)), pos_hg19=as.numeric(as.character(pos_hg19)))%>%
  arrange(chr, pos_hg19)%>%
  mutate(hg19_coordinates=factor(hg19_coordinates, levels=unique(hg19_coordinates)))%>%
  mutate(chr=factor(chr, levels=unique(chr)))

# matrix of snp - trait
pheno2=pheno%>%
  dplyr::select(hgnc, trait, chr, study, sentinelsnp)%>%
  distinct()%>%
  filter(hgnc!="-")%>%
  dplyr::mutate(hgnc=factor(hgnc, levels=unique(hgnc)))%>%
  filter(trait!="<NA>")

# adding categories of traits
# the file with categories of traits has been created from PhenoScanner results and category hand annotation
categories=read.table("./traits_gwas.txt", header=TRUE, sep="\t", fill=TRUE)
categories=categories%>%
  mutate(category=plyr::revalue(category, c("Treatment" = "Other")))%>%
  mutate(category=factor(category,levels=unique(category)))

pheno3=left_join(pheno2, categories)%>%
  mutate(category=factor(category, levels=unique(category)))%>%
  distinct()

# matrix of snp - category
mat=table(pheno3%>%
            dplyr::select(sentinelsnp,category)%>%
            distinct())
mat=mat[, apply(mat,2,sum)!=0]
mat=mat[unique(pheno3$sentinelsnp),]
counts=apply(mat,2,sum)
mat=mat[,order(counts)] # ordering by the number of snps per category
mat=cbind(mat[,which(colnames(mat)%in%"Other"),drop=FALSE],mat[,-which(colnames(mat)%in%"Other"), drop=FALSE]) # setting the "other" category as the last one

#-----------------------
# CHORD DIAGRAM
#-----------------------

# colors
palette=c("#62c1c5","#89cf8b","#fdd247","#f7a60b","#3f3f3f","#b79f4e","#0086fd","#6c2b81","#e1b5ff","#cf202f","#fff42f","#549858","#c1007d","#b7b7b7")
grid.col = c(rev(palette),
             rep("grey", nrow(mat)))
names(grid.col)[1:ncol(mat)]=colnames(mat)
names(grid.col)[(ncol(mat)+1):(ncol(mat)+nrow(mat))]=rownames(mat)

# https://jokergoo.github.io/circlize_book/book/advanced-usage-of-chorddiagram.html

size=0.6

# plotting
pdf(file = "~/figures/Phenoscanner.pdf", height = 7, width = 7)
circos.par(gap.after=2)
chordDiagram(t(mat), grid.col = grid.col, annotationTrack = "grid", annotationTrackHeight=0.02, big.gap = 30, transparency = 0.20,scale = FALSE, 
             preAllocateTracks = list(list(track.height = size*max(strwidth(pheno3$hgnc))), # allocate track 1 height for gene names
                                      list(track.height = size*max(strwidth(rownames(mat))))) # allocate track 2 height for snp ids
)

# we go back to the first track and customize sector labels
# adding snp ids
circos.track(rownames(mat),
             ylim = c(0, 1),
             track.index = 2,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter,
                           CELL_META$ylim[1],
                           CELL_META$sector.index, 
                           facing = "clockwise",
                           niceFacing = TRUE,
                           adj = c(0, 0.5),
                           font=1, 
                           cex=size)
             }, bg.border = NA)

# adding phenotypic category ids
circos.track(colnames(mat),ylim = c(0, 1),
             track.index = 2,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter,
                           CELL_META$ylim[1],
                           CELL_META$sector.index, 
                           facing = "clockwise",
                           niceFacing = TRUE,
                           adj = c(0, 0.5),
                           font=1,
                           cex=size)
             }, bg.border = NA)

# adding gene names
circos.track(rownames(mat),
             ylim = c(0, 1),
             track.index = 1,
             panel.fun = function(x, y) {
               sub=pheno3%>%
                 filter(sentinelsnp==CELL_META$sector.index)%>%
                 dplyr::select(hgnc)%>%
                 distinct()
               
               if(CELL_META$sector.index!="rs72934505"){
                 genename=paste(as.character(sub$hgnc),collapse="\n")
               }else{
                 genename=c("CYP20A1\nICA1L\n...")
               }
               circos.text(CELL_META$xcenter, CELL_META$ylim[1]+(diff(CELL_META$ylim)/5), CELL_META$sector.index,
                           facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), font=3,cex=size,
                           labels=genename)
               circos.lines(c(CELL_META$xlim[1], CELL_META$xlim[2]), c(0, 0), col = "black")
             }, bg.border = NA)

circos.clear()

dev.off()
