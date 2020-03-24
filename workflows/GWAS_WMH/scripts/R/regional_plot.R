args<-commandArgs(trailingOnly=TRUE)

print(args)

tiff(paste(args[1],".tiff",sep=""),2*2**10,2*2**10,pointsize=10,res=300,compression="lzw")
par(mai=c(0,0.75,0.5,0.75))
par(mfrow=c(2,1))
par(xpd=NA)->xpd

ld<-read.table(args[3],header=T,strings=F,colClasses=c("numeric","numeric","character","numeric","numeric","character","numeric"))
SNP.name<-ld[1,3]
names(ld)[6]<-"SNP"
data<-read.table(args[1],header=F)
names(data)[1]<-"SNP"
data<-merge(data,ld,by="SNP")
cols<-heat.colors(5)
BP_min<-min(data$V3)
BP_max<-max(data$V3)
range<-BP_max-BP_min
if(range<100000){BP_min<-BP_min-250000;BP_max<-BP_max+250000}

ymax<-max(-log10(data$V9)+1.5)
topsnp<-data[which(data$SNP==SNP.name),4]
topp<-data[which(data$SNP==SNP.name),9]

data_nor2<-subset(data,!(R2>0.2))
data_r2<-subset(data,(R2>0.2))
CHR<-args[2]


#generate plot in two steps; 1/ SNPs not in LD with top SNP, 2/ SNPs in LD with top SNP
plot(frame=F,data_nor2$V3,-log10(data_nor2$V9),bg=ifelse(data_nor2$R2>0.8,cols[1],ifelse(data_nor2$R2>0.6,cols[2],ifelse(data_nor2$R2>0.4,cols[3],ifelse(data_nor2$R2>0.2,cols[4],"grey")))),pch=21, col=ifelse(data_nor2$R2>0.2,"black","grey"),ylab="",xlab="",ylim=c(0,ymax),xlim=c(BP_min,BP_max),xaxt='n',yaxt='n')
par(new=T)
plot(frame=F,data_r2$V3,-log10(data_r2$V9),
	bg=ifelse(data_r2$R2>0.8,cols[1],
	ifelse(data_r2$R2>0.6,cols[2],
	ifelse(data_r2$R2>0.4,cols[3],
	ifelse(data_r2$R2>0.2,cols[4],"grey")))),
	pch=ifelse(data_r2$SNP==SNP.name,23,21), 
	col=ifelse(data_r2$R2>0.2,"black","grey"),
	ylab="-log10(p-value)",xlab=paste("Genomic position on Chromosome ",args[2],sep=""),
	ylim=c(0,ymax),xlim=c(BP_min,BP_max))
if(topp<1e-7){segments(x0=BP_min,x1=BP_max,y0=-log10(5e-8),lty=3,col="seashell3")}#abline(h=-log10(5e-8),lty=3)}	

#add legend
legend(BP_max-75000,ymax,title=expression(paste("r"^"2")),c("> 0.8","> 0.6","> 0.4","> 0.2","< 0.2"),pch=21,pt.bg=c(cols[1:4],"grey"),col=c(rep("black",4),"grey"))

#add top SNP information
text((BP_min+BP_max)/2,ymax-0.5,paste(topsnp,", p=",formatC(topp,2),sep=""))

#par(mai=c(0,0.75,0.75,0.25))
rectarrows <- function(x0,y0,x1,y1,height,length,...) {
  lwd=par("lwd")
  l0=height*(y1-y0)/sqrt((x1-x0)^2+(y1-y0)^2)
  l1=height*(x1-x0)/sqrt((x1-x0)^2+(y1-y0)^2)
  d0=length*(y1-y0)/sqrt((x1-x0)^2+(y1-y0)^2)
  d1=length*(x1-x0)/sqrt((x1-x0)^2+(y1-y0)^2)
  polygon(x=c(x0+l0,x1+l0-d1,x1,x1-l0-d1,x0-l0),y=c(y0-l1,y1-l1-d0,y1,y1+l1-d0,y0+l1),...)
}

genemap<-read.table("/mnt/lustre/users/k1774406/projects/UKB_imaging_analyses/workflows/smartmeta/glist-hg19",header=F,colClasses=c("character","integer","integer","character"))
allgenes<-subset(genemap,V1==CHR & (V2>BP_min | V3<BP_max))
allgenes<-allgenes[order(allgenes$V2,decreasing=T),]

genes <- allgenes[,4]
gene.start<- allgenes[,2] 
gene.end<- allgenes[,3] 
genes.pos<- as.data.frame(cbind(genes, gene.start, gene.end))

genes.pos[,1]=as.character(genes.pos[,1])
genes.pos[,2]=as.numeric(as.character(genes.pos[,2]))
genes.pos[,3]=as.numeric(as.character(genes.pos[,3]))
plot(frame=F,0,0,type="n",xlim=c(BP_min,BP_max),ylim=c(-1,6),yaxt="n",ylab="",xlab="",xaxt='n')


lapply(seq(length(genes)), function(i) {
	if( genes.pos$gene.start[i]<BP_max && genes.pos$gene.end[i]>BP_min){
	rectarrows( max(genes.pos$gene.start[i],BP_min),i%%6,min(genes.pos$gene.end[i],BP_max),i%%6,col="tomato",height=0.1,length=ifelse(genes.pos$gene.end[i]-genes.pos$gene.start[i]>5000,0,0)) # change to length=1000 or length=10000 for large x-axes, as this is the length of the 'arrow' part of the rectangle
	}
	if( ( mean(c(genes.pos$gene.start[i],genes.pos$gene.end[i]) ) > (BP_max-10000) ) && ( genes.pos$gene.start[i] < BP_max ))
		{ text(BP_max+10000,i%%6+0.4,genes.pos$genes[i],cex=0.85,pos=2,font=3) }
	else if( ( mean(c(genes.pos$gene.start[i],genes.pos$gene.end[i]) ) < (BP_min+10000) ) && ( genes.pos$gene.end[i] > BP_min ))
		{ text(BP_min-10000,i%%6+0.4,genes.pos$genes[i],cex=0.85,pos=4,font=3) }
  	else if( genes.pos$gene.end[i] < BP_min ){}
	else if( genes.pos$gene.start[i] > BP_max){}
	else {text(mean(c(genes.pos$gene.start[i],genes.pos$gene.end[i])),i%%6+0.4,genes.pos$genes[i],cex=0.85,font=3)}
})

dev.off()

