args=commandArgs(trailingOnly=T)
input=args[1]
output=args[2]

data=read.table(gzfile(input), header=T, fill=T)
data2=data[!is.na(data$Z),]
write.table(data2, output, quote=F, row.names=F, col.names=T,sep="\t" )
