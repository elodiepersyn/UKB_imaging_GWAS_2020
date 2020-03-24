args=commandArgs(trailingOnly=T)


cell_type_dir=args[1]
input=args[2]
output=args[3]

data=read.table(input, header=F)
data$V3=paste(cell_type_dir, data[,2], sep="/")
write.table(data[,c(1,3)],output, quote=F, row.names=F, col.names=F)
