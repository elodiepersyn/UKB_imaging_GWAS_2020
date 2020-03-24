args=commandArgs(trailingOnly=TRUE)

# Load libraries
library(data.table)
library(dplyr)

# Column names
ID=c("MARKERNAME", "ID", "SNPID", "SNP", "RSID")
CHR=c("CHR", "CHROM")
POS=c("POS", "BP")
A1=c("ALLELE1", "A1")
A2=c("ALLELE2", "AX", "A2")
P=c("P", "P-VALUE")

# Read files
sumstats=fread(paste("zcat ", args[1]), header=TRUE)
pval_file=fread(args[2], header=TRUE)

# Select interesting columns in the file
colnames(pval_file)=toupper(colnames(pval_file))
pval_file_id=which(colnames(pval_file)%in%ID)[1]
pval_file_a1=which(colnames(pval_file)%in%A1)[1]
pval_file_a2=which(colnames(pval_file)%in%A2)[1]
pval_file_pos=which(colnames(pval_file)%in%POS)[1]
pval_file_chr=which(colnames(pval_file)%in%CHR)[1]
pval_file_p=which(colnames(pval_file)%in%P)[1]

pval_file_select=pval_file%>%
	select(c(
		pval_file_id,
		pval_file_a1,
		pval_file_a2,
		pval_file_pos,
		pval_file_chr,
		pval_file_p
	))

colnames(pval_file_select)=c("SNP", "A1", "A2", "POS", "CHROM", "P")

pval_file_select=pval_file_select%>%
	mutate(A1=toupper(A1), 
			A2=toupper(A2))

# Merge files
merged_file=inner_join(sumstats, pval_file_select)

merged_file1=merged_file%>%
	rename(SNP2=SNP)%>%
	mutate(SNP=paste(CHROM, POS, A1, A2, sep="_"))%>%
	select(SNP, A1, A2, POS, CHROM, P, SNP2, N)

merged_file2=merged_file%>%
    rename(SNP2=SNP)%>%
    mutate(SNP=paste(CHROM, POS, A2, A1, sep="_"))%>%
    select(SNP, A1, A2, POS, CHROM, P, SNP2, N)

merged_file=bind_rows(merged_file1, merged_file2)

write.table(merged_file, paste(args[1], ".temp", sep=""), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

