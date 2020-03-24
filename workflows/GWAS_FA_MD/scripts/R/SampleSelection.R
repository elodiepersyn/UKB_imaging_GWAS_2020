args=commandArgs(trailingOnly=TRUE)

DIR_DOWNL_DATA=args[1]
DIR_PROCESSED_DATA=args[2]
ROOTNAME=args[3]
TRAIT=args[4]
FIELDS=args[5]

library(ukbtools)
library(data.table)
library(dplyr)
library(yaml)

setwd(paste(DIR_PROCESSED_DATA, "/sub_data/", sep=""))

#---------------------------------------------
# PROCESS TAB DATA
#---------------------------------------------
if(!paste(DIR_DOWNL_DATA,"/",ROOTNAME,"_",TRAIT,".r", sep="")%in% dir(".")){
  system(paste("cp ",DIR_DOWNL_DATA,"/",ROOTNAME,".r ",DIR_PROCESSED_DATA,"/sub_data/",ROOTNAME,"_",TRAIT,".r", sep=""))
}

if(!paste(DIR_DOWNL_DATA,"/",ROOTNAME,"_",TRAIT,".html", sep="")%in% dir(".")){
  system(paste("cp ",DIR_DOWNL_DATA,"/",ROOTNAME,".html ",DIR_PROCESSED_DATA,"/sub_data/",ROOTNAME,"_",TRAIT,".html", sep=""))
}

# Fields
res_fields=ukb_df_field(paste(ROOTNAME, "_",TRAIT,sep=""))
fwrite(res_fields, paste("fields_recoded.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)

# Phenotype data
res_ukb=ukb_df(paste(ROOTNAME, "_",TRAIT,sep=""))
save(res_ukb, file=paste("res_ukb.rda", sep=""))

#---------------------------------------------
# PROCESS PHENOTYPE DATA
#---------------------------------------------
FILTER_SUMMARY=list()


FILTER_SUMMARY[["Initial_number_of_individuals"]]=nrow(res_ukb)

# Selection of unsuspicious individuals

# Brain MRI measurement completed
# grep(colnames(res_ukb), pattern="f12188")
res_ukb$brain_mri_measurement_completeduses_datacoding_21_f12188_2_0=as.factor(res_ukb$brain_mri_measurement_completeduses_datacoding_21_f12188_2_0)
res_ukb$brain_mri_measuring_methoduses_datacoding_470_f12187_2_0=as.factor(res_ukb$brain_mri_measuring_methoduses_datacoding_470_f12187_2_0)
res_ukb$believed_safe_to_perform_brain_mri_scanuses_datacoding_634_f12139_2_0=as.factor(res_ukb$believed_safe_to_perform_brain_mri_scanuses_datacoding_634_f12139_2_0)

summary(res_ukb$brain_mri_measurement_completeduses_datacoding_21_f12188_2_0)
summary(res_ukb$brain_mri_measuring_methoduses_datacoding_470_f12187_2_0)
summary(res_ukb$believed_safe_to_perform_brain_mri_scanuses_datacoding_634_f12139_2_0)

#res_ukb=res_ukb%>%
#		filter(!brain_mri_measurement_completed_f12188_2_0=="no" & brain_mri_measuring_method_f12187_2_0=="Direct entry" & !believed_safe_to_perform_brain_mri_scan_f12139_2_0=="No" & is.na(reason_believed_unsafe_to_perform_brain_mri_f12663_2_0) & is.na(reason_brain_mri_not_completed_f12704_2_0) & is.na(reason_brain_mri_not_performed_f12652_2_0))


# List of imaging fields to process

fields=fread(FIELDS, header=TRUE)
fields_to_keep=fields$col.name


fields_to_remove=list(
  noncancer="f20002_", 
  ICD9=c("f41203_", "f41205_"),
  ICD10=c("f41202_", "f41204_")
)

codes_to_remove=list(
  noncancer=list(
    stroke=c("1081","1086","1491","1583"),
    multiple_sclerosis=c("1261"),
    parkinson=c("1262"),
    dementia=c("1263"),
    neurodegenerative=c("1397")
  ),
  ICD9=list(
    stroke=c("430","431","434","436"),
    multiple_sclerosis=c("340"),
    parkinson=c("332"),
    dementia=c("290"),
    neurodegenerative=c("341")
  ),
  ICD10=list(
    stroke=c("I60","I61","I63","I64"),
    multiple_sclerosis=c("G35"),
    parkinson=c("G20"),
    dementia=c("F00", "F01", "F02", "F03"),
    neurodegenerative=c("G30", "G31", "G32", "G36", "G37")
  )
)

# Selection of individuals
columns_to_keep=unlist(lapply(fields_to_keep, function(x) grep(colnames(res_ukb), pattern=x)))
ind_to_keep=lapply(columns_to_keep, function(x) which(!is.na(res_ukb[,x])))
ind_to_keep= Reduce(intersect, ind_to_keep)
FILTER_SUMMARY[["Number_of_individuals_with_all_imaging_information"]]=length(ind_to_keep)

columns_remove=unlist(lapply(unlist(fields_to_remove), function(x) grep(colnames(res_ukb), pattern=x)))
ind_to_remove=lapply(columns_remove, function(x) which(res_ukb[,x] %in% unlist(codes_to_remove)))
ind_to_remove=Reduce(union, ind_to_remove)

for(i in columns_remove){
  res_ukb[,i]=factor(res_ukb[,i])
}

ind_to_keep2=which(! ind_to_keep %in% ind_to_remove)

FILTER_SUMMARY[["Number_of_individuals_with_at_least_one_code_to_remove"]]=length(ind_to_keep)-length(ind_to_keep2)
FILTER_SUMMARY[["Number_of_individuals_after_removing_codes"]]=length(ind_to_keep2)

# Explore removed samples
removed=res_ukb[intersect(ind_to_remove, ind_to_keep),columns_remove]
for(i in 1:length(columns_remove)){
  removed[,i]=as.factor(removed[,i])
}

replace_codes=function(x){
  x[!x %in% unlist(codes_to_remove)]=NA
  x=factor(x)
  return(x)
}

removed2=apply(removed, 2, function(x) replace_codes(x))
NA_count=apply(removed2, 2, function(x) sum(is.na(x)))
removed3=removed2[,NA_count!=nrow(removed2)]

new_column=function(x){
  sub=strsplit(x, split="_")
  new=paste(sub[[1]][1:(length(sub[[1]])-1)], collapse="_")
  return(new)
}


big_columns=unique(unlist(lapply(colnames(removed3), function(x) new_column(x))))
print(big_columns)
level_counts=list()
for(i in big_columns){
  ind=grep(colnames(removed3), pattern=i)
  tab_from_table=apply(removed3[,ind],2 ,function(x) data.frame(table(x)))
  merge_table=Reduce(function(...) merge(..., by="x", all=T),tab_from_table)
  merge_table2=data.frame(code=merge_table$x, Freq=apply(merge_table[,-1], 1, function(x) sum(x, na.rm=T)))
  level_counts[[i]]=merge_table2
}

total_removed=sum(unlist(lapply(level_counts, function(x) sum(x$Freq))))-sum(apply(removed3, 1, function(x) sum(!is.na(x))-1))
CODE_PER_INDIVIDUAL_SUMMARY=summary(apply(removed3, 1, function(x) sum(!is.na(x))))  # Number of disease per person


# Removing outliers
if (!"outlier"%in% dir(".")){
  dir.create("outlier")
}


list_ind=list()
sub=res_ukb[ind_to_keep[ind_to_keep2],columns_to_keep]
#ks_p=rep(NA, ncol(sub))
for(i in 1:ncol(sub)){
  sub[,i]=as.numeric(sub[,i])
  #ks_p[i]=ks.test(sub[,i], "pnorm", mean(sub[,i]), sd(sub[,i]))$p.value
  scalesub=scale(sub[,i], center = TRUE, scale = TRUE)
  list_ind[[i]]=which(scalesub>=-6 & scalesub<=6)
}

nind_list_ind=unlist(lapply(list_ind, function(x) length(x)))
intersect_list_ind=Reduce(intersect, list_ind)
FILTER_SUMMARY[["Number_of_outliers"]]=length(intersect_list_ind)-length(ind_to_keep2)
FILTER_SUMMARY[["Number_of_individuals_after_removing_outliers"]]=length(intersect_list_ind)

# Exploring description fields
sub2=res_ukb[ind_to_keep[ind_to_keep2][intersect_list_ind],]%>%
  select(brain_mri_measurement_completeduses_datacoding_21_f12188_2_0, brain_mri_measuring_methoduses_datacoding_470_f12187_2_0, believed_safe_to_perform_brain_mri_scanuses_datacoding_634_f12139_2_0)
DESCRIPTION_SUMMARY=summary(sub2)


#---------------------------------------------
# OUTPUTS
#---------------------------------------------
write_yaml(as.yaml(level_counts), "removed_codes_summary.yaml")
write.table(as.matrix(CODE_PER_INDIVIDUAL_SUMMARY), "number_of_codes_in_removed_individuals.txt")
write_yaml(as.yaml(FILTER_SUMMARY), "filter_summary.txt")
write.table(as.matrix(DESCRIPTION_SUMMARY), "description_summary.txt")
ids=matrix(res_ukb[ind_to_keep[ind_to_keep2][intersect_list_ind],]$eid, ncol=1)
write.table(ids, paste(DIR_PROCESSED_DATA, "/samples/", ROOTNAME, "_",TRAIT,".samples.txt",sep=""), sep="\t", row.names=F, col.names=F, quote=F)
