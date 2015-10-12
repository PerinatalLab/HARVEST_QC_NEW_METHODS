#!/usr/bin/Rscript

# this script is designed for HARVEST genotyping QC
# purpose: compare genotype frequencies in a study data with 
# European genotype frequencies from 1000 Genomes Project
#  by Jonas Bacelis. 2015 Oct 11-12

chr = 22

# define the locations and names
plink = "/home/jonasbac/results/moba24-reference-script_5010jb/plink"

# inputs
study_data_dir = "/home/jonasbac/results/moba24-reference-script_5010jb/data/"
study_data_root = "merge-qc"
study_data_fil = paste(study_data_dir,study_data_root,sep="")
study_dupl_ind = paste(study_data_dir,"duplicated_individuals.txt",sep="")

# outputs
working_dir = "/home/jonasbac/results/moba24-reference-script_5010jb/data/WORK/"
study_data_dpl = paste(working_dir,study_data_root,"_duplicates_chr",chr,sep="")
study_data_oth = paste(working_dir,study_data_root,"_all-others_chr",chr,sep="")
temp_file_rix = paste(working_dir,"tempFile_whichRowsToKeep.txt",sep="")
temp_file_hap = paste(working_dir,"tmpFile_extracted_1kGhap.txt",sep="")
reff_freqs =  paste(working_dir,"1kGreference_gntpFrequencies_chr",chr,".txt",sep="")
        
# read-ins
study_data_dpl_map = paste(working_dir,study_data_root,"_duplicates_chr",chr,".map",sep="")
study_data_dpl_ped = paste(working_dir,study_data_root,"_duplicates_chr",chr,".ped",sep="")
study_data_oth_map = paste(working_dir,study_data_root,"_all-others_chr",chr,".map",sep="")
study_data_oth_ped = paste(working_dir,study_data_root,"_all-others_chr",chr,".ped",sep="")

# 1000 genomes data location and file names
reff_data_dir = "/home/jonasbac/1000G/1000GP_Phase3/"
reff_data_legend = paste(reff_data_dir,"1000GP_Phase3_chr",chr,".legend",sep="") # without .gz suffix
reff_data_sample = paste(reff_data_dir,"1000GP_Phase3.sample",sep="")
reff_data_haplot = paste(reff_data_dir,"1000GP_Phase3_chr",chr,".hap",sep="") # without .gz suffix


########### START

# extract relevant study data
cmnd1 = paste(plink, "--bfile",study_data_fil,"--chr",chr,"--keep",study_dupl_ind,"--recode12 --out",study_data_dpl,sep=" ")
cmnd2 = paste(plink, "--bfile",study_data_fil,"--chr",chr,"--remove",study_dupl_ind,"--recode12 --out",study_data_oth,sep=" ")
system(cmnd1,intern = F)
system(cmnd2,intern = F)



# load study marker list for genotyping study (one particular chromosome)
map = read.table(study_data_dpl_map ,stringsAsFactors = F) #dim(map); head(map)


# decompress needed reference files
cmnd3 = paste("gzip -d ",reff_data_legend,".gz",sep="")
cmnd4 = paste("gzip -d ",reff_data_haplot,".gz",sep="")
system(cmnd3,intern = F)
system(cmnd4,intern = F)


# load 1000G marker list and info (only specific columns)
colcl_leg = rep("NULL",11); colcl_leg[c(1,2,5,9)]=NA  # in order to increase the speed
leg = read.table(reff_data_legend,colClasses=colcl_leg,stringsAsFactors = F,h=T)
leg = leg[grep("biallelic",leg$TYPE,ignore.case=T),]
head(leg); dim(leg)


# check whether coordinates do match
#m = merge(map,leg,by.x="V2",by.y="id",all=F)
#dim(m); m[1:10,]
#sum(m$V4==m$position)  # .. they do

# match markers based on their position (better overlap)
dim(map)
sum(map$V4 %in% leg$position)
sum(duplicated(leg$position))
m = merge(map,leg,by.x="V4",by.y="position",all=F)
dim(m); m[1:10,]
sum(duplicated(m$V2)) #  are ther eany duplicated SNPs| names after merging ? 
sum(duplicated(m$id)) #  are ther eany duplicated SNPs| names after merging ? 
# cleanup duplicates
m = m[which(!duplicated(m$id)),]
dim(m)

# reduce row-wise the very large haplotype file
row_indexes = which(leg$id %in% m$id)
write.table(row_indexes,temp_file_rix,row.names=F,col.names=F,quote=F,sep="\t")
#cmnd5 = paste("awk 'NR==FNR{a[$0]=1;next}a[FNR]' ",temp_file_rix," ",reff_data_haplot," > ",temp_file_hap,sep="") #
cmnd5 = paste("awk 'NR==FNR{a[$0]=1;next} FNR in a' ",temp_file_rix," ",reff_data_haplot," > ",temp_file_hap,sep="") # by Julius
system(cmnd5,intern = F)


# load file describing populations and select only European individuals
smp = read.table(reff_data_sample,stringsAsFactors = F,h=T)
dim(smp); head(smp); table(smp$POP)
# select population
ix = which(smp$GROUP=="EUR")
# determine which columns should be read
cix1 = 1+(ix-1)*2  # column indexes
cix2 = 2+(ix-1)*2
cixs = c(cix1,cix2)
colcl_hap = rep("NULL",nrow(smp)*2)
colcl_hap[cixs] = "integer"

# read haplotypes
hap = read.table(temp_file_hap,colClasses=colcl_hap,stringsAsFactors = F,h=F)
leg = leg[row_indexes,]
dim(hap); hap[1:10,1:10]
dim(leg); head(leg)

##############

n_inds = ncol(hap)/2


cixs1 = 1+(seq(n_inds)-1)*2
cixs2 = 2+(seq(n_inds)-1)*2
s1 = hap[,cixs1]
s2 = hap[,cixs2]

df = s1+s2+1

KGgntpFrq = matrix(NA,nr=nrow(df),nc=3)
for ( i in 1:nrow(df)) {
        tmp = as.numeric(table(factor(df[i,],levels=c(1,2,3))))
        KGgntpFrq[i,] = tmp ; rm(tmp)
}
colnames(KGgntpFrq)=c("AA","AB","BB")
head(KGgntpFrq)


# enforce that minor 
if (sum(KGgntpFrq[,1] > KGgntpFrq[,3])==nrow(KGgntpFrq)) {
# inversion of genotypes counts (to make it comparable with PLINK's output)
# rix_inv = which(leg$EUR<0.5) # rows that need genotype inversion (only for March 2012 version)
tmp = KGgntpFrq
KGgntpFrq[,1] = tmp[,3]
KGgntpFrq[,3] = tmp[,1]
rm(tmp)
}

print(paste("AA>BB ",sum(KGgntpFrq[,1]>KGgntpFrq[,3]),sep=""))
print(paste("BB>AA ",sum(KGgntpFrq[,1]<KGgntpFrq[,3]),sep=""))
print(paste("BB>AB ",sum(KGgntpFrq[,2]<KGgntpFrq[,3]),sep=""))
print(paste("AA<AB ",sum(KGgntpFrq[,1]<KGgntpFrq[,2]),sep=""))

# combine marker info with genotype counts
reff  = data.frame(leg,KGgntpFrq)
head(reff); dim(reff)

write.table(reff,reff_freqs,row.names=F,col.names=T,quote=F,sep="\t")
