#!/usr/bin/Rscript

# this script is designed for HARVEST genotyping QC
# purpose: compare genotype frequencies in a study data with 
# European genotype frequencies from 1000 Genomes Project
#  by Jonas Bacelis. 2015 Oct 11-12

# select the chromosome to work on
chr = 1

### define the locations and name of the PLINK program
plink = "/home/jonasbac/results/moba24-reference-script_5010jb/plink"

### inputs
study_data_dir = "/home/jonasbac/results/moba24-reference-script_5010jb/data/"
# the study file which is analysed
study_data_root = "merge-qc"
study_data_fil = paste(study_data_dir,study_data_root,sep="")
# a list of all individuals in duplicated samples (two-column format: FID, IID)
study_dupl_ind = paste(study_data_dir,"duplicated_individuals.txt",sep="")
# what are the duplicate pairs (which individual belongs to which individual)
study_dupl_match = "/media/local-disk/common/gsexport/moba_24v10_n12874inc135regdup9dualdup/SentrixIDs_moba24_135regdup.txt"
# other files for updating info
recode_famids = "/media/local-disk/common/gsexport/recode-famid-total-moba.fam"
recode_gender = "/media/local-disk/common/gsexport/recode-sex-total-moba.fam"
recode_parent = "/media/local-disk/common/gsexport/recode-parents-total-moba.fam"

### outputs
working_dir = "/home/jonasbac/results/moba24-reference-script_5010jb/data/WORK/"
# temporary outputs
study_data_dpl = paste(working_dir,study_data_root,"_duplicates_chr",chr,sep="")
study_data_oth = paste(working_dir,study_data_root,"_all-others_chr",chr,sep="")
temp_file_rix = paste(working_dir,"tempFile_whichRowsToKeep.txt",sep="")
temp_file_hap = paste(working_dir,"tempFile_extracted_1kGhap.txt",sep="")
temp_file_rndPhe =  paste(working_dir,"tempFile_randomPhenotype_nonDuplSmpls.txt",sep="")
temp_file_gntpCnts = paste(working_dir,"tempFile_gntpCounts_nonDuplInds_onlyFounders",sep="")
temp_genet = paste(working_dir,"tempFile_updatingGeneticFile",sep="")
# result outputs
reff_freqs              = paste(working_dir,study_data_root,"_1000Greference_gntpFrequencies_chr",chr,".txt",sep="")
study_dupl_freqs        = paste(working_dir,study_data_root,"_studyDuplSamples_concordances_chr",chr,".txt",sep="")
study_other_freqs       = paste(working_dir,study_data_root,"_studyNonDuplSamples_gntpFrequencies_chr",chr,".txt",sep="")
study_allele_frqs       = paste(working_dir,study_data_root,"_studyNonDuplSamples_alleleFrequencies_chr",chr,sep="")

# read-ins
study_data_dpl_map = paste(working_dir,study_data_root,"_duplicates_chr",chr,".map",sep="")
study_data_dpl_ped = paste(working_dir,study_data_root,"_duplicates_chr",chr,".ped",sep="")
#study_data_oth_map = paste(working_dir,study_data_root,"_all-others_chr",chr,".map",sep="")
#study_data_oth_ped = paste(working_dir,study_data_root,"_all-others_chr",chr,".ped",sep="")

# 1000 genomes data location and file names
reff_data_dir = "/home/jonasbac/1000G/1000GP_Phase3/"
reff_data_legend = paste(reff_data_dir,"1000GP_Phase3_chr",chr,".legend",sep="") # without .gz suffix
reff_data_sample = paste(reff_data_dir,"1000GP_Phase3.sample",sep="")
reff_data_haplot = paste(reff_data_dir,"1000GP_Phase3_chr",chr,".hap",sep="") # without .gz suffix


########### START

# update genetic file and extract relevant study data to a non-binary format
cmnd00d = paste(plink," --bfile ",study_data_fil," --keep "  ,study_dupl_ind," --chr ",chr," --make-bed --out ",temp_genet,"_0dup",sep="")
cmnd00o = paste(plink," --bfile ",study_data_fil," --remove ",study_dupl_ind," --chr ",chr," --make-bed --out ",temp_genet,"_0oth",sep="")
cmnd01d = paste(plink," --bfile ",temp_genet,"_0dup --update-sex ",recode_gender," --make-bed --out ",temp_genet,"_1dup",sep="")
cmnd01o = paste(plink," --bfile ",temp_genet,"_0oth --update-sex ",recode_gender," --make-bed --out ",temp_genet,"_1oth",sep="")
#cmnd02d = paste(plink," --bfile ",temp_genet,"_1dup --update-ids ",recode_famids," --make-bed --out ",temp_genet,"_2dup",sep="")
cmnd02o = paste(plink," --bfile ",temp_genet,"_1oth --update-ids ",recode_famids," --make-bed --out ",temp_genet,"_2oth",sep="")
#cmnd03d = paste(plink," --bfile ",temp_genet,"_2dup --update-parents ",recode_parent," --make-bed --out ",temp_genet,"_3dup",sep="")
#cmnd04d = paste(plink," --bfile ",temp_genet,"_3dup --filter-founders --recode12 --out ",study_data_dpl,sep="")
cmnd03ddd = paste(plink," --bfile ",temp_genet,"_1dup --recode12 --out ",study_data_dpl,sep="")

cmnd03o = paste(plink," --bfile ",temp_genet,"_2oth --update-parents ",recode_parent," --make-bed --out ",study_data_oth,sep="") # this can stay binary
system(cmnd00d,intern = F); system("wait")
system(cmnd00o,intern = F); system("wait")
system(cmnd01d,intern = F); system("wait")
system(cmnd01o,intern = F); system("wait")
#system(cmnd02d,intern = F); system("wait",intern = F)
system(cmnd02o,intern = F); system("wait",intern = F)
#system(cmnd03d,intern = F); system("wait",intern = F)
#system(cmnd04d,intern = F); system("wait",intern = F)
system(cmnd03o,intern = F); system("wait",intern = F)
system(cmnd03ddd,intern = F); system("wait",intern = F)

# estimate summary statistic
cmnd1 = paste(plink," --bfile ",study_data_oth," --filter-founders --freq --out ",study_allele_frqs,sep="")
system(cmnd1,intern = F); system("wait",intern = F)

#### create a genotype count table for every marker in DUPLICATED samples

# read-in marker data
map = read.table(study_data_dpl_map , stringsAsFactors = F)
# read-in the family and genotype data
tmp = read.table(study_data_dpl_ped , stringsAsFactors = F)
# save family information separately
fam = tmp[,1:6]
# convert ped file format into genotype matrix
n_snps = (ncol(tmp)-6)/2 # is also  = nrow(map)
# create a matrix of genotype values coded as 1,2,3 (one cell per individual*SNP). one row = one individual
mtr = matrix(NA,nr=nrow(tmp),nc=n_snps) 
for (snpix in 1:n_snps) {
        a1 = tmp[,6 + 1 + (snpix-1)*2 ]  # allele A
        a2 = tmp[,6 + 2 + (snpix-1)*2 ]  # allele B
        gntp = a1*a2 # these are not doses, just number-converted allele names!
        gntp[which(gntp==4)]=3
        gntp[is.na(gntp)]=NA
        mtr[,snpix] = gntp
        rm(a1,a2,gntp)
}
#dim(mtr); mtr[1:10,1:10]

# reconcile the ped file 
ped = data.frame(fam,mtr,stringsAsFactors = F)
#dim(ped); ped[1:10,1:10]

# read the duplicate-match info. two columns
tmp = read.table(study_dupl_match,stringsAsFactors = F)
tmp1 = data.frame(ID = tmp$V1,sq=seq(nrow(tmp)),stringsAsFactors = F)
tmp2 = data.frame(ID = tmp$V2,sq=seq(nrow(tmp)),stringsAsFactors = F)
dup = rbind(tmp1,tmp2)  # information on what who is in pair with whom is still here
#dim(dup); head(dup)

m0 = merge(dup,ped,by.x="ID",by.y="V1",all=T) # arguable = FALSE
#dim(m0); m0[1:10,1:10]
m1 = m0[which(m0$ID %in% tmp$V1),] # genetic info for the "A" individuals from pairs
m2 = m0[which(m0$ID %in% tmp$V2),] # genetic info for the "B" individuals from pairs
m1 = m1[order(m1$sq),] # preserve the original matching order (for comparability)
m2 = m2[order(m2$sq),]
m1 = m1[,-c(1:7)] # note that m1 has one fam column more than a pure fam (due to sq column)
m2 = m2[,-c(1:7)]

# extract the genotype counts
rez = matrix(NA,nr=n_snps,nc=9)
for (j in 1:ncol(m1)) {
        gr1 = factor(m1[,j],levels = c(1,2,3)) # missing genotypes not included
        gr2 = factor(m2[,j],levels = c(1,2,3))
        tbl = table(gr1,gr2)
        rez[j,] = c( as.numeric(tbl)[c(1,5,9)], # diagonal (3 types of matches)
                     sum(as.numeric(tbl)[c(2,4)]), # heteroz = homozAA  (4col)
                     sum(as.numeric(tbl)[c(6,8)]), # heteroz = homozBB  (5col)
                     sum(as.numeric(tbl)[c(3,7)]), # homozAA == homozBB (6col)
                     sum(is.na(gr1)),sum(is.na(gr2)), # missing at one individual from pair
                     sum( (is.na(gr1))&(is.na(gr2)))) # missing match  (9col)
        rm(gr1,gr2,tbl)
}
colnames(rez) = c("AAok","ABok","BBok","err_homHET","err_HOMhet","err_homHOM","mis1","mis2","mis12")
rez = data.frame(map,rez,stringsAsFactors = F)
write.table(rez, study_dupl_freqs, row.names=F, col.names=T, quote=F, sep="\t")


#### create a genotype count table for every marker in NONDUPLICATED (remaining) samples
# only founders are used

# generate a random phenotype (necessary to run plink command)
fam_oth = read.table(paste(study_data_oth,".fam",sep=""),h=F,stringsAsFactors = F)
rnd_phe = rnorm(nrow(fam_oth))
rnd_phe_df = data.frame(fam_oth[,c(1,2)],rPHE=rnd_phe,stringsAsFactors = F)
colnames(rnd_phe_df)=c("FID","IID","rPHE")
write.table(rnd_phe_df, temp_file_rndPhe ,row.names=F,col.names=T,quote=F,sep="\t")

# run plink
cmnd3 = paste(plink,"--bfile",study_data_oth,"--no-pheno --pheno",temp_file_rndPhe,
              "--assoc qt-means --allow-no-sex --filter-founders --out", temp_file_gntpCnts,sep=" ")
                        # allow-no-sex is necessary, otherwise all phenotypes are ignored..
system(cmnd3,intern = F); system("wait",intern = F)

# extract and reformat the genotype counts from PLINK output
tmp0 = read.table(paste(temp_file_gntpCnts,".qassoc.means",sep=""),h=T)
tmp1 = tmp0[which(tmp0$VALUE=="COUNTS"),] # only relevant rows
tmp2 = tmp1[,c("SNP","G11","G12","G22")]
colnames(tmp2) = c("SNP","AA","AB","BB")
bim = read.table(paste(study_data_oth,".bim",sep=""),stringsAsFactors = F,h=F)
othInd_gntpCnts = merge(tmp2,bim,by.x="SNP",by.y="V2",all.x=T)
write.table(othInd_gntpCnts, study_other_freqs, row.names=F,col.names=T,quote=F,sep="\t")


############  continue extracting reference data (1000 Genomes)


# load study marker list for genotyping study (one particular chromosome)
map = read.table(study_data_dpl_map ,stringsAsFactors = F) #dim(map); head(map)


# decompress needed reference files
cmnd4 = paste("gzip -d ",reff_data_legend,".gz",sep="")
cmnd5 = paste("gzip -d ",reff_data_haplot,".gz",sep="")
system(cmnd4,intern = F); system("wait",intern = F)
system(cmnd5,intern = F); system("wait",intern = F)


# load 1000G marker list and info (only specific columns)
colcl_leg = rep("NULL",11); colcl_leg[c(1,2,5,9)]=NA  # in order to increase the speed
leg = read.table(reff_data_legend,colClasses=colcl_leg,stringsAsFactors = F,h=T)
valid_indxs = grep("biallelic",leg$TYPE,ignore.case=T) # to be used later
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
#dim(m)

# reduce row-wise the very large haplotype file
row_indexes = which(leg$id %in% m$id)
row_indexes = row_indexes[which(row_indexes %in% valid_indxs)] # leave only those that are biallelic
write.table(row_indexes,temp_file_rix,row.names=F,col.names=F,quote=F,sep="\t")
#cmnd6 = paste("awk 'NR==FNR{a[$0]=1;next}a[FNR]' ",temp_file_rix," ",reff_data_haplot," > ",temp_file_hap,sep="") #
cmnd6 = paste("awk 'NR==FNR{a[$0]=1;next} FNR in a' ",temp_file_rix," ",reff_data_haplot," > ",temp_file_hap,sep="") # by Julius
system(cmnd6,intern = F); system("wait",intern = F)


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
#dim(hap); hap[1:10,1:10]
#dim(leg); head(leg)

##############

n_inds = ncol(hap)/2


cixs1 = 1+(seq(n_inds)-1)*2
cixs2 = 2+(seq(n_inds)-1)*2
s1 = hap[,cixs1]
s2 = hap[,cixs2]

df = s1+s2+1 # 1 = for compatability reasons (gntps = 1,2,3)

KGgntpFrq = matrix(NA,nr=nrow(df),nc=3)
for ( i in 1:nrow(df)) {
        tmp = as.numeric(table(factor(df[i,],levels=c(1,2,3))))
        KGgntpFrq[i,] = tmp ; rm(tmp)
}
colnames(KGgntpFrq)=c("AA","AB","BB")
#head(KGgntpFrq)


# inversion of genotypes counts (to make it comparable with PLINK's output enforce that minor woud be on the left)
rix_inv = which(leg$EUR<0.5) # rows that need genotype inversion
tmp = KGgntpFrq
KGgntpFrq[rix_inv,1] = tmp[rix_inv,3]
KGgntpFrq[rix_inv,3] = tmp[rix_inv,1]
rm(tmp)

print(paste("AA>BB ",sum(KGgntpFrq[,1]>KGgntpFrq[,3]),sep=""))
print(paste("BB>AA ",sum(KGgntpFrq[,1]<KGgntpFrq[,3]),sep=""))
print(paste("BB>AB ",sum(KGgntpFrq[,2]<KGgntpFrq[,3]),sep=""))
print(paste("AA<AB ",sum(KGgntpFrq[,1]<KGgntpFrq[,2]),sep=""))

# combine marker info with genotype counts
reff  = data.frame(leg,KGgntpFrq)
#head(reff); dim(reff)

write.table(reff,reff_freqs,row.names=F,col.names=T,quote=F,sep="\t")




# cleanup
cmnd7 = paste("rm ",temp_file_rix," ",temp_file_hap," ",temp_file_rndPhe," ",temp_file_gntpCnts,"* ",
              study_data_dpl,"* ",study_data_oth,"* ",temp_genet,"* ",sep="")
system(cmnd7,intern = F); system("wait",intern = F)


# compress 1000G to its original state
cmnd8 = paste("gzip ",reff_data_legend,sep="")
cmnd9 = paste("gzip ",reff_data_haplot,sep="")
system(cmnd8,intern = F); system("wait",intern = F)
system(cmnd9,intern = F); system("wait",intern = F)


