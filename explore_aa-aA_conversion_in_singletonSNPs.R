#!/usr/bin/Rscript

# this script is designed for HARVEST genotyping QC to test Stefan's idea:
# does the caller (zcall) have a tendency to call heterozygote as minor homozygote
# when singleton SNPs (per whole batch) are investigated. Script test the ratio between 
# AA and AB (where A is minor allele). We expect to find AA > AB (in case of such bias).
#
#  by Jonas Bacelis. 2015 Oct 15

### define the locations and name of the PLINK program
plink = "/home/jonasbac/results/moba24-reference-script_5010jb/plink"

# the study data directory
study_data_dir = "/home/jonasbac/results/moba24-reference-script_5010jb/data/"
# results directory
working_dir = "/home/jonasbac/results/moba24-reference-script_5010jb/data/zcallBias_AAAB/"

# the study data file which is analysed
datasets = c("gencall-raw", "gencall-qc", "zcall-raw", "zcall-qc", "merge-raw","merge-qc")
for (study_data_root in datasets) {
        print(study_data_root)
        
        # temporary files
        temp_file_phenotypes =  paste(working_dir,"tempFile_randomPhenotps_ALLsamples.txt",sep="")
        temp_file_gentypefrq =  paste(working_dir,"tempFile_genotypeCounts_ALLsamples",sep="")
        # output file
        study_ALL_freqs  = paste(working_dir,study_data_root,"_Stefan_AAAB-zcallBias_genotypeFrequencies",sep="")
        
        # generate a random phenotype (necessary to run plink command)
        study_data_fil = paste(study_data_dir,study_data_root,sep="")
        fam_all = read.table(paste(study_data_fil,".fam",sep=""),h=F,stringsAsFactors = F)
        all_phe = rnorm(nrow(fam_all))
        all_phe_df = data.frame(fam_all[,c(1,2)],rPHE=all_phe,stringsAsFactors = F)
        colnames(all_phe_df)=c("FID","IID","rPHE")
        write.table(all_phe_df, temp_file_phenotypes ,row.names=F,col.names=T,quote=F,sep="\t")
        
        # run plink
        cmnd0 = paste(plink,"--bfile",study_data_fil,"--no-pheno --pheno",temp_file_phenotypes,
                      "--assoc qt-means --nonfounders --allow-no-sex --out", temp_file_gentypefrq,sep=" ") # in absolutely all individuals
        # allow-no-sex is necessary, otherwise all phenotypes are ignored..
        system(cmnd0,intern = F); system("wait",intern = F)
        
        # extract and reformat the genotype counts from PLINK output
        tmp0 = read.table(paste(temp_file_gentypefrq,".qassoc.means",sep=""),h=T); head(tmp0)
        tmp1 = tmp0[which(tmp0$VALUE=="COUNTS"),] # only relevant rows
        tmp2 = tmp1[,c("SNP","G11","G12","G22")]
        colnames(tmp2) = c("SNP","AA","AB","BB")
        
        # add information about position etc
        bim = read.table(paste(study_data_fil,".bim",sep=""),stringsAsFactors = F,h=F)
        ALLInd_gntpCnts = merge(tmp2,bim,by.x="SNP",by.y="V2",all.x=T)
        
        write.table(ALLInd_gntpCnts, study_ALL_freqs , row.names=F,col.names=T,quote=F,sep="\t")
}


# cleanup
cmnd1 = paste("rm ",temp_file_phenotypes," ",temp_file_gentypefrq,"* ","*log *nosex",sep="")
system(cmnd1,intern = F)



