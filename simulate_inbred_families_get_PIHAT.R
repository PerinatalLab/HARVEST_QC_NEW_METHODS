#!/usr/bin/Rscript

# jonas b.
# 2015 Nov 10

#####   SIMULATE THE PI_HAT VALUE FOR INBRED-FAMILY MEMBERS:
#  scenario for siblings with parents who are siblings

# we will investigate a single family with 6 members: 
# 2 grandparents, 2 their children who inbreed, 2 resulting grandchildren

n_snps = 1e4  # number of markers
ltrs = c("a","b","c","d") # for quadru-allelic marker

# generate data for grandparents who are unrelated
PAT = data.frame(A1=sample(ltrs,n_snps,replace = T),
      A2 = sample(ltrs,n_snps,replace = T),stringsAsFactors=F)
MAT = data.frame(A1=sample(ltrs,n_snps,replace = T),
      A2 = sample(ltrs,n_snps,replace = T),stringsAsFactors=F)

# simulate transmission to father
pat = PAT # only temporary, data will be replaced
ix = sample(c(T,F),n_snps,replace = T)  # decision vector for transmission or non-transmission
pat[ix,1] = PAT[ix,1] # transmitted A1 from father
pat[!ix,1] = PAT[!ix,2] # transmitted A2 from father
ix = sample(c(T,F),n_snps,replace = T)  # decision vector for transmission or non-transmission
pat[ix,2] = MAT[ix,1] # transmitted A1 from mother
pat[!ix,2] = MAT[!ix,2]  # transmitted A2 from mother

# simulate transmission to mother
mat = PAT # only temporary, data will be replaced
ix = sample(c(T,F),n_snps,replace = T)   # decision vector for transmission or non-transmission
mat[ix,1] = PAT[ix,1]
mat[!ix,1] = PAT[!ix,2]
ix = sample(c(T,F),n_snps,replace = T)   # decision vector for transmission or non-transmission
mat[ix,2] = MAT[ix,1]
mat[!ix,2] = MAT[!ix,2]

# simulate transmission to child 1
ch1 = mat # only temporary, data will be replaced
ix = sample(c(T,F),n_snps,replace = T)   # decision vector for transmission or non-transmission
ch1[ix,1] = pat[ix,1]
ch1[!ix,1] = pat[!ix,2]
ix = sample(c(T,F),n_snps,replace = T)   # decision vector for transmission or non-transmission
ch1[ix,2] = mat[ix,1]
ch1[!ix,2] = mat[!ix,2]

# simulate transmission to child 2
ch2 = mat # only temporary, data will be replaced
ix = sample(c(T,F),n_snps,replace = T)  # decision vector for transmission or non-transmission
ch2[ix,1] = pat[ix,1]
ch2[!ix,1] = pat[!ix,2]
ix = sample(c(T,F),n_snps,replace = T)   # decision vector for transmission or non-transmission
ch2[ix,2] = mat[ix,1]
ch2[!ix,2] = mat[!ix,2]

# estimate the instances (markers) with 0,1,2 shared alleles
fun = function(x) sum(ch1[x,] %in% ch2[x,])
ibd = unlist(lapply(1:n_snps,fun)) 

# summarize IBD values for markers
tbl = table(ibd)/length(ibd)

# estiamate PI_HAT (relatedness) value
tbl[2]*0.5+tbl[3]  # PI_HAT = P(IBD=2)+0.5*P(IBD=1)   (proportion IBD)

# = 0.766

