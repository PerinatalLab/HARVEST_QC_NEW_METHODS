#!/usr/bin/Rscript

# jonas b.
# 2015 Nov 10

#####   SIMULATE THE PI_HAT VALUE FOR INBRED-FAMILY MEMBERS:
#  scenario for families with parents who are siblings

# we will investigate a single family with 8 members: 
# 3 grandparents, 3 their children 2 of which are inbreed, 2 resulting grandchildren

n_snps = 1e4  # number of markers
ltrs = c("a","b","c","d") # genetic pool for quadru-allelic markers

############# GENERATION 1

# generate data for grandparents who are unrelated
PAT = data.frame(A1=sample(ltrs,n_snps,replace = T),
      A2 = sample(ltrs,n_snps,replace = T),stringsAsFactors=F) # grandfather
MAT = data.frame(A1=sample(ltrs,n_snps,replace = T),
      A2 = sample(ltrs,n_snps,replace = T),stringsAsFactors=F) # grandmother 1
OTH = data.frame(A1=sample(ltrs,n_snps,replace = T),
      A2 = sample(ltrs,n_snps,replace = T),stringsAsFactors=F) # grandmother 2

############# GENERATION 2

# simulate transmission to inbred father (son of grandfather and grandmother 1)
pat = PAT # only temporary, data will be replaced
ix = sample(c(T,F),n_snps,replace = T)  # decision vector for transmission or non-transmission
pat[ix,1] = PAT[ix,1] # transmitted A1 from father
pat[!ix,1] = PAT[!ix,2] # transmitted A2 from father
ix = sample(c(T,F),n_snps,replace = T)  # decision vector for transmission or non-transmission
pat[ix,2] = MAT[ix,1] # transmitted A1 from mother
pat[!ix,2] = MAT[!ix,2]  # transmitted A2 from mother

# simulate transmission to inbred mother (daughter of grandfather and grandmother 1)
mat = PAT # only temporary, data will be replaced
ix = sample(c(T,F),n_snps,replace = T)   # decision vector for transmission or non-transmission
mat[ix,1] = PAT[ix,1]
mat[!ix,1] = PAT[!ix,2]
ix = sample(c(T,F),n_snps,replace = T)   # decision vector for transmission or non-transmission
mat[ix,2] = MAT[ix,1]
mat[!ix,2] = MAT[!ix,2]

# # simulate transmission to half-sib mother (daughter of grandfather and grandmother 2)
ix = sample(c(T,F),n_snps,replace = T)
hsm = PAT
hsm[ix,1] = PAT[ix,1]
hsm[!ix,1] = PAT[!ix,2]
ix = sample(c(T,F),n_snps,replace = T)
hsm[ix,2] = OTH[ix,1]
hsm[!ix,2] = OTH[!ix,2]

#############  GENERATION 3

# simulate transmission to child 1 of inbred parents
ch1 = mat # only temporary, data will be replaced
ix = sample(c(T,F),n_snps,replace = T)   # decision vector for transmission or non-transmission
ch1[ix,1] = pat[ix,1]
ch1[!ix,1] = pat[!ix,2]
ix = sample(c(T,F),n_snps,replace = T)   # decision vector for transmission or non-transmission
ch1[ix,2] = mat[ix,1]
ch1[!ix,2] = mat[!ix,2]

# simulate transmission to child 2 of inbred parents
ch2 = mat # only temporary, data will be replaced
ix = sample(c(T,F),n_snps,replace = T)  # decision vector for transmission or non-transmission
ch2[ix,1] = pat[ix,1]
ch2[!ix,1] = pat[!ix,2]
ix = sample(c(T,F),n_snps,replace = T)   # decision vector for transmission or non-transmission
ch2[ix,2] = mat[ix,1]
ch2[!ix,2] = mat[!ix,2]

# simulate transmission to child 1 of half-sib mother 
hsmch1 = pat # only temporary, data will be replaced
ix = sample(c(T,F),n_snps,replace = T)   # decision vector for transmission or non-transmission
hsmch1[ix,1] = pat[ix,1]
hsmch1[!ix,1] = pat[!ix,2]
ix = sample(c(T,F),n_snps,replace = T)   # decision vector for transmission or non-transmission
hsmch1[ix,2] = hsm[ix,1]
hsmch1[!ix,2] = hsm[!ix,2]

# simulate transmission to child 2 of half-sib mother
hsmch2 = pat # only temporary, data will be replaced
ix = sample(c(T,F),n_snps,replace = T)   # decision vector for transmission or non-transmission
hsmch2[ix,1] = pat[ix,1]
hsmch2[!ix,1] = pat[!ix,2]
ix = sample(c(T,F),n_snps,replace = T)   # decision vector for transmission or non-transmission
hsmch2[ix,2] = hsm[ix,1]
hsmch2[!ix,2] = hsm[!ix,2]

##########################################  PI_HAT values

# for two full-sibs of inbred parents
fun = function(x) sum(ch1[x,] %in% ch2[x,])
ibd = unlist(lapply(1:n_snps,fun))
tbl = table(ibd)/length(ibd)
tbl[2]*0.5+tbl[3]
# PI_HAT = 0.76

# for two sibs of half-sib mother
fun = function(x) sum(hsmch1[x,] %in% hsmch2[x,])
ibd = unlist(lapply(1:n_snps,fun))
tbl = table(ibd)/length(ibd)
tbl[2]*0.5+tbl[3]
# PI_HAT = 0.74

# for two inbred half-sibs who are also inbred cousins 
fun = function(x) sum(ch1[x,] %in% hsmch1[x,])
ibd = unlist(lapply(1:n_snps,fun))
tbl = table(ibd)/length(ibd)
tbl[2]*0.5+tbl[3]
# PI_HAT = 0.675


