#!/usr/bin/Rscript

# by: jonas b.
# 2015 Nov 12

#####  three strategies on how to utilize triplicates when checking for sample concordance
#####  in HARVEST data (especially Tromso samples)

# the problems/questions: 
# is the estimator biased if all three edges/connections are used per each triplicate
# is the confidence of this estimator biased
# is there anything to gain when estimating three edges instead of two


est_err_normal = sds_err_normal = NULL  # esimated genotyping error prob. when pure duplicates are used (as in MoBa)
est_err_jonas = sds_err_jonas = NULL   # esimated genotyping error prob. when 2 edges from triplicate are useed
est_err_reidar =  sds_err_reidar = NULL  # esimated genotyping error prob. when 3 edges from triplicate are useed

p_errs = seq(0.01,0.5,0.01)  #  probability of an error in the genotype calling

for (p_err in p_errs) {
print(p_err)
# SCENARIO 1: using only true duplicates (as in MoBa)
e_err  = NULL # estimated error 
n_pairs = 1e4  #  number of duplicates (independent pairs)
for (i in 1:1000) { # number of simulations
ind1 = sample(c(1,2),n_pairs,prob = c(1-p_err,p_err),replace = T) # first individuals from pairs
ind2 = sample(c(1,2),n_pairs,prob = c(1-p_err,p_err),replace = T) # second individuals from pairs
errs = sum(ind1!=ind2) # number of mismatches (errors)
total = length(ind1)   # maximal possible number of mismatches
prop_err = errs/total  # estimated error rate per pair
prop_err = prop_err/2  #  estimated error rate per individual
e_err = c(e_err, prop_err) # accumulate the estimates
}
#hist(e_err,breaks=100,col="grey")
est_err_normal = c(est_err_normal,median(e_err))
sds_err_normal = c(sds_err_normal,sd(e_err))
rm(e_err)

# SCENARIO 2: using Jonas's method (using three individuals and only two edges per triplicate)
e_err = NULL
n_pairs = 1e4
for (i in 1:1000) { # number of simulations
        ind1 = sample(c(1,2),n_pairs,prob = c(1-p_err,p_err),replace = T) # first individuals from triplicate
        ind2 = sample(c(1,2),n_pairs,prob = c(1-p_err,p_err),replace = T) # seecond individuals from triplicate
        ind3 = sample(c(1,2),n_pairs,prob = c(1-p_err,p_err),replace = T) # third individuals from triplicate
        err1 = sum(ind1!=ind2) # number of errors in the first edge
        err2 = sum(ind2!=ind3) # number of errors in the second edge
        errs = sum(err1,err2)  # total number of errors
        total = length(ind1)*2 # maximal number of possible errors
        prop_err = errs/total  # estimated error rate per edge (pair)
        prop_err = prop_err/2  # estimated error rate per individual
        e_err = c(e_err, prop_err)
}
est_err_jonas = c(est_err_jonas,median(e_err))
sds_err_jonas = c(sds_err_jonas,sd(e_err))
rm(e_err)

# SCENARIO 3: using Reidar's method (using three individuals and all three edges per triplicate)
e_err = NULL
n_pairs = 1e4 
for (i in 1:1000) {  # number of simulations
        ind1 = sample(c(1,2),n_pairs,prob = c(1-p_err,p_err),replace = T) # first individuals from triplicate
        ind2 = sample(c(1,2),n_pairs,prob = c(1-p_err,p_err),replace = T) # second individuals from triplicate
        ind3 = sample(c(1,2),n_pairs,prob = c(1-p_err,p_err),replace = T) # third individuals from triplicate
        err1 = sum(ind1!=ind2) # number of errors in the first edge (1-2)
        err2 = sum(ind2!=ind3) # number of errors in the second edge (2-3)
        err3 = sum(ind1!=ind3) # number of errors in the third edge (3-1)
        errs = sum(c(err1,err2,err3)) # total number of errors
        total = length(ind1)*3 # maximal number of errors
        prop_err = errs/total  # estimated error rate per edge (pair)
        prop_err = prop_err/2  # estimated error rate per individual
        e_err = c(e_err, prop_err)
}
est_err_reidar = c(est_err_reidar,median(e_err))
sds_err_reidar =  c(sds_err_reidar,sd(e_err))
rm(e_err)

}  # end of cycling through various genotype-call error rates

ylim = max(c(est_err_normal,est_err_jonas,est_err_reidar))
plot(p_errs,est_err_normal,type="l",col="black",ylim = c(0,ylim),
     main = "estimated and true error rate",
     xlab = "the true genotyping error rate for a sample",
     ylab = "estimated genotyping error rate based on duplicate concordance")
abline(0,1,lty=2,col="grey")
points(p_errs,est_err_jonas,type="l",col="red")
points(p_errs,est_err_reidar,type="l",col="blue")

ylim = max(c(sds_err_normal,sds_err_jonas,sds_err_reidar))
plot(p_errs,sds_err_normal,type="l",col="black",ylim = c(0,ylim),
     main = "simulation",
     xlab = "the true genotyping error rate for a sample",
     ylab = "standard deviation of estimator (concordance)")
points(p_errs,sds_err_jonas,type="l",col="red")
points(p_errs,sds_err_reidar,type="l",col="blue")
legend(0.27,0.0008,legend = c("duplicates, 1 edge","triplicates, 2 edges","triplicates, 3 edges"),
       lwd = 1, col = c("black","red","blue"),cex=0.9)



