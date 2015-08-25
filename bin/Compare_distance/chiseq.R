#d <- matrix(c(11,92,2,448), ncol=2)
#only 419 mPings are present in RILs with more than 10% frequency, we take these as parential mPing. 51 of 419 are within 100 kb
d <- matrix(c(11,51,2,368), ncol=2)
fisher.test(d)
chisq.test(d)


