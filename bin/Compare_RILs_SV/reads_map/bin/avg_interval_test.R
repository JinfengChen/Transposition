x <- read.table('mPing_boundary.linked_50Mb_debug2.table_clean.avg.list')
y1 <- x[x[,2]<100000,]
y2 <- x[x[,2]<200000 & x[,2]>100000,]
y3 <- x[x[,2]<300000 & x[,2]>200000,]
y4 <- x[x[,2]<400000 & x[,2]>300000,]
y5 <- x[x[,2]<500000 & x[,2]>400000,]
y6 <- x[x[,2]<600000 & x[,2]>500000,]
y7 <- x[x[,2]<700000 & x[,2]>600000,]
y8 <- x[x[,2]<800000 & x[,2]>700000,]
y9 <- x[x[,2]<900000 & x[,2]>800000,]
#y10 <- x[x[,2]<1000000 & x[,2]>900000,]
#y11 <- x[x[,2]<1100000 & x[,2]>1000000,]
#y12 <- x[x[,2]<1200000 & x[,2]>1100000,]
print('two tailed test: p-value')
wilcox.test(y1[,3], y2[,3])
wilcox.test(y1[,3], y3[,3])
wilcox.test(y1[,3], y4[,3])
wilcox.test(y1[,3], y5[,3])
wilcox.test(y1[,3], y6[,3])
wilcox.test(y1[,3], y7[,3])
wilcox.test(y1[,3], y8[,3])
wilcox.test(y1[,3], y9[,3])
print('one tailed test: p-value')
wilcox.test(y1[,3], y2[,3], alternative="greater")
wilcox.test(y1[,3], y3[,3], alternative="greater")
wilcox.test(y1[,3], y4[,3], alternative="greater")
wilcox.test(y1[,3], y5[,3], alternative="greater")
wilcox.test(y1[,3], y6[,3], alternative="greater")
wilcox.test(y1[,3], y7[,3], alternative="greater")
wilcox.test(y1[,3], y8[,3], alternative="greater")
wilcox.test(y1[,3], y9[,3], alternative="greater")
#wilcox.test(y1[,3], y10[,3])
#wilcox.test(y1[,3], y11[,3])
#wilcox.test(y1[,3], y12[,3])
#t.test(y2[,11], y3[,11])
#t.test(y2[,11], y4[,11])
#t.test(y2[,11], y5[,11])
