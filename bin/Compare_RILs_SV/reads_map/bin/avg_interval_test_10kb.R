x <- read.table('mPing_boundary.linked_50Mb_debug2.table_clean.avg.list')
y1 <- x[x[,2]<10000,]
y2 <- x[x[,2]<20000 & x[,2]>10000,]
y3 <- x[x[,2]<30000 & x[,2]>20000,]
y4 <- x[x[,2]<40000 & x[,2]>30000,]
y5 <- x[x[,2]<50000 & x[,2]>40000,]
y6 <- x[x[,2]<60000 & x[,2]>50000,]
y7 <- x[x[,2]<70000 & x[,2]>60000,]
y8 <- x[x[,2]<80000 & x[,2]>70000,]
y9 <- x[x[,2]<90000 & x[,2]>80000,]
y10 <- x[x[,2]<100000 & x[,2]>90000,]
print('two tailed test: p-value')
wilcox.test(y1[,3], y2[,3])
wilcox.test(y1[,3], y3[,3])
wilcox.test(y1[,3], y4[,3])
wilcox.test(y1[,3], y5[,3])
wilcox.test(y1[,3], y6[,3])
wilcox.test(y1[,3], y7[,3])
wilcox.test(y1[,3], y8[,3])
wilcox.test(y1[,3], y9[,3])
wilcox.test(y1[,3], y10[,3])
print('one tailed test: p-value')
wilcox.test(y1[,3], y2[,3], alternative='greater')
wilcox.test(y1[,3], y3[,3], alternative='greater')
wilcox.test(y1[,3], y4[,3], alternative='greater')
wilcox.test(y1[,3], y5[,3], alternative='greater')
wilcox.test(y1[,3], y6[,3], alternative='greater')
wilcox.test(y1[,3], y7[,3], alternative='greater')
wilcox.test(y1[,3], y8[,3], alternative='greater')
wilcox.test(y1[,3], y9[,3], alternative='greater')
wilcox.test(y1[,3], y10[,3], alternative='greater')


#t.test(y2[,11], y3[,11], alternative='greater')
#t.test(y2[,11], y4[,11], alternative='greater')
#t.test(y2[,11], y5[,11], alternative='greater')
