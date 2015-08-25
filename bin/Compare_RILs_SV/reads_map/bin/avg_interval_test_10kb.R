x <- read.table('mPing_boundary.linked_50Mb_debug2.table_clean.avg.list')
y1 <- x[x[,2]<10000,]
y2 <- x[x[,2]<20000 & x[,2]>10000,]
y3 <- x[x[,2]<30000 & x[,2]>20000,]
y4 <- x[x[,2]<40000 & x[,2]>30000,]
y5 <- x[x[,2]<50000 & x[,2]>40000,]
t.test(y1[,11], y2[,11], alternative='greater')
t.test(y1[,11], y3[,11], alternative='greater')
t.test(y1[,11], y4[,11], alternative='greater')
t.test(y1[,11], y5[,11], alternative='greater')
t.test(y2[,11], y3[,11], alternative='greater')
t.test(y2[,11], y4[,11], alternative='greater')
t.test(y2[,11], y5[,11], alternative='greater')
