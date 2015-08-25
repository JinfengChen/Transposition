x <- read.table('mPing_boundary.linked_50Mb_debug2.table_clean.avg.list')
y1 <- x[x[,2]<100000,]
y2 <- x[x[,2]<200000 & x[,2]>100000,]
y3 <- x[x[,2]<300000 & x[,2]>200000,]
y4 <- x[x[,2]<400000 & x[,2]>300000,]
y5 <- x[x[,2]<500000 & x[,2]>400000,]
t.test(y1[,11], y2[,11])
t.test(y1[,11], y3[,11])
t.test(y1[,11], y4[,11])
t.test(y1[,11], y5[,11])
t.test(y2[,11], y3[,11])
t.test(y2[,11], y4[,11])
t.test(y2[,11], y5[,11])
