x <- read.table("MSU7.SlidingWin.HEG4vsExcision.table")
cor.test(x[,7], x[,14], method="pearson")
#r=0.41 and p-value < 2.2e-16

