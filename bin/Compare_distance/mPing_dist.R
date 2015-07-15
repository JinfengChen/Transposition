pdf("mPing_dist.pdf")
x <- read.table("mPing_dist.txt")
hist(x[,3]/1000, breaks=200, xlim=c(0, 1000), xlab="Distance (kb)", ylab="#mPing", main="")
dev.off()
