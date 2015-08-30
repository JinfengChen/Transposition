pdf("HEG4_mPing_AF0.1_Distance.pdf")
x <- read.table("mPing_dist2.txt")
#all HEG4 mPing distance plot
dist <- x[,3]/1000000
dist1 <- dist[dist<1.2]
breaks <- c(seq(0,1.21,by=0.1))
xhist <- hist(dist1, breaks=breaks, xlim=c(0,1.21), labels=FALSE, col="lightblue", main="", xlab="Distance to closest mPing (Mb)", ylab="mPing insertions")

#within 100 kb
dist <- x[,3]/1000
dist2 <- dist[dist<100]
breaks <- c(seq(0, 101, by=10))
xhist <- hist(dist2, breaks=breaks, xlim=c(0,101), labels=FALSE, col="lightblue", main="", xlab="Distance to closest mPing (kb)", ylab="mPing insertions")

dev.off()