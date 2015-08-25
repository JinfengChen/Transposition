pdf("mPing_dist_landrace.pdf")
breaks=c(0,1,2,3,4,5,6,7,8,9,10,50)
#HEG4
x <- read.table("HEG4.mPing.distance")
hist(x[,3]/1000, breaks=200, xlim=c(0, 1000), xlab="Distance (kb)", ylab="#mPing", main="HEG4")
hist(x[,3]/1000000, breaks=breaks, xlim=c(0, 10), xlab="Distance (Mb)", ylab="#mPing", main="")

#EG4
x <- read.table("EG4.mPing.distance") 
hist(x[,3]/1000, breaks=200, xlim=c(0, 1000), xlab="Distance (kb)", ylab="#mPing", main="EG4")
hist(x[,3]/1000000, breaks=breaks, xlim=c(0, 10), xlab="Distance (Mb)", ylab="#mPing", main="")

#A123
x <- read.table("A123.mPing.distance") 
hist(x[,3]/1000, breaks=200, xlim=c(0, 1000), xlab="Distance (kb)", ylab="#mPing", main="A123")
hist(x[,3]/1000000, breaks=breaks, xlim=c(0, 10), xlab="Distance (Mb)", ylab="#mPing", main="")

#A119
x <- read.table("A119.mPing.distance") 
hist(x[,3]/1000, breaks=200, xlim=c(0, 1000), xlab="Distance (kb)", ylab="#mPing", main="A119")
hist(x[,3]/1000000, breaks=breaks, xlim=c(0, 10), xlab="Distance (Mb)", ylab="#mPing", main="")



dev.off()

