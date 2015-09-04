pdf("Excision_distance_HEG4.matrix_events.plot.pdf")
par(mar=c(5,5,5,2))
read.table("Excision_distance_HEG4.matrix_events.1.txt", header=T) -> x
plot(x[,3]/1000, x[,2], type="p", pch=18, col="black", xlab="Distance (kb)", ylab="Number of excision", xlim=c(0, 1200000/1000), ylim=c(0, 20), cex.axis=1.4, cex.lab=1.4)
x <- x[x[,3]<200000,]
plot(x[,3]/1000, x[,2], type="p", pch=18, col="black", xlab="Distance (kb)", ylab="Number of excision", xlim=c(0, 200000/1000), ylim=c(0, 20), cex.axis=1.4, cex.lab=1.4)
#axis(1,at=atx,labels=atx)
#axis(2,at=aty,labels=aty)
dev.off()

