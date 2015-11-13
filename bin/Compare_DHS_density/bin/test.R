
error.bar <- function(x, y, upper, lower=upper, color,length=0.06,...){
     if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
     stop("vectors must be same length")
     arrows(x,y+upper, x, y-lower, col=color,angle=90, code=3, length=length, ...)
}

pdf("test.pdf")

par(mar=c(6,4,4,2), cex=1.2)
dist <- read.table("test.sum")
plot(rev(dist[,2]), type='b', pch= 1,lwd = 2 , col="aquamarine3", xaxt='n', frame.plot = FALSE, ylim=c(0,100), ylab="Proportion", xlab="")
#lines(rev(sim[,2]), type='b',pch= 20, cex=0.2,lwd = 2 , col="dim gray")
error.bar(1:length(dist[,2]), rev(dist[,2]), rev(dist[,3]), rev(dist[,3]), 'dim gray')
axis(1,seq(1:length(dist[,1])),line=0, labels=rep("",length(dist[,1])))
text(seq(1:length(dist[,1])),rep(-0.04,7), cex=1, offset=2,labels=rev(dist[,1]/1000),srt=0,xpd=TRUE)
mtext("Distance to TSD (kb)", side=1, cex=1.2, at=6, line=3)
legend('topright', bty='n', border='NA', lty= c(1,2), pch = c(1,20), cex=1 , lwd = 2 ,col=c("aquamarine3", "dim gray"), c("Unique", "Control"))

dist <- read.table("test_all.sum")
dist <- subset(dist, V1<=2500 & V1>=500)
plot(rev(dist[,3]/1000000), type='l', pch= 1,lwd = 2 , col="aquamarine3", xaxt='n', frame.plot = FALSE, ylim=c(5000/1000000,15000/1000000), ylab="Normalized DNase reads", xlab="")
#lines(rev(sim[,2]), type='b',pch= 20, cex=0.2,lwd = 2 , col="dim gray")
#error.bar(1:length(dist[,2]), rev(dist[,2]), rev(dist[,3]), rev(dist[,3]), 'dim gray')
axis(1,seq(0, 2000, by=200),line=0, labels=rep("",11))
text(seq(0, 2000, by=200), rep(0.004, 11), cex=1, offset=2,labels=seq(-1, 1, by=0.2), srt=0,xpd=TRUE)
mtext("Distance to TSD (kb)", side=1, cex=1.2, at=1000, line=3)
legend('topright', bty='n', border='NA', lty= c(1,2), pch = c(1,20), cex=1 , lwd = 2 ,col=c("aquamarine3", "dim gray"), c("Unique", "Control"))

dev.off()

