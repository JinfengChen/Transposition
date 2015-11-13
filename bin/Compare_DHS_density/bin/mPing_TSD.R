
pdf("mPing_TSD.pdf")
par(mar=c(6,4,4,2), cex=1.2)
dist <- read.table('mPing_TSD.profile.sum')
dist <- subset(dist, V1<=2500 & V1>=500)
plot(rev(dist[,3]/1000000), type='l', pch= 1,lwd = 2 , col="aquamarine3", xaxt='n', frame.plot = FALSE, ylim=c(10000/1000000,20000/1000000), ylab="Normalized DNase reads", xlab="")
#lines(rev(sim[,2]), type='b',pch= 20, cex=0.2,lwd = 2 , col="dim gray")
#error.bar(1:length(dist[,2]), rev(dist[,2]), rev(dist[,3]), rev(dist[,3]), 'dim gray')
axis(1,seq(0, 2000, by=200),line=0, labels=rep("",11))
text(seq(0, 2000, by=200), rep(0.009, 11), cex=1, offset=2,labels=seq(-1, 1, by=0.2), srt=0,xpd=TRUE)
mtext("Distance to TSD (kb)", side=1, cex=1.2, at=1000, line=3)
legend('topright', bty='n', border='NA', lty= c(1,2), pch = c(1,20), cex=1 , lwd = 2 ,col=c("aquamarine3", "dim gray"), c("Unique", "Control"))
dev.off()

