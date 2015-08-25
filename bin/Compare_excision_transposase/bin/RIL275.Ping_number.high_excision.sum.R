pdf("RIL275.Ping_number.high_excision.sum.pdf")
par(mar=c(6,4,4,2), cex=1.2)

x <- read.table("RIL275.Ping_number.high_excision.sum.txt")
data <- rbind(x[,2], x[,3], x[,16])
xx <- barplot(data,beside=TRUE, xlab="Number of Ping", ylab="Number of RILs",cex.names=1,cex.axis=1,border=FALSE,axes=FALSE,ylim=c(0, 60),col=c("black", "dim gray", 'gray'))

axis(1,c(0.5,max(xx)+0.5),line=0,labels=c("",""))
axis(2,seq(0, 60, by=20),line=0,labels=seq(0, 60, by=20))
text(xx[1,]+0.2,rep(-4,7), cex=1, offset=2,labels=x[,1],srt=0,xpd=TRUE)
legend('topright', bty='n', border='NA', lty=c(0,0),cex=1 ,c("RILs_ALL", 'RILs_high_excision_mPing', 'RILs_Chr8:24622994'), fill=c("black", "dim gray", 'gray'))
dev.off()

