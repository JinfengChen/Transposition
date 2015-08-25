pdf("mPing_class_proportion.pdf", height=7, width=4)
x <- read.table("RIL275_RelocaTEi.CombinedGFF.characterized.clean.class.summary", row.names=1, header=T)
x[,2] <- x[,2] + x[,3]
x[,3] <- NULL
x <- t(x)
data <- cbind(x[,1]/sum(x[,1]), x[,2]/sum(x[,2]), x[,3]/sum(x[,3]))
colnames(data) <- colnames(x)
barplot(data, ylab="Proportion", col=c("cornflowerblue", "bisque4"))
legend(1, 1.15,c("Homozygous", "Heterozygous"),bty="n", xpd=TRUE, border="NA",lty=c(0,0),cex=1,fill=c("cornflowerblue", "bisque4"))
dev.off()

