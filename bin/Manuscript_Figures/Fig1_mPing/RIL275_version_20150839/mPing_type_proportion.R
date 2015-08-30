pdf("mPing_type_proportion.pdf", height=7, width=2)
x <- read.table("RIL275_RelocaTEi.CombinedGFF.characterized.clean.type.summary", row.names=1)
x <- t(x)
data <- rbind(x[,1]/sum(x[1,]), x[,2]/sum(x[1,]), x[,3]/sum(x[1,]))
rownames(data) <- colnames(x)
barplot(data, ylab="Proportion", space=0.1, col=c("darkseagreen", "lightpink3", "mediumorchid4"))
legend(-0.1, 1.15,c("Parental", "Shared", "Unique"), bty="n", xpd=TRUE, border="NA",lty=c(0,0),cex=1,fill=c("mediumorchid4", "lightpink3", "darkseagreen"))

dev.off()

