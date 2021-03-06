pdf("mPing_type_proportion.pdf", height=7, width=3)
par(mar=c(5,5,4,5))
x <- read.table("RIL230_RelocaTEi.CombinedGFF.characterized.clean.type.summary", row.names=1)
x <- t(x)
data <- rbind(x[,1]/sum(x[1,]), x[,2]/sum(x[1,]), x[,3]/sum(x[1,]))
rownames(data) <- colnames(x)
#no data labels
barplot(data, ylab="Proportion", space=0.1, col=c("darkseagreen", "lightpink3", "mediumorchid4"))
legend(-0.1, 1.15,c("Parental", "Shared", "Unique"), bty="n", xpd=TRUE, border="NA",lty=c(0,0),cex=1,fill=c("mediumorchid4", "lightpink3", "darkseagreen"))
#labels percentage: paste(round(data[1,1]*100), "%", sep="")
barplot(data, ylab="Proportion", space=0.1, col=c("darkseagreen", "lightpink3", "mediumorchid4"))
legend(-0.1, 1.15,c("Parental", "Shared", "Unique"), bty="n", xpd=TRUE, border="NA",lty=c(0,0),cex=1,fill=c("mediumorchid4", "lightpink3", "darkseagreen"))
text(0.6, 0.4, cex=0.8, paste(round(data[1,1]*100), "%", sep=""))
text(0.6, 0.9, cex=0.8, paste(round(data[2,1]*100), "%", sep=""))
text(1.4, 0.98, cex=0.8, paste(round(data[3,1]*100), "%", sep=""), xpd=TRUE)

#labels percentage and number: paste(x[1], "(", round(data[1,1]*100), "%)", sep="")
par(mar=c(5,5,4,5))
barplot(data, ylab="Proportion", space=0.1, col=c("darkseagreen", "lightpink3", "mediumorchid4"), cex.axis=1.2, cex.names=1.2, cex.lab=1.2)
legend(-0.1, 1.17,c("Parental", "Shared", "Unique"), bty="n", xpd=TRUE, border="NA",lty=c(0,0), cex=1.2, fill=c("mediumorchid4", "lightpink3", "darkseagreen"))
text(0.6, 0.355, cex=1.4, paste("(", x[1], ")", sep=""), col='white')
text(0.6, 0.4, cex=1.4, paste(round(data[1,1]*100), "%", sep=""), xpd=TRUE, col='white')
text(0.6, 0.895, cex=1.4, paste("(", x[2], ")", sep=""), col='white')
text(0.6, 0.94, cex=1.4, paste(round(data[2,1]*100), "%", sep=""), xpd=TRUE, col='white')
text(1.5, 0.96, cex=1.4, paste("(", x[3], ")", sep=""), col='black', xpd=TRUE)
text(1.5, 1, cex=1.4, paste(round(data[3,1]*100), "%", sep=""), xpd=TRUE, col='black')
text(0.3, -0.08, cex=1.2, 'Non-Ref', xpd=TRUE)
text(1.11, -0.082, cex=1.2, 'mPing', font=3, xpd=TRUE)
text(0.6, -0.13, cex=1.2, '(n = 14303)', xpd=TRUE)
dev.off()

