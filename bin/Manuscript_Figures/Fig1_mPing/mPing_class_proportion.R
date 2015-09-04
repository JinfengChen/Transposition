pdf("mPing_class_proportion.pdf", height=7, width=4)
x <- read.table("RIL230_RelocaTEi.CombinedGFF.characterized.clean.class.summary", row.names=1, header=T)
x[,2] <- x[,2] + x[,3]
x[,3] <- NULL
x <- t(x)
data <- cbind(x[,1]/sum(x[,1]), x[,2]/sum(x[,2]), x[,3]/sum(x[,3]))
colnames(data) <- colnames(x)
#all three, parental, shared and unique
#no data
barplot(data, ylab="Proportion", col=c("cornflowerblue", "bisque4"))
legend(1, 1.15,c("Homozygous", "Heterozygous"),bty="n", xpd=TRUE, border="NA",lty=c(0,0),cex=1,fill=c("cornflowerblue", "bisque4"))
#percentage
barplot(data, ylab="Proportion", col=c("cornflowerblue", "bisque4"))
legend(1, 1.15,c("Homozygous", "Heterozygous"),bty="n", xpd=TRUE, border="NA",lty=c(0,0),cex=1,fill=c("cornflowerblue", "bisque4"))
text(0.7, 0.4, cex=0.7, paste(round(data[1,1]*100), "%", sep=""))
text(0.7, 1.02, cex=0.7, paste(round(data[2,1]*100, 2), "%", sep=""), xpd=TRUE)
text(1.9, 0.4, cex=0.7, paste(round(data[1,2]*100), "%", sep=""))
text(1.9, 0.95, cex=0.7, paste(round(data[2,2]*100), "%", sep=""))
text(3.1, 0.4, cex=0.7, paste(round(data[1,3]*100), "%", sep=""))
text(3.1, 0.9, cex=0.7, paste(round(data[2,3]*100), "%", sep=""))
#number and percentage
barplot(data, ylab="Proportion", col=c("cornflowerblue", "bisque4"))
legend(1, 1.15,c("Homozygous", "Heterozygous"),bty="n", xpd=TRUE, border="NA",lty=c(0,0),cex=1,fill=c("cornflowerblue", "bisque4"))
text(0.7, 0.4, cex=0.7, paste(x[1,1], "(", round(data[1,1]*100), "%)", sep=""))
text(0.7, 1.02, cex=0.7, paste(x[2,1], "(", round(data[2,1]*100,2), "%)", sep=""), xpd=TRUE)
text(1.9, 0.4, cex=0.7, paste(x[1,2], "(", round(data[1,2]*100), "%)", sep=""))
text(1.9, 0.95, cex=0.7, paste(x[2,2], "(", round(data[2,2]*100), "%)", sep=""))
text(3.1, 0.4, cex=0.7, paste(x[1,3], "(", round(data[1,3]*100), "%)", sep=""))
text(3.1, 0.9, cex=0.7, paste(x[2,3], "(", round(data[2,3]*100), "%)", sep=""))

dev.off()

pdf("mPing_class_proportion_unique.pdf", height=7, width=3)
#unique only
par(mar=c(5,5,4,5))
data1 <- data[,3]
y <- rbind(data1[1], data1[2])
rownames(y) <- c("Homozygous", "Heterzygous")
colnames(y) <- ""
#no data label
barplot(y, xlab="", ylab="", col=c("cornflowerblue", "bisque4"), axes=FALSE)
legend(0, 1.15,c("Homozygous", "Heterozygous"),bty="n", xpd=TRUE, border="NA",lty=c(0,0),cex=1,fill=c("cornflowerblue", "bisque4"))
#percentage
barplot(y, xlab="", ylab="", col=c("cornflowerblue", "bisque4"), axes=FALSE)
legend(0, 1.15,c("Homozygous", "Heterozygous"),bty="n", xpd=TRUE, border="NA",lty=c(0,0),cex=1,fill=c("cornflowerblue", "bisque4"))
text(0.7, 0.4, cex=0.8, paste(round(y[1]*100), "%", sep=""))
text(0.7, 0.9, cex=0.8, paste(round(y[2]*100), "%", sep=""))
#percentage and number
barplot(y, xlab="", ylab="", col=c("cornflowerblue", "bisque4"), axes=FALSE, cex.axis=1.4, cex.lab=1.4)
legend(0, 1.15,c("Homozygous", "Heterozygous"),bty="n", xpd=TRUE, border="NA",lty=c(0,0),cex=1.2,fill=c("cornflowerblue", "bisque4"))
#text(0.7, 0.4, cex=0.8, paste(x[1,3], "(", round(y[1]*100), "%)", sep=""))
#text(0.7, 0.9, cex=0.8, paste(x[2,3], "(", round(y[2]*100), "%)", sep=""))
#text(0.7, 0.4, cex=1.4, x[1,3], col='white')
#text(1.5, 0.4, cex=1.4, paste(round(y[1]*100), "%", sep=""), xpd=TRUE)
#text(0.7, 0.9, cex=1.4, x[2,3], col='white')
#text(1.5, 0.9, cex=1.4, paste(round(y[2]*100), "%", sep=""), xpd=TRUE)
text(0.7, 0.355, cex=1.4, paste("(", x[1,3], ")", sep=""), col='white')
text(0.7, 0.4, cex=1.4, paste(round(y[1]*100), "%", sep=""), xpd=TRUE, col='white')
text(0.7, 0.855, cex=1.4, paste("(", x[2,3], ")", sep=""), col='white')
text(0.7, 0.9, cex=1.4, paste(round(y[2]*100), "%", sep=""), xpd=TRUE, col='white')

text(0.7, -0.1, cex=1.2, 'Unique', xpd=TRUE)
#
dev.off()

