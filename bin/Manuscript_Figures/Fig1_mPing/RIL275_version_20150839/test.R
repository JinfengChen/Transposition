pdf("mPing_Call_Correlation_unique_clean_RIL230.pdf")
relocate <- read.table("RIL275_RelocaTE.sofia.sorted.table", skip=1)
relocatei <- read.table("RIL275_RelocaTEi.summary_clean.unique.table", skip=1)
temp <- read.table("RIL275_TEMP.summary.table", skip=1)

#blackout low confidence RILs
bl <-read.table("Bam.Core.blacklist")
sel <- match(relocate[,1], bl[,1])
relocate <- relocate[is.na(sel),]
relocatei <- relocatei[is.na(sel),]
temp <- temp[is.na(sel),]

##RelocaTEi and TEMPv1.02
plot(relocatei[,7], temp[,7], xlim=c(0,600),ylim=c(0, 600),type="p",pch=20,col="cadetblue",xlab="RelocaTEi",ylab="TEMPv1.02")
points(relocatei[,8], temp[,8], type="p",pch=20,col="chocolate")
segments(-100, -100, 1100, 1100, col='black', lty=1,lwd = 0.5)
legend('topright',bty="n",lty=c(0,0), pch=c(20, 20), col=c("cadetblue","chocolate"), cex=1.2,c("ALL insertion", "Confident insertion"))

##RelocaTEi and RelocaTE
plot(relocatei[,7], relocate[,3], xlim=c(0,1000),ylim=c(0, 1000),type="p",pch=20,col="cadetblue",xlab="RelocaTEi",ylab="RelocaTE")
points(relocatei[,8], relocate[,3], xlim=c(0,1000),ylim=c(0, 1000),type="p",pch=20,col="chocolate")
segments(-100, -100, 1100, 1100, col='black', lty=1,lwd = 0.5)
legend('topright',bty="n",lty=c(0,0), pch=c(20, 20), col=c("cadetblue","chocolate"), cex=1.2,c("ALL insertion", "Confident insertion"))


##set two plots in one figure
par(mfrow=c(2,1))
par(mar=c(5,4,2,2))
##RelocaTEi and insertion size
relocatei_sortbysize <- relocatei[order(relocatei[,2]),]
relocatei_calls <- rbind(relocatei_sortbysize[,8], relocatei_sortbysize[,9])
barx <- barplot(relocatei_calls, ylim=c(0, 800), border=FALSE, col=c("cadetblue","chocolate"), xlab='Insertion Size (bp)', ylab='mPing insertions', main="RelocaTE2")
step <- (max(barx)+0.6-0.9)/275
pos <- c(seq(1, 275, by=23), 275)
pos_marks <- step*pos
axis(1,pos_marks, tck=-0.04, line=0.2,labels=rep("",13))
axis(1,c(0.9,pos_marks[length(pos)]), line=0.2,labels=c("",""))
marks <- as.integer(relocatei_sortbysize[,2][pos])
text(pos_marks, rep(-100, 13),offset=2,labels=as.integer(marks),srt=0,xpd=TRUE)
legend("topright",c("Confident insertion", "Potential insertion"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("cadetblue","chocolate"))

#step <- (max(barx)+0.6-0.9)/10
#axis(1,seq(0.9, max(barx)+0.6, by=step), tck=-0.04, line=0.2,labels=rep("",11))
#axis(1,c(0.9,max(barx)+0.6), line=0.2,labels=c("",""))
#insertstep <- (max(relocatei[,2])-min(relocatei[,2]))/10
#marks <- seq(min(relocatei[,2]), max(relocatei[,2]), by = insertstep)
#text(seq(0.9, max(barx)+0.6, by=step), rep(-100, 11),offset=2,labels=as.integer(marks),srt=0,xpd=TRUE)

##RelocaTEi and depth
relocatei_sortbydepth <- relocatei[order(relocatei[,5]),]
relocatei_calls <- rbind(relocatei_sortbydepth[,8], relocatei_sortbydepth[,9])
barx <- barplot(relocatei_calls, ylim=c(0, 800), border=FALSE, col=c("cadetblue","chocolate"), xlab='Sequence Depth (X)', ylab='mPing insertions')
step <- (max(barx)+0.6-0.9)/230
pos <- c(seq(1, 230, by=23), 230)
pos_marks <- step*pos
axis(1,pos_marks, tck=-0.04, line=0.2,labels=rep("",11))
axis(1,c(0.9,pos_marks[length(pos)]), line=0.2,labels=c("",""))
marks <- as.integer(relocatei_sortbydepth[,5][pos])
text(pos_marks, rep(-100, 11),offset=2,labels=as.integer(marks),srt=0,xpd=TRUE)
legend("topright",c("Confident insertion", "Potential insertion"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("cadetblue","chocolate"))

dev.off()
