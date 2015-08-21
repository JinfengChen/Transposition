error.bar <- function(x, y, upper, lower=upper, color,length=0.06,...){
     if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
     stop("vectors must be same length")
     arrows(x,y+upper, x, y-lower, col=color,angle=90, code=3, length=length, ...)
 }

pdf("mping_intergenic_5distance_withsim.pdf")

par(mar=c(6,4,4,2), cex=1.2)
som5 <- read.table("../mPing_distr/Somatic.mRNA.5primer.distance.distr")
str5 <- read.table("../mPing_distr/Strains.mRNA.5primer.distance.distr")
ril5 <- read.table("../mPing_distr/RIL.mRNA.5primer.distance.distr")
sim5 <- read.table("../mPing_distr/results_simulationV2_TSD9mer_somatic/Simulate.TSD9mer.SomaticMat.mRNA.5primer.distance.distr")

som5 <- som5[-1,]
str5 <- str5[-1,]
ril5 <- ril5[-1,]
sim5 <- sim5[-1,]

som5 <- som5[-length(som5[,1]),]
str5 <- str5[-length(str5[,1]),]
ril5 <- ril5[-length(ril5[,1]),]
sim5 <- sim5[-length(sim5[,1]),]

plot(rev(som5[,4]), type='b', pch= 1,lwd = 2 , col="aquamarine3", xaxt='n', frame.plot = FALSE, ylim=c(0,0.2), ylab="Proportion", xlab="")
lines(rev(ril5[,4]), type='b',pch= 2,lwd = 2 , col="steelblue2")
lines(rev(str5[,4]), type='b',pch= 3,lwd = 2 , col="sandybrown")
lines(rev(sim5[,4]), type='b',pch= 20, cex=0.2,lwd = 2 , col="dim gray")
error.bar(1:length(sim5[,4]), rev(sim5[,4]), rev(sim5[,7]-sim5[,4]), rev(sim5[,7]-sim5[,4]), 'dim gray')

#yaxis <- seq(1:length(som5[,1])+0.5
axis(1,seq(1:length(som5[,1])),line=0, labels=rep("",length(som5[,1])))
text(seq(1:length(som5[,1][-1]))+0.5,rep(-0.02,7), cex=1, offset=2,labels=rev(som5[,1]*500/-1000)[-1],srt=55,xpd=TRUE)

legend('topright', bty='n', border='NA', lty= c(1,2,3,4), pch = c(1,2,3,20), cex=1 , lwd = 2 ,col=c("aquamarine3", "steelblue2", "sandybrown", "dim gray"), c("Somatic", "RIL", "Strains", "Simulation"))
mtext("Distance to TSS (kp)", side=1,cex=1.2, at=9,line=3)



dev.off()

