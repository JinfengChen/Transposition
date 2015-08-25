error.bar <- function(x, y, upper, lower=upper, color,length=0.06,...){
     if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
     stop("vectors must be same length")
     arrows(x,y+upper, x, y-lower, col=color,angle=90, code=3, length=length, ...)
 }


pdf("mping_intergenic_3distance_withsim.pdf")

par(mar=c(6,4,4,2), cex=1.2)
som3 <- read.table("../mPing_distr/Somatic.mRNA.3primer.distance.distr")
str3 <- read.table("../mPing_distr/Strains.mRNA.3primer.distance.distr")
ril3 <- read.table("../mPing_distr/RIL.mRNA.3primer.distance.distr")
sim3 <- read.table("../mPing_distr/results_simulationV2_TSD9mer_somatic/Simulate.TSD9mer.SomaticMat.mRNA.3primer.distance.distr")

som3 <- som3[-1,]
str3 <- str3[-1,]
ril3 <- ril3[-1,]
sim3 <- sim3[-1,]

som3 <- som3[-length(som3[,1]),]
str3 <- str3[-length(str3[,1]),]
ril3 <- ril3[-length(ril3[,1]),]
sim3 <- sim3[-length(sim3[,1]),]


plot(som3[,4], type='b', pch= 1,lwd = 2 , col="aquamarine3", xaxt='n', frame.plot = FALSE, ylim=c(0,0.2), ylab="Proportion", xlab="")
lines(ril3[,4], type='b',pch= 2,lwd = 2 , col="steelblue2")
lines(str3[,4], type='b',pch= 3,lwd = 2 , col="sandybrown")
lines(sim3[,4], type='b',pch= 20, cex=0.2, lwd = 2 , col="dim gray")
error.bar(1:length(sim3[,4]), sim3[,4], sim3[,7]-sim3[,4], sim3[,7]-sim3[,4], 1)

axis(1,seq(1:length(som3[,1])),line=0, labels=rep("",length(som3[,1])))
#text(som3[,1],rep(-0.03,7), cex=1, offset=2,labels=som3[,2],srt=55,xpd=TRUE)
text(seq(1:length(som3[,1][-length(som3[,1])]))+0.5,rep(-0.02,7), cex=1, offset=2,labels=som3[,1][-length(som3[,1])]*500/-1000,srt=55,xpd=TRUE)

legend('topright', bty='n', border='NA', lty= c(1,2,3,4), pch = c(1,2,3,20), cex=1 , lwd = 2 ,col=c("aquamarine3", "steelblue2", "sandybrown", "dim gray"), c("Somatic", "RIL", "Landrace", "Simulation"))
mtext("Distance to TTS (kp)", side=1,cex=1.2, at=9,line=3)

dev.off()

