error.bar <- function(x, y, upper, lower=upper, color,length=0.06,...){
     if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
     stop("vectors must be same length")
     arrows(x,y+upper, x, y-lower, col=color,angle=90, code=3, length=length, ...)
 }


pdf("simulate_intron_distance_point_withsim.pdf")

par(mar=c(6,4,4,2), cex=1.2)
sim <- read.table("../mPing_distr/results_simulationV2_TSD9mer_somatic/Simulate.TSD9mer.SomaticMat.intron.distance.distr")
som <- read.table("../mPing_distr/Somatic.intron.distance.distr")
str <- read.table("../mPing_distr/Strains.intron.distance.distr")
ril <- read.table("../mPing_distr/RIL.intron.distance.distr")


plot(som[,4], type='b', pch= 1,lwd = 2 , col="aquamarine3", xaxt='n', frame.plot = FALSE, ylim=c(0,0.4), ylab="Proportion", xlab="")
lines(ril[,4], type='b',pch= 2,lwd = 2 , col="steelblue2")
lines(str[,4], type='b',pch= 3,lwd = 2 , col="sandybrown")
lines(sim[,4], type='b',pch= 20, cex=0.2, lwd = 2 , col="dim gray")
error.bar(1:length(sim[,4]), sim[,4], sim[,7]-sim[,4], sim[,7]-sim[,4], 1)

axis(1,sim[,1],line=0, labels=rep("",length(sim[,1])))
text(sim[,1],rep(-0.06,7), cex=1, offset=2,labels=sim[,2],srt=55,xpd=TRUE)

legend('topright', bty='n', border='NA', lty= c(1,2,3,4), pch = c(1,2,3,20), cex=1 , lwd = 2 ,col=c("aquamarine3", "steelblue2", "sandybrown", "dim gray"), c("Somatic", "RIL", "Landrace","Simulation"))
mtext("Distance to Exon (bp)", side=1,cex=1.2, at=4.5,line=4)



dev.off()

