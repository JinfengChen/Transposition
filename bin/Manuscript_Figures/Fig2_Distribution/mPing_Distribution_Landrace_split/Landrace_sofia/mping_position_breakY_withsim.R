error.bar <- function(x, y, upper, lower=upper, color,length=0.06,...){
     if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
     stop("vectors must be same length")
     arrows(x,y+upper, x, y-lower, col=color,angle=90, code=3, length=length, ...)
 }


pdf("mping_position_breakY_withsim.pdf")
library("plotrix")
par(mar=c(6,4,4,2), cex=1.2)
x1 <- read.table("../mPing_distr/Somatic.position.distr")
x2 <- read.table("../mPing_distr/RIL.position.distr")
x3 <- read.table("../mPing_distr/Strains.position.distr")
x4 <- read.table("../mPing_distr/results_simulationV2_TSD9mer_somatic/Simulate.TSD9mer.SomaticMat.position.distr")

x1 <- x1[-length(x1[,1]),]
x2 <- x2[-length(x1[,1]),]
x3 <- x3[-length(x1[,1]),]
x4 <- x4[-length(x1[,1]),]

x1p <- x1[,3]/sum(x1[,3])
x2p <- x2[,3]/sum(x2[,3])
x3p <- x3[,3]/sum(x3[,3])
x4p <- x4[,3]/sum(x4[,3])

cut <- 0.5
x1p[1] <- x1p[1] - cut
x2p[1] <- x2p[1] - cut
x3p[1] <- x3p[1] - cut
x4p[1] <- x4p[1] - cut

data <- rbind(x1p, x2p, x3p, x4p)
xx <- barplot(data,beside=TRUE,ylab="Proportion",cex.names=1,cex.axis=1,border=FALSE,axes=FALSE,ylim=c(0,0.3),col=c("aquamarine3", "steelblue2" ,"sandybrown", "dim gray"))
error.bar(xx[4,], x4p, x4[,7]-x4[,4], x4[,7]-x4[,4], 1)

b <- 0.19

axis(1,c(0.5,max(xx)+0.5),line=0,labels=c("",""))
axis(2,c(0,0.1,0.2,0.3),line=0,labels=c(0,0.1,0.7,0.8))
axis.break(2, b,style="slash")
rect(0.5,b,max(xx)+0.5,b+0.005,border=FALSE,col='white')
text(xx[1,]+0.2,rep(-0.045,7), cex=1, offset=2,labels=x1[,2],srt=55,xpd=TRUE)
legend('topright', bty='n', border='NA', lty=c(0,0),cex=1 ,c("Somatic", "RIL","Strains", "Simulation"),fill=c("aquamarine3", "steelblue2" ,"sandybrown", "dim gray"))



dev.off()

