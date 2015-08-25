error.bar <- function(x, y, upper, lower=upper, color,length=0.06,...){
     if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
     stop("vectors must be same length")
     arrows(x,y+upper, x, y-lower, col=color,angle=90, code=3, length=length, ...)
 }


pdf("mping_position_barplot_withsim.pdf")

par(mar=c(6,4,4,2), cex=1.2)
som3 <- read.table("../Somatic.position.distr")
str3 <- read.table("../Strains.position.distr")
ril3 <- read.table("../RIL.position.distr")
sim3 <- read.table("../Simulate.TSD9mer.SomaticMat.position.distr")

som3 <- som3[-1,]
str3 <- str3[-1,]
ril3 <- ril3[-1,]
sim3 <- sim3[-1,]

som3 <- som3[-length(som3[,1]),]
str3 <- str3[-length(str3[,1]),]
ril3 <- ril3[-length(ril3[,1]),]
sim3 <- sim3[-length(sim3[,1]),]


data <- rbind(som3[,4], str3[,4], ril3[,4], sim3[,4])
xx <- barplot(data,beside=TRUE,ylab="Proportion",cex.names=1,cex.axis=1,border=FALSE,ylim=c(0,0.2),col=c("aquamarine3", "steelblue2" ,"sandybrown","dim gray"))

error.bar(xx[4,], sim3[,4], sim3[,7]-sim3[,4], sim3[,7]-sim3[,4], 1)

axis(1,c(0.5,max(xx)+0.5),line=0,labels=c("",""))
text(xx[1,]+0.4,rep(-0.033,7), cex=1, offset=2,labels=som3[,2],srt=55,xpd=TRUE)
legend('topright', bty='n', border='NA', lty=c(0,0),cex=1 ,c("Somatic", "Strains", "RIL", "Simulation"),fill=c("aquamarine3", "steelblue2" ,"sandybrown", "dim gray"))



dev.off()

