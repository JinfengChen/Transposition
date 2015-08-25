error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
###plot 100 kb
pdf("mPing_boundary.linked_50Mb_debug2.table_clean.pdf")
data =read.table("mPing_boundary.linked_50Mb_debug2.table_clean.sum100kb.txt")
data <- data[data[,1]<=1200000,]
excision = data[,4]
std = data[,5]
barx <- barplot(excision, col=c("black"), ylim=c(0, 12), border=F, axis.lty=1, xlab='', ylab='')
error.bar(barx, excision, std)
axis(1,c(0.1, max(barx)+0.6),line=0,labels=c("",""))
text(barx, rep(-0.5, 8),offset=2,labels=data[,1]/1000,srt=0,xpd=TRUE)
#legend("topright",c("HEG4","Nipponbare"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("blue","orange"))
xpos <- 6.5
ypos <- 5
#mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1, col="black")
mtext("Distance (kb)", side=1,font=1, at=xpos+1.3,line=3, cex=1, col="black")
#mtext("New", side=2,font=1, at=ypos,line=3, cex=1, col="black")
mtext("Excision", side=2,font=1, at=ypos+1,line=3, cex=1, col="black")
mtext("Number", side=2,font=1, at=ypos+2.5,line=3, cex=1, col="black")

###plot 10 kb
data =read.table("mPing_boundary.linked_50Mb_debug2.table_clean.sum10kb.txt")
data <- data[data[,1]<=100000,]
excision = data[,4]
std = data[,5]
barx <- barplot(excision, col=c("black"), ylim=c(0, 20), border=F, axis.lty=1, xlab='', ylab='')
error.bar(barx, excision, std)
axis(1,c(0.1, max(barx)+0.6),line=0,labels=c("",""))
text(barx, rep(-1, 8),offset=2,labels=data[,1]/1000,srt=0,xpd=TRUE)
#legend("topright",c("HEG4","Nipponbare"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("blue","orange"))
xpos <- 5
ypos <- 9
#mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1, col="black")
mtext("Distance (kb)", side=1,font=1, at=xpos+1.3,line=3, cex=1, col="black")
#mtext("New", side=2,font=1, at=ypos,line=3, cex=1, col="black")
mtext("Excision", side=2,font=1, at=ypos+1,line=3, cex=1, col="black")
mtext("Number", side=2,font=1, at=ypos+3.5,line=3, cex=1, col="black")

dev.off()


