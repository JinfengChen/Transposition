error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
error.text <- function(x, y, n){
    text(x, y, labels=paste("n=", n, sep=''), xpd=TRUE)
}

pvalue <- function(x1, y1, x2, y2, top, p){
     star <- '*'
     if (p > 0.05) {star <- 'n.s.'}
     if (p < 0.001){ star <- '**'}
     if (p < 0.0001){ star <- '***'}
     segments(x1,y1,x1,top)
     segments(x1,top,x2,top)
     segments(x2,top,x2,y2)
     #segments(x1-0.2,y1,x1+0.2,y1)
     #segments(x2-0.2,y2,x2+0.2,y2)
     xt <- min(x1,x2)+abs(x2-x1)/2
     yt <- top*1.1
     text(xt,yt,star, cex=1.5)
} 


###plot 100 kb
pdf("mPing_boundary.linked_50Mb_debug2.table_clean.pdf")
data =read.table("mPing_boundary.linked_50Mb_debug2.table_clean.sum100kb.txt")
data <- data[data[,1]<=1200000,]
excision = data[,4]
std = data[,5]
events = data[,2]
barx <- barplot(excision, col=c("black"), ylim=c(0, 12), border=F, axis.lty=1, xlab='', ylab='')
error.bar(barx, excision, std)
error.text(barx, excision+std+0.5, events)
axis(1,c(0.1, max(barx)+0.6),line=0,labels=c("",""))
text(barx, rep(-0.5, 8),offset=2,labels=data[,1]/1000000,srt=0,xpd=TRUE)
#legend("topright",c("HEG4","Nipponbare"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("blue","orange"))
xpos <- 6.5
ypos <- 5
#mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1, col="black")
mtext("Distance (Mb)", side=1,font=1, at=xpos+1.3,line=3, cex=1, col="black")
#mtext("New", side=2,font=1, at=ypos,line=3, cex=1, col="black")
mtext("Excision", side=2,font=1, at=ypos+1,line=3, cex=1, col="black")
mtext("Number", side=2,font=1, at=ypos+2.5,line=3, cex=1, col="black")

###plot 10 kb
data =read.table("mPing_boundary.linked_50Mb_debug2.table_clean.sum10kb.txt")
data <- data[data[,1]<=100000,]
excision = data[,4]
std = data[,5]
events = data[,2]
barx <- barplot(excision, col=c("black"), ylim=c(0, 22), border=F, axis.lty=1, xlab='', ylab='')
error.bar(barx, excision, std)
error.text(barx, excision+std+0.5, events)
axis(1,c(0.1, max(barx)+0.6),line=0,labels=c("",""))
text(barx, rep(-1, 8),offset=2,labels=data[,1]/1000,srt=0,xpd=TRUE)
#legend("topright",c("HEG4","Nipponbare"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("blue","orange"))
xpos <- 5
ypos <- 9
#mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1, col="black")
mtext("Distance (kb)", side=1,font=1, at=xpos+1.3,line=3, cex=1, col="black")
#mtext("New", side=2,font=1, at=ypos,line=3, cex=1, col="black")
mtext("Excision", side=2,font=1, at=ypos+1,line=3, cex=1, col="black")
mtext("Number", side=2,font=1, at=ypos+3.8,line=3, cex=1, col="black")

dev.off()


