error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
###plot 100 kb
pdf("mPing_boundary.linked_50Mb_debug2.table_clean.point_view.pdf")
data =read.table("mPing_boundary.linked_50Mb_debug2.excision_distance.list")
excision = data[,3]
dist = data[,2]
plot(dist, excision, col='black', pch=19, cex=1, xlab="Distance (kb)", ylab="#Excision", xlim=c(0, 1000000))
#axis(1,c(0.1, 1000000),line=0,labels=c("",""))
#text(barx, rep(-0.5, 8),offset=2,labels=data[,1]/1000,srt=0,xpd=TRUE)
#legend("topright",c("HEG4","Nipponbare"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("blue","orange"))
xpos <- 6.5
ypos <- 5
#mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1, col="black")
#mtext("Distance (kb)", side=1,font=1, at=xpos+1.3,line=3, cex=1, col="black")
#mtext("New", side=2,font=1, at=ypos,line=3, cex=1, col="black")
#mtext("Excision", side=2,font=1, at=ypos+1,line=3, cex=1, col="black")
#mtext("Number", side=2,font=1, at=ypos+2.5,line=3, cex=1, col="black")


dev.off()


