error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
###plot 100 kb
pdf("mPing_boundary.linked_50Mb_debug2.table_clean.proportion_accumulation.pdf")
data =read.table("mPing_boundary.linked_50Mb_debug2.distance_accumulation_excision.list", header=T)
accumulated_excision_freq = data[,4]
control_freq = data[,7]
dist = data[,1]
plot(dist/1000000, accumulated_excision_freq, col='black', pch=19, cex=1, xlab="Distance (Mb)", ylab="Accumulated Fraction of Excision", xlim=c(0, 2000000/1000000), ylim=c(0, 1))
points(dist/1000000, control_freq, col='blue', pch=19, cex=1)
#axis(1,c(0.1, 1000000),line=0,labels=c("",""))
#text(barx, rep(-0.5, 8),offset=2,labels=data[,1]/1000,srt=0,xpd=TRUE)
legend(1600000/1000000, 0.8, c("Excisions", "Control"), bty="n", border="NA", pch=c(19, 19), cex=1, col=c("black", "blue"))

plot(dist/1000, accumulated_excision_freq, col='black', pch=19, cex=1, xlab="Distance (kb)", ylab="Accumulated Fraction of Excision", xlim=c(0, 500000/1000), ylim=c(0, 1))
points(dist/1000, control_freq, col='blue', pch=19, cex=1)
legend(400000/1000, 0.6, c("Excisions", "Control"), bty="n", border="NA", pch=c(19, 19), cex=1, col=c("black", "blue"))

#xpos <- 6.5
#ypos <- 5
#mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1, col="black")
#mtext("Distance (kb)", side=1,font=1, at=xpos+1.3,line=3, cex=1, col="black")
#mtext("New", side=2,font=1, at=ypos,line=3, cex=1, col="black")
#mtext("Excision", side=2,font=1, at=ypos+1,line=3, cex=1, col="black")
#mtext("Number", side=2,font=1, at=ypos+2.5,line=3, cex=1, col="black")


dev.off()


