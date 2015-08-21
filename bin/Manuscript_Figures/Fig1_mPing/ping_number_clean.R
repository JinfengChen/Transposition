error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x,y+upper, x, y, angle=90, code=1, length=length, ...)
}

pdf("ping_number_clean.pdf")
par(mar=c(5,4,4,5))
data =read.table("RIL275_RelocaTEi.CombinedGFF.characterized.clean.ping_number.summary")
expr = rbind(data[,2], data[,4]*1.71)
std = rbind(data[,3], data[,5]*1.71)
barx <- barplot(expr, beside=TRUE, ylim=c(0,140), border=F, axis.lty=1, xlab='', ylab='', axes=FALSE, col=c("darkseagreen","cornflowerblue"))
error.bar(barx, expr, std)
axis(1,c(0.7, max(barx)+0.8),line=0,labels=c("",""))
axis(2,c(0,20,40,60,80,100,120,140),labels=c(0, 20, 40, 60, 80, 100, 120, 140), col="darkseagreen", col.axis="darkseagreen")
axis(4,c(0,20,40,60,80,100,120,140),labels=c(0, 10, 20, 30, 40, 50, 60, 70), col="cornflowerblue", col.axis="cornflowerblue")

text(seq(2, 25, by=3), rep(-5, 6),offset=2,labels=data[,1],srt=0,xpd=TRUE)
legend("topleft",c("Unique_hom","Unique_som/het"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("darkseagreen","cornflowerblue"))
xpos <- 5.6
ypos_l <- 60
mtext("Ping", side=1,font=3, at=xpos+6.8,line=2, cex=1, col="black")
mtext("Number", side=1,font=1, at=xpos+9.5,line=2, cex=1, col="black")
#left y axis
mtext("Unique_hom", side=2,font=1, at=ypos_l,line=3, cex=1, col="darkseagreen")
mtext("mPing", side=2,font=3, at=ypos_l+20,line=3, cex=1, col="darkseagreen")
mtext("number", side=2,font=1, at=ypos_l+36,line=3, cex=1, col="darkseagreen")
#right y axis
ypos_r <- 60
mtext("Unique_som/het", side=4,font=1, at=ypos_r,line=3, cex=1, col="cornflowerblue")
mtext("mPing", side=4,font=3, at=ypos_r+24,line=3, cex=1, col="cornflowerblue")
mtext("number", side=4,font=1, at=ypos_r+40,line=3, cex=1, col="cornflowerblue")

dev.off()


