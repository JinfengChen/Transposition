error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
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

pdf("LOC_Os01g62630.pdf", width=5, height=7)
#png("LOC_Os01g62630.png", width=500, height=700)
par(mar=c(5,5,4,2))
data =read.table("LOC_Os01g62630.rpkm", skip=1)
expr = data[,1][1:3]
std = data[,1][4:6]
barx <- barplot(expr, space=1, col=c("black"), border=F, ylim=c(0,(max(expr)+max(std))*1.5), cex.axis = 1.5, cex.name = 1.5, cex.lab = 1.5, axis.lty=1, ylab="Expression (RPKM)")
error.bar(barx, expr, std)
axis(1,c(0.2,max(barx)+0.6),line=0,labels=c("",""))
text(barx,rep(-max(expr)*0.25,3),offset=2, cex=1.5, labels=c("NB", "HEG4", "EG4"),srt=0,xpd=TRUE)
#legend("topright",c("HEG4","Nipponbare"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("blue","orange"))

p1 = 0.1
y11 = (data[,1][1] + data[,1][4])*1.1 
x11 = barx[1]
y12 = (data[,1][2] + data[,1][5])*1.1
x12 = barx[2]-0.1
top1 = max(y11, y12)*1.1
p2 = 0.1
y21 = (data[,1][2] + data[,1][5])*1.1
x21 = barx[2] + 0.1
y22 = (data[,1][3] + data[,1][6])*1.1
x22 = barx[3]
top2 = max(y21, y22)*1.1

pvalue(x11,y11,x12,y12,top1, p1)
pvalue(x21,y21,x22,y22,top2, p2)
dev.off()


