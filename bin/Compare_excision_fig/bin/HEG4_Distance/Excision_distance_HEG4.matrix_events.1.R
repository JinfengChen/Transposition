scatterhist = function(x, y, xlab="", ylab=""){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  breaks = c(seq(0, 1200000/1000, by=10000/1000))
  xhist = hist(x, breaks=breaks, plot=FALSE)
  breaks = c(seq(-1, 20, by=1))
  yhist = hist(y, breaks=breaks, plot=FALSE)
  topx = max(c(xhist$counts))*1.1
  topy = max(c(yhist$counts))*1.1
  par(mar=c(3,3,1,1))
  plot(x,y, pch=18)
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, topx), space=0)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, topy), space=0, horiz=TRUE)
  par(oma=c(3,3,0,0))
  #mtext(xlab, side=1, line)
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
        at=1.4 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0,
        at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
}
pdf("Excision_distance_HEG4.matrix_events.1.pdf", width=7, height=7)
#data <- read.table("Excision_distance.matrix.1.txt", header=T)
#data <- read.table("Excision_distance.matrix.2.txt", header=T)
data <- read.table("Excision_distance_HEG4.matrix_events.1.txt", header=T)
data <- data[data[,3]<1200000,]
#scatterhist(data[,3]/1000, data[,2], xlab="Distance (kb)", ylab="Number of Excision")
x <- data[,3]/1000
y <- data[,2]
xlab="Distance (kb)"
ylab="Number of Excision"
zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
breaks = c(seq(0, 1200000/1000, by=10000/1000))
xhist = hist(x, breaks=breaks, plot=FALSE)
breaks = c(seq(-1, 20, by=1))
yhist = hist(y, breaks=breaks, plot=FALSE)
topx = max(c(xhist$counts))*1.1
topy = max(c(yhist$counts))*1.1
par(mar=c(3,3,1,1))
plot(x,y, pch=18)
par(mar=c(0,3,1,1))
barplot(xhist$counts, axes=FALSE, ylim=c(0, topx), space=0)
par(mar=c(3,0,1,1))
barplot(yhist$counts, axes=FALSE, xlim=c(0, topy), space=0, horiz=TRUE)
par(oma=c(3,3,0,0))
#mtext(xlab, side=1, line)
mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
      at=1.5 * (mean(x) - min(x))/(max(x)-min(x)))
mtext(ylab, side=2, line=1, outer=TRUE, adj=0,
      at=(5 * (mean(y,na.rm=TRUE) - min(y,na.rm=TRUE))/(max(y,na.rm=TRUE) - min(y,na.rm=TRUE))))
dev.off()
