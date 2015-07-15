pdf("MSU7.SlidingWin.HEG4_mPing_excision.pdf", width = 10, height = 14)
par(mfrow=c(12, 1))
par(mar=c(2,4,1,1))
x <- read.table("MSU7.SlidingWin.HEG4.density")
y <- read.table("MSU7.SlidingWin.excision.density")
z <- read.table("MSU7.high_excision.list")

for (i in paste("Chr", 1:12, sep="")){
    print(i)
    ##unique/excision mPing
    test <- match(x[,1], i)
    selected <- which(!is.na(test))
    u <- x[selected,]
    e <- y[selected,]
    plot(u[,2], u[,7], col='cadetblue', type='l', main="", ylab="Count", xlim=c(0, 45000000), ylim=c(0, 6), frame.plot=FALSE, axes=FALSE)
    lines(e[,2], e[,7]*2, col='chocolate', xlim=c(0, 45000000), ylim=c(0,6))
    ##y axis
    aty = c(0, 2, 4, 6)
    aty1 = c("/0", "/1", "/2", "/3")
    axis(2,at=aty,labels=rep("",length(aty)))
    text(rep(-3200000,length(aty)), aty, offset=2, labels=aty, col='cadetblue', srt=0,xpd=TRUE)
    text(rep(-2700000,length(aty)), aty, offset=2, labels=aty1, col='chocolate', srt=0,xpd=TRUE)
    ##x axis
    xmax <- ceiling(u[length(u[,1]),][,3]/1000000)
    atx <- c(seq(0, xmax, 2))
    if (atx[length(atx)] <= xmax - 1){
        atx <- append(atx, xmax)
    }
    axis(1, at=atx*1000000, labels=atx, cex=0.3)
    mtext("(Mb)", side=1, line = 1, at = (xmax+1)*1000000, cex=0.6)
    ##title
    mtext(i, side=3, line = -1, at = (18)*1000000, cex=1)
    ##high excision mPing
    test1 <- match(z[,1], i)
    selected1 <- which(!is.na(test1))
    h <- z[selected1,]
    for (j in h[,2]){
        arrows(j+1000000,6,j,5,length=0.1,angle=10,col='red',lwd=0.5)
    }
    ##legend
    if (i == "Chr2"){
        legend(40*1000000, 4, c("Insertion", "Excision"), ,bty="n", border="NA",lty=c(1, 1), cex=1, col=c("cadetblue","chocolate"), text.col=c("cadetblue","chocolate"))
        text(43.8*1000000, 0.5, "High Excision", cex=1, col="red")
        arrows(41.5*1000000, 0.5, 40.5*1000000, 0.5, length=0.1, angle=10, col='red',lwd=0.5)
    }
}

dev.off()
