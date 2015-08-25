pdf("ping_mping_number_excision.pdf", height=7, width=8)
x <- read.table("Excision_transposases.table.txt", header=T)
##ping, confident all non-reference mPing vs. excision
symbols(x[,14], x[,8], circles=x[,15], inches=0.35, fg="white", bg="red", axes=FALSE, xlab="", ylab="", xlim=c(0, 8), ylim=c(0, 600))
axis(1, seq(0, 8, by=1), line=0)
axis(2, seq(0, 601, by=100), line=0)
xpos <- 3.6
ypos <- 220
mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1, col="black")
mtext("Number", side=1,font=1, at=xpos+1,line=3, cex=1, col="black")
mtext("mPing", side=2,font=3, at=ypos+1,line=3, cex=1, col="black")
mtext("Number", side=2,font=1, at=ypos+80,line=3, cex=1, col="black")
x1 <- c(7, 7, 7)
y1 <- c(40, 55, 75)
le <- c("1", "5", "10")
size <- c(0.6, 4, 8)
for (i in c(1,2,3)){
    legend(x=x1[i], y=y1[i], 
        legend="", 
        pch = 1, 
        bty = "n", 
        col = "black",
        xpd = TRUE,
        pt.bg = "red",         
        pt.cex = size[i]
    )
}

text(7.2, 28, labels='1')
text(7.2, 58, labels='5')
text(7.2, 100, labels='10')
text(7.2, 130, labels='Excision')
##ping, confident parental non-reference mPing vs. excision
symbols(x[,14], x[,16], circles=x[,15], inches=0.35, fg="white", bg="red", axes=FALSE, xlab="", ylab="", xlim=c(0, 8), ylim=c(0, 300))
axis(1, seq(0, 8, by=1), line=0)
axis(2, seq(0, 300, by=50), line=0)
xpos <- 3.6
ypos <- 120
mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1, col="black")
mtext("Number", side=1,font=1, at=xpos+1,line=3, cex=1, col="black")
mtext("mPing", side=2,font=3, at=ypos+1,line=3, cex=1, col="black")
mtext("Number", side=2,font=1, at=ypos+40,line=3, cex=1, col="black")
x1 <- c(7, 7, 7)
y1 <- c(25, 32, 41)
le <- c("1", "5", "10")
size <- c(0.6, 4, 8)
for (i in c(1,2,3)){
    legend(x=x1[i], y=y1[i], 
        legend="", 
        pch = 1, 
        bty = "n", 
        col = "black",
        xpd = TRUE,
        pt.bg = "red",         
        pt.cex = size[i]
    )
}

text(7.2, 22, labels='1')
text(7.2, 35, labels='5')
text(7.2, 53, labels='10')
text(7.2, 66, labels='Excision')
##ping, confident parental non-reference mPing 100 kb vs. excision
symbols(x[,14], x[,17], circles=x[,15], inches=0.35, fg="white", bg="red", axes=FALSE, xlab="", ylab="", xlim=c(0, 8), ylim=c(0, 80))
axis(1, seq(0, 8, by=1), line=0)
axis(2, seq(0, 80, by=20), line=0)
xpos <- 3.6
ypos <- 40
mtext("Ping", side=1,font=3, at=xpos+0.3,line=3, cex=1, col="black")
mtext("Number", side=1,font=1, at=xpos+1,line=3, cex=1, col="black")
mtext("mPing", side=2,font=3, at=ypos+1,line=3, cex=1, col="black")
mtext("Number", side=2,font=1, at=ypos+11,line=3, cex=1, col="black")
x1 <- c(7, 7, 7)
y1 <- c(7, 9, 11.5)
le <- c("1", "5", "10")
size <- c(0.6, 4, 8)
for (i in c(1,2,3)){
    legend(x=x1[i], y=y1[i], 
        legend="", 
        pch = 1, 
        bty = "n", 
        col = "black",
        xpd = TRUE,
        pt.bg = "red",         
        pt.cex = size[i]
    )
}

text(7.2, 6, labels='1')
text(7.2, 10, labels='5')
text(7.2, 15, labels='10')
text(7.2, 19, labels='Excision')
dev.off()
