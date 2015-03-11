pdf('mping.excision.frq.pdf')
t <- read.table("mping.excision.table")
f <- t[,3]
f <- f[f>0]
xx <- hist(f,ylab="Insertion sites",xlab="Allele Frequency", col="steelblue2", main="",border=T, plot=F)
y <- xx$counts
x <- xx$mids

xx <- barplot(y,beside=TRUE,ylab="Insertion sites",cex.lab=1.2, cex.names=1.2, cex.axis=1.2,border=FALSE,axes=FALSE,ylim=c(0,100),col="steelblue2")

b <- 220

axis(1,c(0,2.5,4.8,7.3,9.7,11),line=0,c(0.1,0.2,0.3,0.4,0.5,""))
axis(2,c(0,20,40,60,80,100),line=0,labels=c(0,20,40,60,80,100))
rect(0.2,b,max(xx)+0.8,b+5,border=FALSE,col='white')
#text(c(0,2.5,4.8,7.3,9.7,12.1,14.5,16.9),rep(-45,length(y)), cex=1, offset=2,labels=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),srt=0,xpd=TRUE)
mtext("Allele Frequency", side=1,cex=1.2, at=6,line=3)

dev.off()
