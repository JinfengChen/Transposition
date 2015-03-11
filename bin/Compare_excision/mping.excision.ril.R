pdf("mping.excision.ril.transposition.pdf")

par(mar=c(6,4,4,2), cex=1.2)
x1 <- read.table("mping.excision.ril.transposition")

plot(x1[,5],x1[,2],font.lab=1,ylab="Ping Number", xlab="Excision Number", frame.plot=F, pch=16, cex=0.6, xlim=c(0,10), ylim=c(0,7),col=c("steelblue2"))
m1 <- lm(x1[,2]~x1[,5])
eq <- paste('Y', '=', sprintf('%.3f',m1$coefficients[[1]]), '+', paste(sprintf('%.3f',m1$coefficients[[2]]), 'X', sep=''), sep=' ')
abline(m1, col='orange')
text(8,4, labels=eq,xpd=TRUE)
cor0 <- cor.test(x1[,2], x1[,5])
cl <- paste('R =', sprintf('%.2f',cor0$estimate), ',', 'P-value =', sprintf('%.3g',cor0$p.value))
text(8,4.5, labels=cl,xpd=TRUE)

plot(x1[,5],x1[,3],font.lab=1,ylab="mPing Number", xlab="Excision Number", frame.plot=F, pch=16, cex=0.6, xlim=c(0,10), ylim=c(0,500),col=c("steelblue2"))
m1 <- lm(x1[,3]~x1[,5])
eq <- paste('Y', '=', sprintf('%.3f',m1$coefficients[[1]]), '+', paste(sprintf('%.3f',m1$coefficients[[2]]), 'X', sep=''), sep=' ')
abline(m1, col='orange')
text(8,450, labels=eq,xpd=TRUE)
cor0 <- cor.test(x1[,3], x1[,5])
cl <- paste('R =', sprintf('%.2f',cor0$estimate), ',', 'P-value =', sprintf('%.3g',cor0$p.value))
text(8,480, labels=cl,xpd=TRUE)

plot(x1[,5],x1[,4],font.lab=1,ylab="Unique mPing Number", xlab="Excision Number", frame.plot=F, pch=16, cex=0.6, xlim=c(0,10), ylim=c(0,200),col=c("steelblue2"))
m1 <- lm(x1[,4]~x1[,5])
eq <- paste('Y', '=', sprintf('%.3f',m1$coefficients[[1]]), '+', paste(sprintf('%.3f',m1$coefficients[[2]]), 'X', sep=''), sep=' ')
abline(m1, col='orange')
text(8,150, labels=eq,xpd=TRUE)
cor0 <- cor.test(x1[,4], x1[,5])
cl <- paste('R =', sprintf('%.2f',cor0$estimate), ',', 'P-value =', sprintf('%.3g',cor0$p.value))
text(8,165, labels=cl,xpd=TRUE)

dev.off()

