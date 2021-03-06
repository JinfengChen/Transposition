pdf("Plot_Density_Chr4.pdf", width = 10, height = 4)
exon <- read.table("../input/Chromosome4/Chr4.exon.histogram.txt")
dhs  <- read.table("../input/Chromosome4/Chr4.DHS.histogram.txt")
ril  <- read.table("../input/Chromosome4/Chr4.RIL.histogram.txt")
sim  <- read.table("../input/Chromosome4/Chr4.simulation.histogram.txt")

plot(exon[,2]/1000000, exon[,4], col='blue', type='l', main="", xlab="Chromosomal position (Mb)",  ylab="Density (200 kb)", xlim=c(0, 40), ylim=c(0, 30), xaxt='n', frame.plot=FALSE)
lines(dhs[,2]/1000000, dhs[,4]*0.25, col='red', type='l')
lines(ril[,2]/1000000, ril[,4], col='aquamarine3', type='l')
lines(sim[,2]/1000000, sim[,4], col='dim gray', type='l')
axis(1,seq(0, 35, by=5),line=0, labels=seq(0, 35, by=5))
#text(seq(1:length(dist[,1])),rep(-0.04,7), cex=1, offset=2,labels=rev(dist[,2]/1000),srt=0,xpd=TRUE)
legend(2, 32, bty='n', border='NA', lty= c(1,1,1,1), cex=1 , lwd = 1.5, col=c("blue", "red", "aquamarine3", "dim gray"), c("Exon","DHS","Unique", "Control"))


dev.off()
