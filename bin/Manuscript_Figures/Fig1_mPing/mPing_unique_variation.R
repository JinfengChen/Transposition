pdf("unique.mping.range.plot.pdf")
x <- read.table("RIL230_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt", header=T)
hist(x[,5], breaks=100, xlab="Unique mPing number", ylab="Frequency in 230 RILs", main="Variation of unique mPing in 230 RILs")
dev.off()

