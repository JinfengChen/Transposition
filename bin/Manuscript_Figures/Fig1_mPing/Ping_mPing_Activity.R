#pdf("Ping_mPing_Activity.pdf")
x <- read.table("RIL230_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt", header=T)
parental <- x[,2]
hom <- x[,6]
som <- x[,7]+x[,8]
ping <- x[,9]
#fit hom~ping+parental
A <- data.frame(parental, ping, hom)
res <- lm(hom~ping+parental, data=A)
plot(res)
#fit som~ping+parental
B <- data.frame(parental, ping, som)
res <- lm(hom~ping+parental, data=B)
summary(res)
plot(res)
#dev.off()

