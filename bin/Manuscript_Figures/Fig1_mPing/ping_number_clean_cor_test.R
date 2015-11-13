x <- read.table("RIL230_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt", header=T)
#x <- read.table("RIL230_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.narrow_range.txt", header=T)
hom <- x[,6]
som <- x[,7]+x[,8]
ping <- x[,9]
cor.test(ping, hom)
cor.test(ping, som)
