x <- read.table("../input/RIL275_RelocaTE.sofia.ping_code.table", header=T)
y <- read.table("mping.excision.draw.highexcision.ping_code.list")
t.test(x[,1], y[,1])
#pvalue=3.66e-06

