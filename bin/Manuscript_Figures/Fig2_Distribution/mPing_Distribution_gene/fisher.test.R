#             strain   Somatic
#-1501-2000   70       144
#-1001-1500   55       189
#p-value=0.015
test = rbind(c(70, 55), c(144, 189))
chisq.test(test)
fisher.test(test)
