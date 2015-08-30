data =read.table("mPing_boundary.linked_50Mb_debug2.distance_accumulation_excision.list", header=T)
excision = data[,2]
control  = data[,5]
ks.test(excision, control)
