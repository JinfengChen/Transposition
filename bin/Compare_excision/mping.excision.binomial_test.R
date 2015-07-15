#P=909/555/275/2=0.012 we expected only 50% of the RILs have these 555 mPing from HEG4.
#binom.test(5, 136, p=0.012)
#pvalue=0.036 when we observed 5 or more independent events for one locus

#update 20150628
#p=909/(555*(275/2)*2)=0.0059. 275/2 means only 50% of RILs have these 555 mPing from HEG4. *2 means diploid, which have two mping per loci.
#pvalue=0.001 when we observed 5 or more independent events for one locus
binom.test(5, 137, p=0.0059)
#pvalue=0.009 when we observed 5 or more independent events for one locus
binom.test(4, 137, p=0.0059)
