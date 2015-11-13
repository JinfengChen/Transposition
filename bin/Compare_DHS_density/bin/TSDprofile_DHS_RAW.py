import sys
import re
import numpy as np
import HTSeq
import numpy
from collections import defaultdict

def sum_table(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                for i in range(1,len(unit)):
                    data[i].append(int(unit[i]))
    for r in sorted(data.keys(), key=int):
        print '%s\t%s\t%s' %(r, np.mean(data[r]), np.std(data[r]))
    return data


if len(sys.argv) < 3:
    HELP='''
python TSDprofile_DHS.py test.bed test.profile > test.sum &
    '''



#bamfile = HTSeq.BAM_Reader( "DHS.Chr1.unique.bam" )
sortedbamfile = HTSeq.BAM_Reader( "../input/DHS.unique.bam" )
gfffile = sys.argv[1]
tsd_profile = sys.argv[2]

halfwinwidth = 1500
fragmentsize = 120
#total = 60745783.00/1000000 ## nucleosome
total = 23299296/1000000 #DHS unique
gsize = 372000000
readlen = 36

##gff
#Chr1    RelocaTE        mPing   10046640        10046642        .       .       .       Strains=HEG4_2;GT=heterozygous
#Chr1    RelocaTE        mPing   10508960        10508962        .       .       .       Strains=RIL;GT=heterozygous
tsdpos = []
with open (gfffile, 'r') as filehd:
    for line in filehd:
        line = line.rstrip()
        if line.startswith('Chr'): 
            unit = re.split(r'\t',line)
            #print unit[0], unit[3], unit[4]
            tsdpos.append([unit[0], unit[3], unit[4]])

ofile = open(sys.argv[2], 'w')
ofile1= open('test_all.sum', 'w')
profile_ALL = numpy.zeros( 2*halfwinwidth, dtype='i' )
for p in tsdpos:
   #define empty vector of profile
   profile = numpy.zeros( 2*halfwinwidth, dtype='i' )
   #define 1000 bp window around TSD
   window = HTSeq.GenomicInterval( p[0], int(p[1]) - halfwinwidth, int(p[1]) + halfwinwidth, "." )
   #for all the aligned reads in window
   for almnt in sortedbamfile[ window ]:
       #set reads to fragmentsize, which might be 200 or 300
       almnt.iv.length = fragmentsize
       #convert reads coordinate on chromosome to coordinate in window
       start_in_window = almnt.iv.start - int(p[1]) + halfwinwidth
       end_in_window   = almnt.iv.end   - int(p[1]) + halfwinwidth
       #make sure coordinate within the window range becaue we extended reads to 200 bp
       start_in_window = max( start_in_window, 0 )
       end_in_window = min( end_in_window, 2*halfwinwidth )
       if start_in_window >= 2*halfwinwidth or end_in_window < 0:
           continue
       #add depth to profile
       profile[start_in_window : end_in_window] += 1
       profile_ALL[start_in_window : end_in_window] += 1
       
   mping = '%s_%s' %(p[0], p[1])
   depth = '\t'.join(map(str, profile))
   print >> ofile, '%s\t%s' %(mping, depth)

for i in range(len(profile_ALL)):
    print >> ofile1, '%s\t%s\t%s' %(i, profile_ALL[i], float(profile_ALL[i])/total)
ofile.close()
ofile1.close()

sum_table(sys.argv[2])

#print profile        
#for i in profile:
#    print i

#pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile )  
#pyplot.savefig('Nucleosome_HighExp_Chr1.pdf')
#pyplot.savefig('%s.pdf' %(sys.argv[2]))

