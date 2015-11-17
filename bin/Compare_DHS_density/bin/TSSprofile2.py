import sys
import HTSeq
import numpy
import matplotlib as mpl
mpl.use('pdf')
from matplotlib import pyplot

#bamfile = HTSeq.BAM_Reader( "Chr1.unique.bam" )
#bamfile = HTSeq.BAM_Reader( "../input/Nucleosome.unique.bam" )
bamfile = HTSeq.BAM_Reader( "../input/DHS.Chr1.unique.bam" )
#gtffile = HTSeq.GFF_Reader( "../input/MSU7.gene.exon_number.HighExp.gtf" )
gtffile = HTSeq.GFF_Reader(sys.argv[1])

halfwinwidth = 2000
fragmentsize = 150
total = 60745783.00/1000000 ## nucleosome
gsize = 372000000

coverage = HTSeq.GenomicArray( "auto", stranded=False, typecode="i" )
for almnt in bamfile:
   if almnt.aligned:
      #almnt.iv.length = fragmentsize
      #print almnt.iv
      if not almnt.iv.start < 500:
          coverage[ almnt.iv ] += 1

tsspos = set()
for feature in gtffile:
   if feature.type == "exon" and feature.attr["exon_number"] == "-1":
      #print feature.iv.start_d_as_pos.pos
      if feature.iv.start_d_as_pos.pos > 5000:
          tsspos.add( feature.iv.end_d_as_pos )

exon_s = set()
exon_e = set()
for feature in gtffile:
   if feature.type == "exon":
       if feature.iv.start_d_as_pos.pos > 5000:
           if feature.iv.strand == '+':
               #print feature.iv.start_d_as_pos, feature.iv.end_d_as_pos
               exon_s.add(feature.iv.start_d_as_pos)
               exon_e.add(feature.iv.end_d_as_pos)
           else:
               exon_s.add(feature.iv.end_d_as_pos)
               exon_e.add(feature.iv.start_d_as_pos)

profile = numpy.zeros( 2*halfwinwidth, dtype='i' )      
for p in tsspos:
   window = HTSeq.GenomicInterval( p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, "." )
   wincvg = numpy.fromiter( coverage[window], dtype='i', count=2*halfwinwidth )
   if p.strand == "+":
      profile += wincvg
   else:
      profile += wincvg[::-1]

ofile = open('DHS_TTS1.profile', 'w')
for i in range(len(profile)):
    print >> ofile, profile[i]
ofile.close()

pyplot.plot( numpy.arange( -halfwinwidth, halfwinwidth ), profile )  
#pyplot.savefig('Nucleosome_HighExp_Chr1.pdf')
pyplot.savefig('DHS_TTS1.pdf')

