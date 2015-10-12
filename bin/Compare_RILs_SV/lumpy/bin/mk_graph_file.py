import sys
#python mk_graph_file.py sample.1.fq.readdepth

# make a bedgraph file
def mk_graph_file(out_fn):
	bedgraph_fn = out_fn + ".bed"
	f = open(out_fn + '.txt', 'r')
	fdata = f.readlines()
	f.close()
	
	nf = []
	prev_chr = ""
	prev_end = 0
	for x in fdata:
  		pieces = x.split()
  		idx = pieces[1].find(':')
  		idx2 = pieces[1].find('-')
  		chr = pieces[1][0:idx]
		# remove any 'chr' prefixes from the chrom name
		if chr.startswith('chr'): chr = chr[3:]
  		start = int(pieces[1][idx+1:idx2])
  		end = int(pieces[1][idx2+1:])
  		prev_chr, prev_end = chr, end
  		line = "%s\t%d\t%d\t%s" % (chr, start, end, pieces[3])
  		nf.append(line)
	nfstr = '\n'.join(nf)
	
	f = open(bedgraph_fn, 'w')
	f.write(nfstr)
	f.close()
# end of make bedgraph file

mk_graph_file(sys.argv[1])

