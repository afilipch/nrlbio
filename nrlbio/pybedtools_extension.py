# /usr/bin/python
'''extends functionality of pybedtools'''
import sys

import pybedtools

from nrlbio.numerictools import overlap

def generate_overlaping_intervals(bed, distance):
	'''applicable only for sorted bed(with strandness) files'''
	
	first = bed[0];
		
	merged = [first]
	start, end = first.start, first.stop
	rname = (first.chrom, first.strand)
	
	for i in bed[1:]:
		if(rname == (i.chrom, i.strand)):
			s, e = overlap((i.start, i.stop), (start, end))
			if(e - s >= distance):
				merged.append(i);
				end = max(i.stop, end)
			else:
				yield merged;
				start, end = i.start, i.stop;
				merged = [i];
		else:
			yield merged
			start, end = i.start, i.stop
			merged = [i];
			rname = (i.chrom, i.strand);
	yield merged
	
def construct_gff_interval(chrom, start, stop, feature, score='0', strand='.', source='.', frame='.', attrs={}):
	attrs_str = "; ".join( [ "%s=%s" % (str(x[0]), str(x[1])) for x in attrs.items() ] )
	return pybedtools.create_interval_from_list( [chrom, source, feature, str(start), str(stop), score, strand, frame, attrs_str] );
	
	
#testing section
if(__name__ == "__main__"):
	
	bed = pybedtools.BedTool(sys.argv[1])
	
	for m in generate_overlaping_intervals(bed, 1):
		for i in m:
			print i[:6]
			print i[6:]
			ai = pybedtools.Interval(i[6], int(i[7]), int(i[8]), i[9], i[10], i[11], otherfields = i[:6] + i[12:])
			print ai[:6]
			print ai[6:]
			print
		print
		#print
		
	#i = construct_gff_interval('chr1', 10, 100, 'gene', strand='+', attrs = {'type': 'exon', 'Gene': 'Sox2'})
	#b = pybedtools.BedTool([i]);
	#b.saveas("todel.gff")
	#a = pybedtools.BedTool("todel.gff");
	#for i in a:
		#print i.attrs
	#print i.file_type