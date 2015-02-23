# /usr/bin/python
'''extends functionality of pybedtools'''
import sys
from collections import defaultdict

from pybedtools import Interval, BedTool, create_interval_from_list

from nrlbio.numerictools import overlap

def generate_overlaping_intervals(bed, distance):
	'''applicable only for sorted(with strandness) bed files'''
	
	first = bed[0];
		
	merged = [first]
	start, end = first.start, first.stop
	rname = (first.chrom, first.strand)
	
	for i in bed[1:]:
		if(rname == (i.chrom, i.strand)):
			s, e = overlap((i.start, i.stop), (start, end))
			if(e - s >= -1*distance):
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
	
	
def doublebed2dict(bed):
	'''Creates dictionary from doublebed file. keys - ids, values: list of interacting regions
	
		bed pybedtools.BedTool: genomic intervals.
		
		Returns dict: keys - interation ids, values: list of interacting regions
	'''
	d = defaultdict(list)
	for i in bed:
		name = i.name.split("|")[0]
		d[name].append(i);
	return d;	
	
	
def construct_gff_interval(chrom, start, stop, feature, score='0', strand='.', source='.', frame='.', attrs=[]):
	attrs_str = "; ".join(["%s=%s" % (str(x[0]), str(x[1])) for x in attrs])
	return create_interval_from_list( [chrom, source, feature, str(start), str(stop), score, strand, frame, attrs_str] );
	
def list2interval(l):
	'''creates an interval from the list of strings. items in the list should be in the order: chromosome, start, stop, name, score, strand. In case of dummy bed entry(chromosome is equal ".") returns None
	'''
	if(l[0] == '.'):
		return None
	else:	
		return Interval(l[0], int(l[1]), int(l[2]), name=l[3], score=l[4], strand=l[5])	
		
		
		
def interval2seq(interval, reference):
	if(interval.strand == '+'):
		return str(reference[interval.chrom][interval.start:interval.stop].seq.upper())
	elif(interval.strand == '-'):
		return str(reference[interval.chrom][interval.start:interval.stop].seq.reverse_complement().upper())
	else:
		sys.stderr.write("Strand is not defined. Plus strand sequence will be returned\n")
		return str(reference[interval.chrom][interval.start:interval.stop].seq.upper())		
	
	
#testing section
if(__name__ == "__main__"):
	
	bed = BedTool(sys.argv[1])
	
	for m in generate_overlaping_intervals(bed, 1):
		for i in m:
			print i[:6]
			print i[6:]
			ai = Interval(i[6], int(i[7]), int(i[8]), i[9], i[10], i[11], otherfields = i[:6] + i[12:])
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