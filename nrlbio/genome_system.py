'''Collection of classes and functions to work with gene models'''

import sys
from collections import namedtuple, defaultdict
import bisect

import pybedtools






def refseq2interval(refseq):
	for l in refseq:
		a = l.split("\t");
		name = a[1];
		chrom = a[2]
		strand = a[3];
		for i, (s,e) in enumerate(zip(a[9].split(",")[:-1], a[10].split(",")[:-1])):
			yield pybedtools.Interval(chrom, int(s), int(e), name, str(i), strand);
			

def intervals2dict(intervals):
	d = defaultdict(list);
	for i in intervals:
		d[(i.chrom, i.strand)].append(i);
	for v in d.values():
		v.sort(key= lambda x: x.start);
	return d;
	
	
def find_closest(chrom, strand, start, stop, interval_dict, max_distance=4, lookforward = 10):
	l = interval_dict[(chrom, strand)]
	index =  bisect.bisect([x.start for x in l], start)
	
	start_intervals = [];
	for i in range(index, len(l)):
		if(abs(start - l[i].start) <= max_distance):
			start_intervals.append(l[i]);
		else:
			break;
			
	for i in range(index-1, -1, -1):
		if(abs(start - l[i].start) <= max_distance):
			start_intervals.append(l[i]);
		else:
			break;
			
	stop_intervals = [];
	for i in range(index, min(index+lookforward, len(l))):
		if(abs(stop - l[i].stop) <= max_distance):
			stop_intervals.append(l[i]);
		#else:
			#break;
			
	for i in range(index-1, max(-1, index-lookforward-1), -1):
		if(abs(stop - l[i].stop) <= max_distance):
			stop_intervals.append(l[i]);
		#else:
			#break;
			
	return start_intervals, stop_intervals
		
		
	
			




class Gene(object):
	def __init__(self, name, chrom, strand, start, end, score=0, exons=[]):
		self.name = name;
		self.chrom = chrom;
		self.strand = strand;
		self.start = start;
		self.end = end;
		self.score = score;
		self.exons = exons;
		
	@classmethod
	def from_refseq(cls, l):
		a = l.split("\t");
		name = a[1];
		chrom = a[2]
		strand = a[3];
		start = int(a[4]);
		end = int(a[5]);
		score = int(a[8])
		
		exons = [];
		for i, (s,e) in enumerate(zip(a[9].split(","), a[10].split(","))):
			exons.append(pybedtools.Interval(chrom, int(s), int(e), name, str(i), strand));
			
		return cls(name, chrom, strand, start, end, score, exons);
		
	def __str__(self):
		exon_starts = ",".join([str(x.start) for x in self.exons])
		exon_stops = ",".join([str(x.stop) for x in self.exons])
		return "\t".join(self.chrom, str(self.start), str(self.end), self.name, str(self.score), self.strand, self.exon_starts, self.exon_stops)
		
	
	
	
	
	
#testing section
if(__name__ == "__main__"):
	f = open(sys.argv[1]);
	f.readline();
	#for i in refseq2interval(f):
		#print i;
	d = intervals2dict(refseq2interval(f))
	for add in range(1):
		start_intervals, stop_intervals = find_closest("chr3", "+", 1134344+add*500, 1134810+add*500, d, max_distance = 4)
	
	print "\nstartintervals\n"
	for i in start_intervals:
		print i;
	print "\nstopintervals\n"
	for i in stop_intervals:
		print i;		