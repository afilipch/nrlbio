'''Collection of classes and functions to work with gene models'''

import sys
from collections import namedtuple, defaultdict
import bisect

from pybedtools import Interval, BedTool



def seqrecord2seq(seqrecord, start, end, strand=None):
	'''Extracts sequence from Bio.SeqRecord for a given strand and position.
	
		seqrecord Bio.SeqRecord: sequence record to extract sequence from
		start int: start(0-based inclusive) on the seqrecord  the sequence extracted
		end int: end(0-based exclusive) on the seqrecord  the sequence extracted
		strand '+'|'-'|None: strand of the sequence. If it is set to '-', then reverse conmplement of the sequence will be returned
	'''
		
	if(strand == '+'):
		return str(seqrecord[start:end].seq.upper())
	elif(strand == '-'):
		return str(seqrecord[start:end].seq.reverse_complement().upper())
	else:
		sys.stderr.write("WARNING: Strand is not properly defined. Plus strand sequence will be returned\n")
		return str(seqrecord[start:end].seq.reverse_complement().upper())
	
	


def refseq2interval(refseq):
	for l in refseq:
		a = l.split("\t");
		name = a[1];
		chrom = a[2]
		strand = a[3];
		for i, (s,e) in enumerate(zip(a[9].split(",")[:-1], a[10].split(",")[:-1])):
			yield Interval(chrom, int(s), int(e), name, str(i), strand);
			

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
	def __init__(self, name, transcript, exons=[], introns=[], cds=None, utr5=None, utr3=None, type_='', gene_name='', otherfields={}):
		self.name = name
		self.transcript = transcript;
		self.exons = exons;
		self.introns = introns;
		self.cds = cds 
		self.utr5 = utr5 
		self.utr3 = utr3 
		self.type = type_
		self.gene_name = gene_name
		self.otherfields = otherfields
		
	@classmethod
	def from_refseq(cls, a):
		#a = l.split("\t");
		otherfields = {'bin': int(a[0]), 'exons_count': int(a[8])}
		name = a[1];
		chrom = a[2]
		strand = a[3];
		start = int(a[4]);
		end = int(a[5]);
		cds_start = int(a[6])
		cds_end = int(a[7])
		score = a[11]
		gene_name = a[12]
			
		transcript = Interval(chrom, start, end, name, score, strand)
		cds = Interval(chrom, cds_start, cds_end, name, score, strand)
		
		if(strand == '+'):
			utr5 = Interval(chrom, start, cds_start, name, score, strand);
			utr3 = Interval(chrom, cds_end, end, name, score, strand);
		elif(strand == '-'):
			utr3 = Interval(chrom, start, cds_start, name, score, strand);
			utr5 = Interval(chrom, cds_end, end, name, score, strand);
		else:
			utr3 = None
			utr5 = None
		
		exons = [];
		for i, (s,e) in enumerate(zip(a[9].split(","), a[10].split(","))):
			if(s and e):
				exons.append(Interval(chrom, int(s), int(e), name, str(i), strand));
			
		introns = [];
		for i, (left,right) in enumerate(zip(exons[:-1], exons[1:])):
			introns.append(Interval(chrom, left.end, right.start, name, str(i), strand));
			
		return cls(name, transcript, exons=exons, introns=introns, cds=cds, utr5=utr5, utr3=utr3, type_='refseq', gene_name=gene_name, otherfields=otherfields);
		
	
	def to_refseq(self):
		exon_starts = ",".join([str(x.start) for x in self.exons]) + ","
		exon_stops = ",".join([str(x.stop) for x in self.exons]) + ","
		
		return "\t".join((  str(self.otherfields['bin']), self.name, self.transcript.chrom, self.transcript.strand, str(self.transcript.start), str(self.transcript.end),
				   str(self.cds.start), str(self.cds.end), str(self.otherfields['exons_count']),
				   exon_starts, exon_stops, self.transcript.score, str(self.gene_name)  ))
		
		
	def __str__(self):
		return self.to_refseq()
		
		
def generate_genes_from_refseq(refseq):
	'''Yields Gene instances from given refseq(bed12 format) file'''
	with open(refseq) as f:
		for l in f:
			if(l.startswith('#')):
				pass;
			else:
				yield Gene.from_refseq(l.strip().split("\t"));

		
	
	
	
	
#testing section
if(__name__ == "__main__"):
	#f = open(sys.argv[1]);
	#f.readline();
	##for i in refseq2interval(f):
		##print i;
	#d = intervals2dict(refseq2interval(f))
	#for add in range(1):
		#start_intervals, stop_intervals = find_closest("chr3", "+", 1134344+add*500, 1134810+add*500, d, max_distance = 4)
	
	#print "\nstartintervals\n"
	#for i in start_intervals:
		#print i;
	#print "\nstopintervals\n"
	#for i in stop_intervals:
		#print i;
		
	for gene in generate_genes_from_refseq(sys.argv[1]):
		pass