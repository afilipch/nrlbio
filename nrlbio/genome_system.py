'''Collection of classes and functions to work with gene models'''

import sys
from collections import namedtuple, defaultdict, Counter
import bisect

from pybedtools import Interval, BedTool

from nrlbio.pybedtools_extension import get_distance


def exons2introns(exons):
	'''Returns list of introns on basis of given exons'''
	introns = [];
	if(exons and len(exons)>1):
		for i, (left,right) in enumerate(zip(exons[:-1], exons[1:])):
			introns.append(Interval(exons[0].chrom, left.end, right.start, '.', str(i), exons[0].strand));
	return introns
	


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
	'''Finds all interavls within a certain distance for given chrom, strand, start and stop positions:
	
	chrom str: chromosome of interval of interest
	strand str: strand of interval of interest
	start int: start(0-based inclusive) of interval of interest
	stop int: end(0-based exclusive) of interval of interest
	interval_dict dict:
		keys: tuple of chromosome strand
		values: list of intervals on corresponding chromosome and strand, supposed to be sorted;
	max_distance int: max allowed distance between interval of interest and others
	lookforward int: num of trials to look further for intersecting intervals
	
	Returns:
		start_intervals list of pybedtools.Interval: all interval within a max_distance to a start
		stop_intervals list of pybedtools.Interval: all interval within a max_distance to a stop
	'''	
	
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
		
		


class Transcript(object):
	'''Class wraps description and functionality for unprocessed RNA transcripts.
	
	Attributes:
		name str: name of the transcript
		transcript pybedtools.Interval: interval of the transcript itself
		exons list of pybedtools.Interval: exons of the transcript, supposed to be sorted
		introns list of pybedtools.Interval: introns of the transcript, supposed to be sorted
		cds pybedtools.Interval: interval of the transcript coding sequence
		utr5 pybedtools.Interval: interval of the transcript 5'UTR 
		utr3 pybedtools.Interval: interval of the transcript 3'UTR
		biotype str: type of the transcript;
		gene_name str: name of the parent gene
		otherfields dict: additonal miscellanious annotation of hte transcript
	'''
	def __init__(self, name, transcript, exons=[], introns=[], cds=None, utr5=None, utr3=None, biotype='', gene_name='', otherfields={}):
		self.name = name
		self.transcript = transcript;
		self.exons = exons;
		self.introns = introns;
		self.cds = cds 
		self.utr5 = utr5 
		self.utr3 = utr3 
		self.biotype = biotype
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
			
		introns = exons2introns(exons)
			
		return cls(name, transcript, exons=exons, introns=introns, cds=cds, utr5=utr5, utr3=utr3, biotype='refseq', gene_name=gene_name, otherfields=otherfields);
	
	
	
	def annotate_regulation(self, interval):
		'''annotates given interval according to its position on UTR5 or CDS or UTR3'''
		l = (interval.start - interval.end)/2
		if(self.utr3 and get_distance(self.utr3, interval)<l):
			return 'utr3'
		if(self.cds and get_distance(self.cds, interval)<l):
			return 'cds'
		if(self.utr5 and get_distance(self.utr5, interval)<l):
			return 'utr5'
		
		return None
		

	def annotate_transcription(self, interval):
		'''annotates given interval according to its position on UTR5 or CDS or UTR3'''
		l = (interval.start - interval.end)/2
		for exon in self.exons:
			if(exon and get_distance(exon, interval)<l):
				return 'exon'
		return 'intron'
	
	
	
	def to_refseq(self):
		exon_starts = ",".join([str(x.start) for x in self.exons]) + ","
		exon_stops = ",".join([str(x.stop) for x in self.exons]) + ","
		
		return "\t".join((  str(self.otherfields['bin']), self.name, self.transcript.chrom, self.transcript.strand, str(self.transcript.start), str(self.transcript.end),
				   str(self.cds.start), str(self.cds.end), str(self.otherfields['exons_count']),
				   exon_starts, exon_stops, self.transcript.score, str(self.gene_name)  ))
	
	def to_yaml(self):
		obj= {'name': self.name,
		'transcript': (self.transcript.chrom, self.transcript.start, self.transcript.stop, self.transcript.name, self.transcript.score, self.transcript.strand),
		'exons' : [(x.chrom, x.start, x.stop, '.', x.score, x.strand) for x in self.exons],
		'introns' : [(x.chrom, x.start, x.stop, '.', x.score, x.strand) for x in self.introns],
		'biotype' : self.biotype,
		'gene_name' : self.gene_name,
		'otherfields' : self.otherfields}
		
		if(self.cds):
			obj['cds'] = (self.cds.chrom, self.cds.start, self.cds.stop, '.', self.cds.score, self.cds.strand)
		if(self.utr5):
			obj['utr5'] = (self.utr5.chrom, self.utr5.start, self.utr5.stop, '.', self.utr5.score, self.utr5.strand)
		if(self.utr3):
			obj['utr3'] = (self.utr3.chrom, self.utr3.start, self.utr3.stop, '.', self.utr3.score, self.utr3.strand)
		
		return obj;
	
	@classmethod
	def from_yaml(cls, obj):
		exons = [Interval(*x) for x in obj['exons']]
		introns = [Interval(*x) for x in obj['introns']]
		
		utr3 = obj.get('utr3')
		if(utr3):
			utr3 = Interval(*utr3)
		utr5 = obj.get('utr5')
		if(utr5):
			utr5 = Interval(*utr5)
		cds = obj.get('cds')
		if(cds):
			cds = Interval(*cds)
			
		return cls(obj['name'], Interval(*obj['transcript']), exons=exons, introns=introns, cds=cds, utr5=utr5, utr3=utr3, biotype=obj['biotype'], gene_name=obj['gene_name'], otherfields=obj['otherfields']);	
		
		
		
	def __str__(self):
		return "Region: %s\nName: %s\nGene name: %s\n\nTranscript: %s\n5\'UTR: %s\nCDS: %s\n3\'UTR: %s\nexons:\n%sintrons:\n%s" % ("\t".join([ str(x) for x in (self.transcript.chrom, self.transcript.start, self.transcript.end, self.transcript.strand)]), self.name, self.gene_name, self.transcript, self.utr5, self.cds, self.utr3, "".join([str(x) for x in self.exons]), "".join([str(x) for x in self.introns]))
		
		
def generate_transcripts_from_refseq(refseq):
	'''Yields Transcript instances from given refseq(bed12 format) file'''
	with open(refseq) as f:
		for l in f:
			if(l.startswith('#')):
				pass;
			else:
				yield Transcript.from_refseq(l.strip().split("\t"));

		
	
class Gene(object):
	def __init__(self, name, gene, gene_symbol='', transcripts=[]):
		self.name=name
		self.gene_symbol = gene_symbol;
		self.gene = gene
		self.transcripts=transcripts
		
	def annotate_interval(self, interval):
		
		l = (interval.start - interval.end)/2	
		regulation = [];
		transcription = [];
		biotypes = [];
		for transcript in self.transcripts:
			if(get_distance(transcript.transcript, interval)<l):
				regulation.append(transcript.annotate_regulation(interval))
				transcription.append(transcript.annotate_transcription(interval))
				biotypes.append(transcript.biotype)
		return tuple(set(filter(bool, regulation))), tuple(set(filter(bool, transcription))), tuple(set(filter(bool, biotypes)))
		
		
	def to_yaml(self):
		return {'name': self.name, 'gene_symbol': self.gene_symbol, 'gene': (self.gene.chrom, self.gene.start, self.gene.stop, self.gene.name, self.gene.score, self.gene.strand), 'transcripts': [x.to_yaml() for x in self.transcripts]}
	
	@classmethod
	def from_yaml(cls, obj):
		return cls(obj['name'], Interval(*obj['gene']),  gene_symbol = obj['gene_symbol'], transcripts = [Transcript.from_yaml(x) for x in obj['transcripts']])
		
	def __str__(self):
		l = ['Region: %s\nName: %s\nGene_symbol: %s\n' % ("\t".join([ str(x) for x in (self.gene.chrom, self.gene.start, self.gene.end, self.gene.strand)]), self.name, self.gene_symbol), "Gene: %s" % str(self.gene)]
		for transcript in self.transcripts:
			l.append("%s\n%s" % ("_"*140, str(transcript)))
		return "\n".join(l)
	
	
		
		

def _gff3_to_transcript(transcript_interval, transcription_blocks):
	'''creates Transcript object on basis of gff3 specific transcript_interval and transcription_blocks entries. Is used internally in gff3_to_genes'''

	exons = [Interval(x.chrom, x.start, x.stop, x.name, x.attrs['rank'], x.strand) for x in transcription_blocks.get('exon')]
	introns = exons2introns(exons)
	
	utr3_intervals = transcription_blocks.get('three_prime_UTR');
	cds_intervals = transcription_blocks.get('CDS');
	utr5_intervals = transcription_blocks.get('five_prime_UTR');
	if(utr3_intervals):
		utr3=Interval(transcript_interval.chrom, min([x.start for x in utr3_intervals]), max([x.end for x in utr3_intervals]), '.', '0', transcript_interval.strand)
	else:
		utr3=None;
	if(cds_intervals):
		cds=Interval(transcript_interval.chrom, min([x.start for x in cds_intervals]), max([x.end for x in cds_intervals]), '.', '0', transcript_interval.strand)
	else:
		cds=None
	if(utr5_intervals):
		utr5=Interval(transcript_interval.chrom, min([x.start for x in utr5_intervals]), max([x.end for x in utr5_intervals]), '.', '0', transcript_interval.strand)
	else:
		utr5=None
		
	return Transcript(transcript_interval.name.split(":")[1], transcript_interval, exons=exons, introns=introns, cds=cds, utr5=utr5, utr3=utr3, biotype=transcript_interval.attrs['biotype'], gene_name=transcript_interval.attrs['Parent'].split(":")[1])
				
				
def gff3_to_genes(gff3):
	'''Yields genes from given iterable of gff3 intervals'''
	genes = [];
	gene_interval = None
	transcript_interval=None
	transcripts = [];
	transcription_blocks = defaultdict(list);
	#gt = []
	
	for interval in gff3:
		if(interval.attrs.get('ID', 'stub').startswith('gene')):
			if(transcript_interval):
				transcripts.append(_gff3_to_transcript(transcript_interval, transcription_blocks));
			if(gene_interval):
				gene = Gene(gene_interval.name.split(":")[1], gene_interval, gene_interval.attrs.get('Name', ''), transcripts);
				genes.append(gene);
			gene_interval = interval
			transcript_interval = None;
			transcripts = [];
			
		elif(interval.attrs.get('ID', 'stub').startswith('transcript')):
			if(transcript_interval):
				transcripts.append(_gff3_to_transcript(transcript_interval, transcription_blocks))
			transcript_interval = interval;
			transcription_blocks = defaultdict(list);
		else:
			transcription_blocks[interval[2]].append(interval);
			
	else:
		if(transcript_interval):
			transcripts.append(_gff3_to_transcript(transcript_interval, transcription_blocks));
		if(gene_interval):
			gene = Gene(gene_interval.name.split(":")[1], gene_interval, gene_interval.attrs.get('Name', ''), transcripts);
			genes.append(gene);
		
	return genes

def yaml2genes(path):
	import yaml
	with open(path, 'r') as f:
		return [Gene.from_yaml(x) for x in yaml.load(f)]
	
	
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
		
	#for gene in generate_genes_from_refseq(sys.argv[1]):
		#pass
		
	import yaml	
	#with open('todel.yaml', 'w') as f:
		#f.write(yaml.dump(gff3_to_genes(BedTool(sys.argv[1])), default_flow_style=False))	
		
	with open('todel.yml', 'w') as f:
		for c, gene in enumerate(gff3_to_genes(BedTool(sys.argv[1]))):
			print "*"*140
			print gene
			print
			if(c==1):
				f.write(yaml.dump(gene.to_yaml(), default_flow_style=False))
				break;
				
	#with open(sys.argv[2], 'r') as f:
		#d = yaml.load(f)
		#print Gene.from_yaml(d);