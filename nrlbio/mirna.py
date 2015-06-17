'''Collection of function and classes specific for miRNA sequences'''

import sys;
from collections import defaultdict

from Bio.Seq import reverse_complement;
from Bio import SeqIO

from nrlbio.sequencetools import diverge_with_1mm, multifind

def mirbase_conversion(vers1, vers2):
	conv = {};
	d1 = {};
	d2 = {};
	
	for seqrecord in SeqIO.parse(vers1, "fasta"):
		d1[str(seqrecord.seq.upper())] = seqrecord.id;
	for seqrecord in SeqIO.parse(vers2, "fasta"):
		d2[str(seqrecord.seq.upper())] = seqrecord.id;
		
	for seq, mirid in d1.items():
		conv[mirid] = d2.get(seq, None)
		
	return conv	



def get_seed_match(seq, start, stop):
	return reverse_complement(seq[start: stop])

			
class Mirna(object):
	'''Wraps functionality, specific for miRNA sequence
	
	Attributes:
		name str: miRNA id, For example 'cel-let-7'
		seq str: miRNA sequence
		seed_start int: start postion of mirna seed, 0-based and inclusive
		seed_stop int: stop postion of mirna seed, 0-based and exclusive
		seed str: miRNA's seed sequence
		match str: miRNA's seed match sequence
	'''
	def __init__(self, name, seq, seed_start, seed_stop):
		self.name = name;
		self.seq = seq;
		self.seed_start = seed_start
		self.seed_stop = seed_stop;
		
		self.expression = 0;
		self.seed = seq[seed_start: seed_stop];
		self.match = get_seed_match(seq, seed_start, seed_stop);
		
	def set_1mm_matches(self):
		self.matches_1mm = diverge_with_1mm(self.match);
		
		
	def find_matches(self, seq, overlap=False):
		return multifind(seq, self.match, overlap = overlap)
	
	
	def find_1mm_matches(self, seq, overlap=False):
		if (not hasattr(self, 'matches_1mm')):
			self.set_1mm_matches();		
		dm = {};
		for mm1 in self.matches_1mm:
			l = multifind(seq, mm1, overlap = overlap);
			if(l):
				dm[mm1] = l;
		return dm;
		
		
		
def fasta2mirnas(fasta, seed_start, seed_stop):
	mirdict = {};
	for seqrecord in SeqIO.parse(fasta, "fasta"):
		mirdict[seqrecord.id] = Mirna(seqrecord.id, str(seqrecord.seq.upper()), seed_start, seed_stop)
	return mirdict;
	
	
def mirnas2families(mirnas):
	seed2mirna = defaultdict(list);
	for mirna in mirnas:
		seed2mirna[mirna.seed].append(mirna)
		
	return [Family(x) for x in seed2mirna.values()];
	
	
	
def assign_expression(expr_file, mirdict, sep="\t"):
	with open(expr_file) as f:
		for l in f:
			mirid, expr = l.strip().split(sep)
			expr = float(expr);
			mirdict[mirid].expression = expr;
			
			
			
			
class Family():
	def __init__(self, mirnas, name=''):
		if(not mirnas):
			raise ValueError('It is impossimble to construct miRNA family from empty iterable of miRNAs')
		if(name):
			self.name = name;
		else:
			self.name = "|".join([x.name for x in mirnas])
		
		self.seed_start = mirnas[0].seed_start;
		self.seed_stop = mirnas[0].seed_stop;		
		self.seed = mirnas[0].seed
		self.match = mirnas[0].match;
		
		self.expression = sum([x.expression for x in mirnas])
		
	def set_1mm_matches(self):
		self.matches_1mm = diverge_with_1mm(self.match);
		
		
	def find_matches(self, seq, overlap=False):
		return multifind(seq, self.match, overlap = overlap)
	
	
	def find_1mm_matches(self, seq, overlap=False):
		if (not hasattr(self, 'matches_1mm')):
			self.set_1mm_matches();
			
		dm = {};
		for mm1 in self.matches_1mm:
			l = multifind(seq, mm1, overlap = overlap);
			if(l):
				dm[mm1] = l;
		return dm;
		
		
if(__name__ == "__main__"):
	#mir = Mirna('hsa-miR-30a', 'TGTAAACATCCTCGACTGGAAG', seed_start=1, seed_stop=7);
	#mir.set_1mm_matches();
	#print 'GTTTAC'
	#print mir.match;
	#print
	#print mir.find_1mm_matches('GTTTACAGCTGTATACAAGTTTAA')
	
	#for mirid, mirna in fasta2mirnas(sys.argv[1], seed_start=1, seed_stop=8).items():
		#print mirid
		#print mirna.seq;
		#print mirna.seed;
		#print mirna.match;
		#print "_"*160	
		
	#mirdict = fasta2mirnas(sys.argv[1], seed_start=1, seed_stop=7);
	#assign_expression(sys.argv[2], mirdict, sep="\t");
	
	#conv = mirbase_conversion(sys.argv[1], sys.argv[2])
	
	mirdict = fasta2mirnas(sys.argv[1], seed_start=1, seed_stop=7);
	assign_expression(sys.argv[2], mirdict, sep="\t");	
	families = mirnas2families(mirdict.values())
	for fam in families:
		print fam.name
		print fam.seed
		print "_"*160
	
	#print mirdict['hsa-miR-326'].props
	
	#for v in mirdict.values():
		#print v.props
	
	#for kv in conv.items():
		#print "%s\t%s" % kv
		