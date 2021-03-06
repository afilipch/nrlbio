'''Collection of function and classes specific for miRNA sequences'''

import sys;
import math;
from collections import defaultdict

from Bio.Seq import reverse_complement;
from Bio import SeqIO

from nrlbio.sequencetools import diverge_with_1mm, multifind

#constants
MODES_ORDER = ('m29a', 'm28a', 'm27a', 'm29', 'm28', 'm27', 'm38', 'mm28')


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
		
		self.m27 = get_seed_match(seq, 1, 7)
		self.m38 = get_seed_match(seq, 2, 8)
		self.m8 = reverse_complement(seq[7])
		self.m9 = reverse_complement(seq[8])
		self.first = 'A';
		
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
	
	def find_fixed_types(self, seq):
		mm28_count = 0
		for mm1 in diverge_with_1mm(get_seed_match(self.seq, 1, 8)):
			mm28_count += len(multifind(seq, mm1, overlap = False));
			
		
		m27_positions = multifind(seq, self.m27, overlap = False)
		m38_positions = multifind(seq, self.m38, overlap = False)
		mt = {'mm28': mm28_count, 'm38': len(m38_positions), 'm27': len(m27_positions), 'm28': 0, 'm29': 0, 'm27a': 0, 'm28a': 0, 'm29a': 0}
		
		for pos in m27_positions:
			if(len(seq)>pos+6 and seq[pos+6] == self.first):
				mt['m27a']+=1;
				first = True;
			else:
				first = False;
				
			if(pos>0 and seq[pos-1] == self.m8):
				mt['m28'] += 1;
				if(first):
					mt['m28a'] += 1;	
				if(pos>1 and seq[pos-2] == self.m9):
					mt['m29'] += 1;
					if(first):
						mt['m29a'] += 1;
						
		return mt;
	
	def find_fixed_positions(self, seq):
		mm28_positions = set();
		for mm1 in diverge_with_1mm(get_seed_match(self.seq, 1, 8)):
			mm28_positions.update(multifind(seq, mm1, overlap = False))
			
		
		m27_positions = multifind(seq, self.m27, overlap = False)
		m38_positions = multifind(seq, self.m38, overlap = False)
		mt = defaultdict(set)
		mt['mm28']= mm28_positions 
		mt['m27']= set(m27_positions)
		mt['m38']= set(m38_positions)
		
		
		for pos in m27_positions:
			if(len(seq)>pos+6 and seq[pos+6] == self.first):
				mt['m27a'].add(pos);
				first = True;
			else:
				first = False;
				
			if(pos>0 and seq[pos-1] == self.m8):
				mt['m28'].add(pos+1);
				if(first):
					mt['m28a'].add(pos+1);	
				if(pos>1 and seq[pos-2] == self.m9):
					mt['m29'].add(pos+1);
					if(first):
						mt['m29a'].add(pos+1);
						
		return dict([(x[0], list(sorted(list(x[1])))) for x in mt.items()] )
				
				
					
	def __str__(self):
		return "name=%s, seq=%s, start=%d, stop=%d, seed=%s, match=%s, expr=%1.1f" % (self.name, self.seq, self.seed_start, self.seed_stop, self.seed, self.match, self.expression);
		
		
		
def fasta2mirnas(fasta, seed_start=1, seed_stop=7):
	mirdict = {};
	for seqrecord in SeqIO.parse(fasta, "fasta"):
		mirdict[seqrecord.id] = Mirna(seqrecord.id, str(seqrecord.seq.upper()).replace('U', 'T'), seed_start, seed_stop)
	return mirdict;
	
	
def mirnas2families(mirnas):
	seed2mirna = defaultdict(list);
	for mirna in mirnas:
		seed2mirna[mirna.seed].append(mirna)
		
	return [Family(x) for x in seed2mirna.values()];
	
	
def find_family(mirid, families):
	for fam in families:
		if(mirid in fam.names):
			return fam
	else:
		return None;
	
	
	
def assign_expression(expr_file, mirdict, sep="\t"):
	with open(expr_file) as f:
		for l in f:
			mirid, expr = l.strip().split(sep)
			expr = float(expr);
			mirdict[mirid].expression = expr;
			
			
			
			
class Family():
	def __init__(self, mirnas):
		if(not mirnas):
			raise ValueError('It is impossimble to construct miRNA family from empty iterable of miRNAs')
		
		self.seed_start = mirnas[0].seed_start;
		self.seed_stop = mirnas[0].seed_stop;		
		self.seed = mirnas[0].seed
		self.match = mirnas[0].match;
		
		self.expression = sum([x.expression for x in mirnas])
		self.names = [x.name for x in mirnas]
		self.name = "|".join(self.names)
		
	#@property	
	#def name(self):
		#if(self.name):
			#return self.name
		#else:
			#return "|".join(self.names)
		
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
	
	def __str__(self):
		return "mirnas=%s, start=%d, stop=%d, seed=%s, match=%s, expr=%1.1f" % (self.name, self.seed_start, self.seed_stop, self.seed, self.match, self.expression);
	
	
	
	
#####################################################################################################################################################################
#Tools for miRNA destructive sites prediction


def mirfasta2conservation(mirfasta, translational_table=None):
	'''Parses conservation fasta file of miRNAs into a dictionary which connects conserved miRNAs
		
		mirfasta str: path to the conserved miRNAs, fasta file
		translational_table dict: dictionary which maps mirBase specie names to UCSC specie names
		
	Returns dict of dicts: key1 - miRNA id, key2 - specie id, value - miRNA seq 
	'''
	mir2cons = defaultdict(dict)
	refmir = {}
	records = SeqIO.parse(mirfasta, 'fasta');
	
	fsr = records.next();
	rspecie = fsr.id.split("-")[0]
	tkey = fsr.id
	refmir[tkey] = str(fsr.seq.upper()).replace('U', 'T')
	
	mir2cons = defaultdict(dict)
	for seqrecord in records:
		tspecie = seqrecord.id.split("-")[0]
		if(tspecie == rspecie):
			tkey = seqrecord.id
			refmir[tkey] = str(seqrecord.seq.upper()).replace('U', 'T')
		else:
			if(translational_table):
				mir2cons[tkey][translational_table[tspecie]] = str(seqrecord.seq.upper()).replace('U', 'T')
			else:
				mir2cons[tkey][tspecie] = str(seqrecord.seq.upper()).replace('U', 'T')
			
	return mir2cons, refmir



RTT = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
DAWARD = [1, 2, 2, 2, 2, 2, 2, 2, 1.5, 0, 0, 0, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
DPENALTY = [0, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]


def destructive_score(basepairing, award=DAWARD, penalty=DPENALTY):
		scores = [];
		bulge = [];
		pos = 0;
		
		revbase = [list(x)[::-1] for x in basepairing]
		for tu, tm, mm, mu in zip(*revbase):
			if(mm!=' ' and tm!=' '):
				if(RTT[mm] == tm):
					scores.append(award[pos]);
				else:
					scores.append(0.5*award[pos]);
					
				if(pos == 9 or pos == 10 or pos == 11):
					bulge.append(True);
					
			elif(tu!=' ' and mu!=' '):
				scores.append(-1.5*penalty[pos]);
				
			elif(tu!=' ' or mu!=' '):
				scores.append(-2*penalty[pos]);
				
			else:
				sys.exit('WARNING: Something wrong with rnahybrid parsing\n')
				
				
			if(mm!=' ' or mu!=' '):
				pos += 1;
				
			
		#fix last hanging nucleotide issue
		if(basepairing[0][0] != ' ' and basepairing[3][0] == ' '):
			scores[-1] = -1*scores[-1];
			
		#penalize slicing
		score = sum(scores);
		if(len(bulge)==3):
			score -= sum(award[1:7])
			
		return score;
	


def slicing_score(pattern):
	return int(all(pattern[1:20]))
		
		
		
		
#####################################################################################################################################################################
#Testing section
if(__name__ == "__main__"):
	mir = Mirna('hsa-miR-30a', 'TGTACTAGATCCTCGACTGGAAG', seed_start=1, seed_stop=7);
	fullmatch = 'TCTAGTACA'
	match = 'CTAGTACA'
	seq = match + "GTCACACGTG" + fullmatch + "GTCATCCCTTAAG" + match
	print mir.find_fixed_types(seq);
	
	
		