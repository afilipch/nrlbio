'''Collection of function and classes specific for miRNA sequences'''
from Bio.Seq import reverse_complement;
from itertools import combinations, product, combinations_with_replacement;

from nrlbio.sequencetools import diverge_with_1mm

def get_seed_match(seq, start, stop):
	return reverse_complement(seq[start: stop])
	

class Mirna(object):
	'''Wraps functionality specific for miRNA sequence
	
	Attributes:
		name str: miRNA id, For example 'cel-let-7'
		seq str: miRNA sequence
		seed_start int: start postion of mirna seed, 0-based and inclusive
		seed_stop int: stop postion of mirna seed, 0-based and exclusive
		seed str: miRNA's seed sequence
		match str: miRNA's seed match sequence
	'''
	def __init__(self, name, seq, seed_start=1, seed_stop=7):
		self.name = name;
		self.seq = seq;
		self.seed_start = seed_start
		self.seed_stop = seed_stop;
		self.seed = seq[seed_start: seed_stop];
		self.match = get_seed_match(seq, seed_start, seed_stop);
		
	def set_1mm_matches(self):
		self.matches_1mm = diverge_with_1mm(self.match);
		