'''library contains classes and functions to deal with chimeric reads'''

import sys
from copy import copy
from itertools import combinations;


class ChimeraException(Exception):
    pass

    
def _alignment_score(chimera):
	'''score function for Chimera based only on alignment score'''
	return reduce(lambda x, y: x*y, chimera.AS)
	
def _alignment_gap_score(chimera, maxgap=8):
	'''score function for Chimera based on alignment score and gap'''
	return chimera.AS[0]*chimera.AS[1]*(1 - chimera.gap**2/maxgap**2)	
	
	

class Chimera(object):
	'''Chimera represents chimeric read based on two consecutive(with some gap or overlap) mappings for one read.
	
		aligned_reads list: element is pysam.AlignedRead. Two consecutive hits comprising the chimera. Must be in consecutive(from left to right) order
		rnames list: list of hits' reference names. The order must correspond to the order of aligned_reads list
		AS list: list of alignment scores corresponding to the aligned_reads
		self.gap int: gap(if negative, overlap) between hits(aligned_reads)
		self.score float: Chimera score. The bigger the score the more reliable is chimera
	'''	
		
		
	def __init__(self, aligned_reads, samfile, score_function):
		if(len(set([x.qname for x in aligned_reads])) != 1): 
			raise ChimeraException('chimera cannot be made from aligned reads with differenr identifiers')
			return None
			
		self.aligned_reads = copy(aligned_reads);
		self.rnames = [samfile.getrname(x.tid) for x in aligned_reads];
		
		self.AS = [x.opt('AS') for x in aligned_reads]
		self.gap = aligned_reads[1].qstart - aligned_reads[0].qend
		self.score = score_function(self);
		
		
	def __str__(self):
		'''converts Chimera in a bed-like entry(line). First six elements correspond to the first hit, Second six elements to the second one'''
		l = [];
		for a, rname in zip(self.aligned_reads, self.rnames):
			l.append("%s\t%d\t%d\t%s\t%d\t%s\t" % (rname, a.pos, a.aend, a.qname, a.opt('AS'), '+'))
			
		add_info = "%1.5f\t%d" % (self.score, self.gap);
		l.append(add_info)
		
		return "".join(l);
		
		
	def __cmp__(self, other):
		return  cmp(self.score, other.score)
		
		
		
		
def arlist2chimera(arlist, samfile, gap = 1, overlap = 4, score_function = _alignment_score):
	'''Compiles all possible chimeras from hits(pysam aligned reads) provided
	
		arlist iterable: element is pysam.AlignedRead. All hits of the initial read to compile into chimeras
		gap int: maximum gap allowed between two hits
		gap int: maximum gap allowed between two hits
		
		Returns list: all possible chimeras compiled from hits(pysam aligned reads)
	'''	
	chimeras = []

	for a1, a2 in combinations(arlist, 2):

		if(-overlap <= a2.qstart - a1.qend <= gap):
			chimeras.append(Chimera([a1, a2], samfile, score_function))
			
		elif(-overlap <= a1.qstart - a2.qend <= gap):
			chimeras.append(Chimera([a2, a1], samfile, score_function))
			
	return chimeras;		