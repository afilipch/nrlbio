'''library contains classes and functions to deal with chimeric reads'''

import sys
from copy import copy
from itertools import combinations;

from nrlbio import numerictools


class ChimeraException(Exception):
    pass

    
#__________________________________________________________________________________
#demultiplexing API

#__________________________________________________________________________________    
    
def as_score(chimera):
	'''score function for Chimera based only on alignment score'''
	return sum(chimera.AS)
	
def as_gap_score(chimera, maxgap=8):
	'''score function for Chimera based on alignment score and gap'''
	if(chimera.gap>0):
		gap = chimera.gap*2
	else:
		gap = chimera.gap
	
	return sum(chimera.AS)*(1 - chimera.gap**2/maxgap**2)	
	
	

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
		self.control = [x.split("_")[0] == "random" for x in self.rnames]
		
		self.AS = [x.opt('AS') for x in aligned_reads]
		self.gap = aligned_reads[1].qstart - aligned_reads[0].qend
		self.score = score_function(self);
		
		
	def __str__(self):
		'''converts Chimera in a bed-like entry(line). First six elements correspond to the first hit, Second six elements to the second one'''
		l = [];
		for a, rname in zip(self.aligned_reads, self.rnames):
			l.append("%s\t%d\t%d\t%s\t%d\t%s\t" % (rname, a.pos, a.aend, a.qname, a.opt('AS'), '+'))
			
		gen_info = "%1.5f\t%d" % (self.score, self.gap);
		l.append(gen_info)
	
		for a in self.aligned_reads:
			l.append("\t%d\t%d\t%s" % (a.qstart, a.qend, a.query[-1]))
		
		return "".join(l);
		
		
	def __cmp__(self, other):
		return  cmp(self.score, other.score)
		
		
		
		
def arlist2chimera(arlist, samfile, gap = 1, overlap = 4, score_function = as_score):
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



def demultiplex(chimeras):
	'''Demultiplex hits derived from the same read (choose the best ones on basis of key_function). Assignes if the read comes from decoy or nonunique.
	
		arwlist list: ArWrappers of the aligned reads(hits) derived from the same reads
		key_function function: method selects only the best hits on basis of key_function 
		
		Returns tuple: 3-element tuple
			1st ArWrapper|None: the best uniquely aligned read, if the best alignment is nonunique originated from decoy
			2nd list: list of nonunique alignments. List is empty if the best uniquely aligned read present
			3rd ArWrapper|None: the best aligned read originated from decoy, None if it is not the best among all alignments for the read
	'''
	real = filter(lambda x: not any(x.control), chimeras);
	control = filter(lambda x: any(x.control), chimeras);
	
	best_real, max_real = numerictools.maxes(real, lambda x: x.score)
	best_control, max_control = numerictools.maxes(control, lambda x: x.score)
	
	if(max_real > max_control):
		if(len(best_real) == 1):
			return best_real[0], [], None
		else:
			return None, best_real, None
	elif(max_control > max_real):
		return None, best_real, best_control[0];
	else:
		return None, best_real, None
		
		
		
		
#__________________________________________________________________________________
#filtering API
#__________________________________________________________________________________

def get_attributes(l, indices):
	'''Converts line in chimera file in a list ready to be passed to filtering'''
	a = l.strip().split("\t");
	for i in indices:
		a[i] = float(a[i]);
	return a;
	

	
def filter_generator(path, indices):
	'''Yields list of attributes corresponding to the chimera. Each list will be used as entry in further filtering
		
		path str:path to the file with chimeras
		indices list of integers: list of indices important for filtering
	'''
	with open(path) as f:
		for l in f:
			yield get_attributes(l, indices);	
			
			
def apply_filter(path, indices, filter_):
	'''Applies given filter to each entry(aligned) in chimeras
		
		path str:path to the file with chimeras
		indices list of integers: list of indices important for filtering
		filter_ str: rule to filter list corresponding to each chimera
		
	Yields str: ready to print representation of chimera
	'''	
	with open(path) as f:
		for l in f:
			x = get_attributes(l, indices)
			if(eval(filter_)):
				yield l
			else:
				pass;


		