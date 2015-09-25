'''library contains classes and functions to deal with chimeric reads'''

import sys
from copy import copy
from itertools import combinations;

from Bio.Seq import reverse_complement


from nrlbio import numerictools


class ChimeraException(Exception):
    pass

    
#__________________________________________________________________________________
#demultiplexing API

#__________________________________________________________________________________    
    
strand_conv = {True: '-', False: '+'}
    
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
		
		
	def __init__(self, ar_wrappers, score_function):
		if(len(set([x.qname for x in ar_wrappers])) != 1): 
			raise ChimeraException('Chimera cannot be made from aligned reads with differenr identifiers\nFollowing are given:\n%s\n' % "\n".join(["\t%s" % x for x.qname in ar.wrappers]))
			
		self.ar_wrappers = ar_wrappers;
		self.control = [x.rname.split("_")[0] == "random" for x in self.ar_wrappers]
		
		self.gap = ar_wrappers[1].qstart - ar_wrappers[0].qend;		
		self.AS = [x.AS for x in ar_wrappers];
		self.score = score_function(self);
		
		if(ar_wrappers[0].aligned_read.is_reverse):
			self.gap_seq = reverse_complement(ar_wrappers[0].aligned_read.query_sequence)[self.ar_wrappers[1].qstart:self.ar_wrappers[0].qend]
		else:
			self.gap_seq = ar_wrappers[0].aligned_read.query_sequence[self.ar_wrappers[1].qstart:self.ar_wrappers[0].qend]
			
		
		
	def __str__(self):
		'''converts Chimera in a bed-like entry(line). First six elements correspond to the first hit, Second six elements to the second one'''
		l = [];
		for arw in self.ar_wrappers:
			l.append("%s\t%d\t%d\t%s\t%d\t%s\t" % (arw.rname, arw.aligned_read.pos, arw.aligned_read.aend, arw.qname, arw.AS, strand_conv[arw.aligned_read.is_reverse]))
			
		gen_info = "%1.5f\t%d" % (self.score, self.gap);
		l.append(gen_info)
	
		for a in self.aligned_reads:
			l.append("\t%d\t%d\t%s" % (a.qstart, a.qend, a.aligned_read.query[-1]))
		
		return "".join(l);	
		
		
	def __cmp__(self, other):
		return  cmp(self.score, other.score)
		
		
	def doublebed(self):
		'''Converts chimera in two bed entries with the same identifiers'''
		l = [];
		gen_info = "%1.5f\t%d" % (self.score, self.gap);
		
		if(self.gap_seq):
			gap_seq = self.gap_seq;
		else:
			gap_seq = str(self.gap)
			
		for c, arw in enumerate(self.ar_wrappers):
			l.append("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%s" % (arw.rname, arw.aligned_read.pos, arw.aligned_read.aend, "|".join((arw.qname, str(c))), arw.AS, strand_conv[arw.aligned_read.is_reverse], arw.qstart, arw.qend, gap_seq, gen_info));
			
		return "\n".join(l)
		
		
def arwlist2chimera(arwlist, gap = 100, overlap = 100, score_function = as_score):
	'''Compiles all possible chimeras from hits(pysam aligned reads) provided
	
		arlist iterable: element is pysam.AlignedRead. All hits of the initial read to compile into chimeras
		gap int: maximum gap allowed between two hits
		gap int: maximum gap allowed between two hits
		
		Returns list: all possible chimeras compiled from hits(pysam aligned reads)
	'''	
	chimeras = []
	#print "*"*120
	for a1, a2 in combinations(arwlist, 2):
		#print a1
		#print a2
		#print a1.qstart, a1.qend, a2.qstart, a2.qend
		#print "_"*120
		if((-overlap <= a2.qstart - a1.qend <= gap) or (-overlap <= a1.qstart - a2.qend <= gap)):
			if(a1.qstart < a2.qstart):
				chimeras.append(Chimera([a1, a2], score_function))			
			else:
				chimeras.append(Chimera([a2, a1], score_function))
			
	return chimeras;
	



def demultiplex(chimeras, min_difference):
	'''Demultiplex chimeras derived from the same read (choose the best ones on basis of its score). Assignes if the chimera comes from real or decoy refrence, ornonunique.
	
		Chimeras list: list of Chimeras derived from the same reads
		min_difference int: minimal distance allowed between the best and the second best chimera. If the actual distance is less, than hit will be assigned as nonunique
		
		Returns tuple: 3-element tuple
			1st Chimera|None: the best unique chimera, if the best chimera is nonunique originated from decoy
			2nd list: list of nonunique chimeras. List is empty if the best unique chimera present
			3rd Chimera|None: the best chimera originated from decoy, None if it is not the best among all possible chimeras for the read
	'''
	if(not chimeras):
		return None, [], None
	
	real = filter(lambda x: not any(x.control), chimeras);
	control = filter(lambda x: any(x.control), chimeras);
	
	best_real, max_real = numerictools.maxes(real, lambda x: x.score)
	best_control, max_control = numerictools.maxes(control, lambda x: x.score)
	
	if(max_real):
		distances = [max_real - x.score for x in chimeras if x.score!=max_real]
		if(distances):
			dif_from_best = min(distances);
		else:
			dif_from_best = min_difference + 1;
	
	if(max_real > max_control):
		if(len(best_real) == 1 and min_difference<dif_from_best):
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

def get_attributes(a, indices):
	'''Converts line in chimera file in a list ready to be passed to filtering'''
	for i in indices:
		if(isinstance(a[i], int)):
			pass;
		elif(isinstance(a[i], str) and numerictools.isinteger(a[i])):
			a[i] = int(a[i]);
			#sys.stderr.write("%s\t%s\n" % (i, a[i]))
		else:
			a[i] = float(a[i]);
	return a;
	

	
def filter_generator(path, indices):
	'''Yields list of attributes corresponding to the chimera. Each list will be used as entry in further filtering
		
		path str:path to the file with chimeras
		indices list of integers: list of indices important for filtering
	'''
	with open(path) as f:
		for l in f:
			a = l.strip().split("\t");
			yield get_attributes(a, indices);
			
			
def filter_generator_doublebed(path, indices):
	'''USE FOR DOUBLEBED FORMAT. Yields list of attributes corresponding to the chimera. Each list will be used as entry in further filtering
		
		path str:path to the file with chimeras
		indices list of integers: list of indices important for filtering
	'''
	i = 1;
	with open(path) as f:
		for l in f:
			if(i % 2):
				a1 = l.strip().split("\t");
			else:
				a2 = l.strip().split("\t");
				a = a1[0:6] + a2[0:6] + a1[6:11] + a2[6:]
				yield get_attributes(a, indices);	
			i+=1	
			
			
def apply_filter(path, indices, filter_):
	'''Applies given filter to each entry(aligned) in chimeras
		
		path str:path to the file with chimeras
		indices list of integers: list of indices important for filtering
		filter_ str: rule to filter list corresponding to each chimera
		
	Yields str: ready to print representation of chimera
	'''	
	with open(path) as f:
		for l in f:
			a = l.strip().split("\t");
			x = get_attributes(a, indices)
			if(eval(filter_)):
				yield l;
			else:
				pass;
				

def apply_filter_doublebed(path, indices, filter_):
	'''USE FOR DOUBLEBED FORMAT. Applies given filter to each entry(aligned) in chimeras
		
		path str:path to the file with chimeras
		indices list of integers: list of indices important for filtering
		filter_ str: rule to filter list corresponding to each chimera
		
	Yields str: ready to print representation of chimera
	'''	
	for x in filter_generator_doublebed(path, indices):
			if(eval(filter_)):
				s1 = "\t".join([str(e) for e in x[0:6] + x[12:17]]);
				s2 = "\t".join([str(e) for e in x[6:12] + x[17:]]);
				yield "\n".join((s1, s2));
			else:
				pass;		



				
				
#testing section
if(__name__ == "__main__"):
	pass