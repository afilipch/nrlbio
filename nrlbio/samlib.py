# /usr/bin/python
'''collections of classes and functions to deal with sam/bam files'''

import sys;
from collections import namedtuple, Counter;

from nrlbio import numerictools;
from sam_statistics import get_conversions


def key_alignment_score(arw):
	return arw.AS;
	


class ArWrapper(object):
	'''Wrapper for pysam.aligned reads. Adds some additional fields
		
	Attributes:
		aligned_read pysam.aligned_read: read to wrap
		qname str: read id
		rname str: reference id
		AS float: alignment score
		control bool: if True, read comes from decoy
	'''	
		
	def __init__(self, aligned_read, rname, add_nr_tag = False):
		self.aligned_read = aligned_read;
		self.qname = aligned_read.qname;
		self.rname = rname
		
		self.AS = aligned_read.opt("AS")
		if(rname.split("_")[0] == "random"):
			self.control = True;
		else:
			self.control = False;
			
		if(add_nr_tag):
			number_of_reads = int(self.qname.split("_")[-1][1:])
			self.aligned_read.tags = self.aligned_read.tags + [("NR", number_of_reads)];
			
			
	def set_tc(self):
		'''adds number of T->C conversions as a tag 'TC' to the self.aligned_read'''
		conversions = get_conversions(self.aligned_read)
		tc = Counter(conversions)[('T', 'C')]
		self.aligned_read.tags = self.aligned_read.tags + [("TC", tc)];
		
		
		
		
class BackwardWrapper(ArWrapper):		
	def __init__(self, aligned_read, rname):
		super(BackwardWrapper, self).__init__(aligned_read, rname, add_nr_tag = False);


def demultiplex_read_hits(arwlist, key_function):
	'''Demultiplex hits derived from the same read (choose the best ones on basis of key_function). Assignes if the read comes from decoy or nonunique.
	
		arwlist list: ArWrappers of the aligned reads(hits) derived from the same reads
		key_function function: method selects only the best hits on basis of key_function 
		
		Returns tuple: 3-element tuple
			1st ArWrapper|None: the best uniquely aligned read, if the best alignment is nonunique originated from decoy
			2nd list: list of nonunique alignments. List is empty if the best uniquely aligned read present
			3rd ArWrapper|None: the best aligned read originated from decoy, None if it is not the best among all alignments for the read
	'''
			
	real = filter(lambda x: not x.control, arwlist);
	control = filter(lambda x: x.control, arwlist);
	best_real, max_real = numerictools.maxes(real, key_function)
	best_control, max_control = numerictools.maxes(control, key_function)
	if(max_real > max_control):
		if(len(best_real) == 1):
			return best_real[0], [], None
		else:
			return None, best_real, None
	elif(max_control > max_real):
		return None, best_real, best_control[0];
	else:
		return None, best_real, None

		
		
		
def get_attributes(ar, attributes):
	l = []
	for attr in attributes:
		if(hasattr(ar, attr)):
			l.append(getattr(ar, attr));
		else:
			l.append(ar.opt(attr));
	return l;		

	
def filter_generator(samfile, attributes):
	for aligned_read in samfile.fetch(until_eof=True):
		if(not aligned_read.is_unmapped):
			yield get_attributes(aligned_read, attributes);
	samfile.close()		
			
def apply_filter(samfile, attributes, filter_):
	for aligned_read in samfile.fetch(until_eof=True):
		if(not aligned_read.is_unmapped):
			x = get_attributes(aligned_read, attributes)
			if(eval(filter_)):
				yield aligned_read
			else:
				pass;
	samfile.close();			
			
			
			
	
