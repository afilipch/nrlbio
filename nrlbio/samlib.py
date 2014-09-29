# /usr/bin/python
'''collections of classes and functions to deal with sam/bam files'''

import sys;
from collections import namedtuple, Counter;

from nrlbio import numerictools;
from sam_statistics import get_conversions, get_alignment


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
		conversions list of tuples: 1st element in each tuple is nucleotide/gap (string/None) in reference. 2st element in each tuple is nucleotide/gap (string/None) in query.
	'''	
		
	def __init__(self, aligned_read, rname, add_nr_tag = False):
		self.aligned_read = aligned_read;
		self.qname = aligned_read.qname;
		self.rname = rname
		
		self.conversions = get_conversions(self.aligned_read)
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
		tc = Counter(self.conversions)[('T', 'C')]
		self.aligned_read.tags = self.aligned_read.tags + [("TC", tc)];
		
		
		
		
class BackwardWrapper(ArWrapper):		
	def __init__(self, aligned_read, rname, add_nr_tag = False):
		super(BackwardWrapper, self).__init__(aligned_read, rname, add_nr_tag = False);
		l = aligned_read.qname.split("_")
		self.qname = "_".join(l[:-1]);

		if(add_nr_tag):
			number_of_reads = int(self.qname.split("_")[-1][1:])
			self.aligned_read.tags = self.aligned_read.tags + [("NR", number_of_reads)];	
			
		self.tc_pos = int(l[-1].split(":")[-1]);
		if(self.tc_pos >= 0):
			
			if(self.aligned_read.qstart<=self.tc_pos<self.aligned_read.qend):
				pos = self.tc_pos - self.aligned_read.qstart
				self.conversions = [];
				alignment = get_alignment(self.aligned_read);
				p = 0;
				
				for rn, qn in alignment:
					if(p==pos and rn!='C'):
						self.conversions.append((rn, "C"))
					elif(rn!=qn):
						self.conversions.append((rn, qn))
					if(qn):
						p+=1
						
				self.aligned_read.tags = self.aligned_read.tags + [("NT", 1)];		
										
			else:
				self.qname = ''
			
			self.aligned_read.seq = "".join((self.aligned_read.seq[:self.tc_pos], "C", self.aligned_read.seq[self.tc_pos+1:]))
		else:
			pass;
			
		self.aligned_read.qname = self.qname	
			
		#print self.conversions;	
		#print self.aligned_read.seq
		#print
			
		
	
		#for kv in vars(self.aligned_read).items:
			#print "%s\t%s" % kv
		#print 
		#print


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
	'''Converts aligned_read into list corresponding to the attributes provided. We need to do so, since some of attribute of the aligned_read are not accessible via getattr'''
	l = []
	for attr in attributes:
		if(hasattr(ar, attr)):
			l.append(getattr(ar, attr));
		else:
			l.append(ar.opt(attr));
	return l;
	
	
def get_attributes_masked(ar, attributes):
	'''Converts aligned_read into list corresponding to the attributes provided. We need to do so, since some of attribute of the aligned_read are not accessible via getattr'''
	l = []
	for attr in attributes:
		if(hasattr(ar, attr)):
			if(attr == 'qstart'):
				l.append(ar.qstart -  ar.seq.rfind('N'))
			else:	
				l.append(getattr(ar, attr));
		else:
			l.append(ar.opt(attr));
	return l;
	

	
def filter_generator(samfile, attributes, ga = get_attributes):
	'''Yields list of attributes corresponding to the aligned_reads in samfile provided. Each list will be used as entry in further filtering
		
		samfile pysam.Samfile: samfile to generate lists for further filtering
		attributes list: list of attributes important for filtering
	'''
	for aligned_read in samfile.fetch(until_eof=True):
		if(not aligned_read.is_unmapped):
			yield ga(aligned_read, attributes);
	samfile.close()		
			
			
def apply_filter(samfile, attributes, filter_, ga = get_attributes):
	'''Applies given filter to each entry(aligned) in the samfile
		
		samfile pysam.Samfile: samfile to generate lists for further filtering
		attributes list: list of attributes important for filtering. Important: it must be the same as attributes argument in filter_generator
		filter_ str: rule to filter list corresponding to each aligned_read
		
	Yields pysam.AlignedRead: sam entry passed the filtering	
	'''	
	for aligned_read in samfile.fetch(until_eof=True):
		if(not aligned_read.is_unmapped):
			x = ga(aligned_read, attributes)
			if(eval(filter_)):
				yield aligned_read
			else:
				pass;
	samfile.close();			
			
			
			
	
