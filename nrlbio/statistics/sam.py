# /usr/bin/python
'''collections of classes and functions to produce various statistics of mapping'''

import os
import sys
import re
from collections import defaultdict;

import Bio;
import pysam;
#import jinja2;

from nrlbio.config.config import load_config
from nrlbio import pyplot_extension
#from nrlbio import html;



'''pattern to get matched, mismatched and del/ins stretches to a reference from MD field in sam/bam file file'''
MD_pattern = '([0-9]+)([A-Z]+)*(\^[A-Z]+)*';
MD_compiled_pattern = re.compile(MD_pattern);

def intermediate_alignment(ar):
	'''function to get all mismatches/deletions at defined positions of the aligned read based on the aligned sequence and MD field in sam/bam file. The output is used further to define mismatches/deletions/insertions or/and alignment

	ar pysam.AlignedRead: aligned read to be analyzed

	Return tuple. 1st element is mismatch dictionary (Key: position of mismatch in query, Value: corresponding nucleotide in reference). 2nd element is listiterator (Element: Deleted nucleotides in an order they appear in reference)
	'''	
	mdlist = MD_compiled_pattern.findall(ar.opt('MD'));
	mmdict = {};
	dellist = [];
	p = 0;
	for l in mdlist:
		p += int(l[0]);
		for n in l[1]:
			mmdict[p] = n;
			p+=1;			
		for n in l[2][1:]:
			dellist.append(n);		
		
	deliter = iter(dellist);
	return mmdict, deliter

def get_alignment(ar):	
	'''function to get actual alignment for given pysam.AlignedRead
	
	ar pysam.AlignedRead: aligned read to be analyzed

	Return list of tuples. 1st element in each tuple is nucleotide/gap (string/None) in reference. 2st element in each tuple is nucleotide/gap (string/None) in query.
	'''	
	alignment = [];
	mmdict, deliter = intermediate_alignment(ar);
	ins_adjust = 0;	
	for i, j in ar.aligned_pairs:
		if(i != None):
			if(j != None):
				if(i - ins_adjust in mmdict):
					apair = (mmdict[i-ins_adjust], ar.query[i]); 
				else:
					apair = (ar.query[i], ar.query[i]);
			else:
				ins_adjust += 1;
				apair = (None, ar.query[i])
		else:
			apair = (next(deliter), None)
		alignment.append(apair);	
			
	return alignment;
	
def get_conversions(ar):
	'''function to get all mismatches/deletions/insertions for given pysam.AlignedRead
	
	ar pysam.AlignedRead: aligned read to be analyzed

	Return list of tuples. 1st element in each tuple is nucleotide/gap (string/None) in reference. 2st element in each tuple is nucleotide/gap (string/None) in query.
	'''	
	conversions = [];
	mmdict, deliter = intermediate_alignment(ar);
	ins_adjust = 0;	
	for i, j in ar.aligned_pairs:
		if(i != None):
			if(j != None):
				if(i - ins_adjust in mmdict):
					apair = (mmdict[i-ins_adjust], ar.query[i]); 
				else:
					continue;
			else:
				ins_adjust += 1;
				apair = (None, ar.query[i])
		else:
			apair = (next(deliter), None)
		conversions.append(apair);	
			
	return conversions;
	
	
def conv2labels(conversions):
	_order = [(None, None), ('A', 'C'), ('A', 'T'), ('A', 'G'), ('C', 'A'), ('C', 'T'), ('C', 'G'), ('T', 'A'), ('T', 'C'), ('T', 'G'), ('G', 'A'), ('G', 'C'), ('G', 'T'), ('G', 'G'), ('A', None), ('C', None), ('G', None), ('T', None), (None, 'A'), (None, 'C'), (None, 'T'), (None, 'G')]
	
	labels = [];
	ncounter = {};
	for c, mtype in enumerate(_order):
		ncounter[c+1] = conversions[mtype];
		
		if(mtype[0] and mtype[1]):
			labels.append("%s->%s" % mtype);
		elif((not mtype[0]) and (not mtype[1])):
			labels.append("perf")
		elif(not mtype[0]):
			labels.append("ins %s" % mtype[1]);
		elif(not mtype[1]):
			labels.append("del %s" % mtype[0]);
			
	return ncounter, labels		
	
	
class Stat(object):
	'''Stat class is a container for a various statistics of mapping. Each type of statistics/feature(conversions, alignment score etc.) is represented as collections.defaultdict (Key: certain value of feature e.g. alignment_score=50, Value: number of aligned_reads which have feature value equal to the Key)
	
	Attributes:
		name string: name of Stat object, may correspond to a set of pysam.AlignedRead analyzed
        query_start collections.defaultdict: Key: start position of the aligned portion of the read sequence (0-based, inclusive) on read sequence, Value: number of aligned reads corresponding the Key value;
        query_end collections.defaultdict: Key: end position of the aligned portion of the read sequence (0-based, exclusive) on read sequence, Value: number of aligned reads corresponding the Key value;
        conv collections.defaultdict: Key: type of conversion tuple(nt in reference, nt in read sequence), Value: number of conversions corresponding the Key value;
        conv_weighted collections.defaultdict: Key: type of conversion tuple(nt in reference, nt in read sequence), Value: number of aligned reads corresponding the Key value;
        conv_number collections.defaultdict: Key: number of conversions, Value: number of aligned reads corresponding the Key value;       
        as_score collections.defaultdict: Key: alignment score, Value: number of aligned reads corresponding the Key value;
        clipped_length_right collections.defaultdict: Key: number of nucleotides soft clipped from the right, Value: number of aligned reads corresponding the Key value;
        clipped_seq_left collections.defaultdict: Key: sequence soft clipped from the left, Value: number of aligned reads corresponding the Key value;
        clipped_seq_right collections.defaultdict: Key: sequence soft clipped from the right, Value: number of aligned reads corresponding the Key value;   
        ref_start collections.defaultdict: Key: start position of the aligned portion of the read sequence (0-based, inclusive) on reference, Value: number of aligned reads corresponding the Key value;
        ref_end collections.defaultdict: Key: end position of the aligned portion of the read sequence (0-based, exclusive) on reference, Value: number of aligned reads corresponding the Key value;
        query_ref_start collections.defaultdict: Key: tuple (start position of the aligned portion of the read sequence (0-based, exclusive) on read sequence, start position of the aligned portion of the read sequence (0-based, inclusive) on reference), Value: number of aligned reads corresponding the Key value;
        _counters_names list of str: names(names of object attributes) of the counters(statistics attributes).
        _counters_titles list of str: titles of the counters(statistics attributes) which can be used as table/plot titles
	'''
	def __init__(self, name = None):
		self.name = name;
		self.query_start = defaultdict(float)#is equal to clipped_length_left, which is therefore redundant and therefore is not used
		self.query_end = defaultdict(float)  
		self.conv = defaultdict(float)
		self.conv_weighted = defaultdict(float) 
		self.conv_number = defaultdict(float)
		self.ascore = defaultdict(float)
		self.clipped_length_right = defaultdict(float)	
		
		#these attributes are used in detailed analysis, could triger memory overload in case of frequent soft-clipping
		self.clipped_seq_left = defaultdict(int)
		self.clipped_seq_right = defaultdict(int) 
		
		#these attributes are used for only the analysis of mapping to the short reference (miRNA, piRNA etc.)
		self.ref_start = defaultdict(float)  
		self.ref_end = defaultdict(float)    
		self.query_ref_start = defaultdict(float) 
		#self.cut_right = defaultdict(float) 
		
		self._counters_names = ["ascore", "query_start", "clipped_length_right", "query_end", "conv", "conv_weighted", "conv_number",	"clipped_seq_left", "clipped_seq_right", "query_ref_start", "ref_start", "ref_end"]
		self._counters_titles = ["Alignment Score", "Start position of the match in query", "number of nucleotides soft clipped downstream", "End position of the match in query", "Type of conversion", "Type of conversion weighted to a number of conversions in one read", "Number of conversion per read", "Soft clipped upstream sequence", "Soft clipped downstream sequence", "Start position of the match in query and reference", "Start position of the match in reference", "End position of the match in reference"]
		
	def increment_basic(self, ar, conversions=None):
		self.query_start[ar.qstart] += 1;
		self.query_end[ar.qend] += 1;
		self.ascore[ar.opt('AS')] += 1;
		self.clipped_length_right[ar.rlen - ar.qend] += 1;
		
		if(not conversions):
			conversions = get_conversions(ar)
		cn = len(conversions);
		self.conv_number[cn]+=1;
		if(cn):
			for c in conversions:
				self.conv[c] += 1;
				self.conv_weighted[c] += 1.0/cn;
		else:	
			self.conv[(None, None)] += 1;
			self.conv_weighted[(None, None)] += 1.0;
		return True;
				
				
	def increment_short(self, ar):
		self.ref_start[ar.pos] += 1;
		self.ref_end[ar.aend] += 1;
		self.query_ref_start[(ar.qstart, ar.pos)] += 1;
		return True;
		

	def increment_detailed(self, ar):
		self.clipped_seq_left[ar.seq[:ar.qstart]] +=1;
		self.clipped_seq_right[ar.seq[ar.qend:]] += 1;
		return True;

	def increment_all(self, ar):
		self.increment_basic(ar)
		self.increment_short(ar)
		self.increment_detailed(ar)	
		return True
		

			
	def fill_stat_sheet(self, ar_iter, short_reference = False, detailed = False, sparse_coefficient = 1):
		'''Extracts statistics of provided iterable containg aligned reads and stores it in the attributes of the class

		ar_iter iterable: any iterable of pysam.AlignedRead. In the most case an output of pysam.Samfile.fetch()
		short_reference bool: if True additional statistics is collected. Makes sense for a mapping to the reference composed of short reads(piRNA, miRNA, ncRNA etc.)
		detailed bool: if True collects advanced statistics, which can require more memory
		sparse_coefficient int: analyses only each {sparse_coefficient}th pysam.AlignedRead;

		Return True if no exception raised
		'''			
		if(short_reference):
			if(detailed):
				for i, ar in enumerate(ar_iter):
					if (i%sparse_coefficient == 0):
						self.increment_all(ar)
			else:		
				for i, ar in enumerate(ar_iter):
					if (i%sparse_coefficient == 0):
						self.increment_basic(ar)
						self.increment_short(ar)
		elif(detailed):
			for i, ar in enumerate(ar_iter):
				if (i%sparse_coefficient == 0):
					self.increment_basic(ar)
					self.increment_detailed(ar)
		else:
			for i, ar in enumerate(ar_iter):
				if (i%sparse_coefficient == 0):
					self.increment_basic(ar)
		return True
		
	def generate_plots(self, configuration, output=''):
		n=5
		config_ = load_config(configuration)
		for attr_name, title in zip(self._counters_names[:n], self._counters_titles[:n]):
			attr = getattr(self, attr_name)
			if(attr):
				kwargs = config_[attr_name];
				kwargs['output'] = os.path.join(output, "%s.pdf" % attr_name);
				kwargs['label'] = self.name;
				if(attr_name in ["conv", "conv_weighted"]):
					ncounter, xticklabels = conv2labels(attr)
					kwargs['xticklabels'] = xticklabels
					kwargs['xticks'] = sorted(ncounter.keys())
					pyplot_extension.histogram(ncounter, **kwargs)
				else:	
					pyplot_extension.histogram(attr, **kwargs)

		

if (__name__=='__main__'):
	samfile = pysam.Samfile(sys.argv[1])
	samstat = Stat('test');	
	samstat.fill_stat_sheet(samfile.fetch(until_eof=True), short_reference = False, detailed = False, sparse_coefficient = 1)
	samstat.generate_plots('samstat')
	#print samstat.conv
	