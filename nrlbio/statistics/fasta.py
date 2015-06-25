# /usr/bin/python
'''collections of classes and functions to get statistics of a sequence file(fasta/fastq/genbank)'''

import os
import sys
from collections import defaultdict, Counter;

from nrlbio.config.config import load_config
from nrlbio.numerictools import dict2simple_stat

class Stat(object):
	'''Stat class is a container for a various statistics of fasta file. Each type of statistics/feature(length, nucleotide composition, etc.) is represented as collections.defaultdict (Key: certain value of feature e.g. length=50, Value: number of entries which have feature value equal to the Key)
	
	Attributes:
		name string: name of Stat object, may correspond to a fasta file name
        length collections.defaultdict: Key: length of the fasta entry, Value: number of entries 
        nucleotides collections.defaultdict: Key: nucleotide, Value: number of the nucleotide in fasta entries
	'''
	
	_hist_names = ["length"]
	
	def __init__(self, name = None):
		if(name):
			self.name = name;
		else:
			self.name = 'unknown'
		self.length = defaultdict(int)
		self.nucleotides = Counter(); 
		

	def increment(self, seq):
		self.length[len(seq)] += 1;
		self.nucleotides += Counter(seq)


			
	def fill_stat_sheet(self, seq_iter, sparse_coefficient=1):
		'''Get statistics of provided iterable of sequences and stores it in the attributes of the object

		ar_iter iterable: any iterable of pysam.AlignedRead. In the most case an output of pysam.Samfile.fetch()
		sparse_coefficient int: statistics is generated only for each {sparse_coefficient}th entry;

		'''
		if(sparse_coefficient==1):
			for seq in seq_iter:
				self.increment(seq);
				
		else:
			for c, interval in enumerate(bed_iter):
				if (c%sparse_coefficient == 0):
					self.increment(interval, attributes=attributes)
					
					

	@classmethod				
	def from_file(cls, path, ftype, sparse_coefficient=1):
		'''Creates Stat instance and fills it using sequences from provided file(fasta/fastq/genbank). NOTE: Stat name is set to os.path.basename(path)
			
			path str: path to fasta/fastq/genbank file to get statistics on
			ftype str: format of files with seqrecords('genbank', 'fastq', 'fasta', ets.)
			sparse_coefficient int: statistics is generated only for each {sparse_coefficient}th entry;
		'''
		from nrlbio.generators import generator_seq
		obj = cls(os.path.basename(path))
		obj.fill_stat_sheet(generator_seq([path], ftype), sparse_coefficient)
		return obj;	 
				
				
	def generate_plots(self, configuration='fastastat', output=''):
		pass;
		
		
	
	def __str__(self):
		nentries, length_total, length_mean, length_median, length_min, length_max = dict2simple_stat(self.length);
		
		sl = "Statistics on sequence file: %s\n\nNumber of sequences: %d\nTotal length: %d\nMean length: %1.1f\nMedian length: %1.1f\nMin length: %d\nMax length: %d\n" % (self.name, nentries, length_total, length_mean, length_median, length_min, length_max);
		
		sn = "\n".join(["%s: %d" % (x, self.nucleotides.get(x,0)) for x in "ACTGN"])
		
		return "%s%s\n\nNucleotide Composition:\n%s" % ("_"*140, sl, sn)
		