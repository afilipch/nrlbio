# /usr/bin/python
'''collections of classes and functions to get statistics of a bed/gff file'''

import os
import sys
import re
from collections import defaultdict;


class Stat(object):
	'''Stat class is a container for a various statistics of bed/gff file. Each type of statistics/feature(length, score, etc.) is represented as collections.defaultdict (Key: certain value of feature e.g. length=50, Value: number of entries which have feature value equal to the Key)
	
	Attributes:
		name string: name of Stat object, may correspond to a bed/gff file name
        length collections.defaultdict: Key: length of the bed interval, Value: number of entries 
        score collections.defaultdict: Key: score of the bed interval, Value: number of entries
	'''
	def __init__(self, name = None):
		self.name = name;
		self.query_start = defaultdict(float)
		self.query_end = defaultdict(float)  

		

	def increment(self, interval):
		self.interval[interval.length] += 1;
		try:
			self.score[interval.score] += 1;
		except:
			self.score['unassigned'] += 1;

			
	def fill_stat_sheet(self, bed_iter, sparse_coefficient = 1):
		'''Extracts statistics of provided iterable containg intervals and stores it in the attributes of the class

		ar_iter iterable: any iterable of pysam.AlignedRead. In the most case an output of pysam.Samfile.fetch()
		sparse_coefficient int: analyses only each {sparse_coefficient}th pysam.AlignedRead;

		Return True if no exception raised
		'''			
		for c, interval in enumerate(bed_iter):
			if (c%sparse_coefficient == 0):
				self.increment(interval)
		
