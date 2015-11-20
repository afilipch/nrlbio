# /usr/bin/python
'''collections of classes and functions to get statistics of a bed/gff file'''

import os
import sys
from collections import defaultdict;

from nrlbio.config.config import load_config

class Stat(object):
	'''Stat class is a container for a various statistics of bed/gff file. Each type of statistics/feature(length, score, etc.) is represented as collections.defaultdict (Key: certain value of feature e.g. length=50, Value: number of entries which have feature value equal to the Key)
	
	Attributes:
		name string: name of Stat object, may correspond to a bed/gff file name
        length collections.defaultdict: Key: length of the bed interval, Value: number of entries 
        score collections.defaultdict: Key: score of the bed interval, Value: number of entries
	'''
	_hist_names = ["length", "score"]
	
	def __init__(self, name = None):
		self.name = name;
		self.length = defaultdict(int)
		self.score = defaultdict(float)  
		self.attrs = defaultdict(lambda: defaultdict(int))
		

	def increment(self, interval, attributes=[]):
		self.length[interval.length] += 1;
		try:
			self.score[float(interval.score)] += 1;
		except:
			pass;
		for a in attributes:
			try:
				ma = filter(bool, interval.attrs[a].split(","))
				for sa in ma:
					self.attrs[a][sa] += 1.0/len(ma);
			except:
				self.attrs[a]['unassigned'] += 1;

			
	def fill_stat_sheet(self, bed_iter, attributes=[], sparse_coefficient = 1):
		'''Get statistics of provided iterable of intervals and stores it in the attributes of the object

		ar_iter iterable: any iterable of pysam.AlignedRead. In the most case an output of pysam.Samfile.fetch()
		sparse_coefficient int: statistics is generated only for each {sparse_coefficient}th entry;
		'''			
		for c, interval in enumerate(bed_iter):
			if (c%sparse_coefficient == 0):
				self.increment(interval, attributes=attributes)
				

				
				
	def generate_plots(self, configuration='bedstat', output=''):
		from nrlbio import pyplot_extension
		config_ = load_config(configuration)
		for attr_name in self._hist_names:
			kwargs = config_[attr_name];
			kwargs['output'] = os.path.join(output, "%s.png" % attr_name);
			kwargs['label'] = self.name;
			attr = getattr(self, attr_name)
			if(attr):
				#print attr_name;
				pyplot_extension.histogram(attr, **kwargs);
			else:
				pass;
		for attr_name, attr in self.attrs.iteritems():
			kwargs = config_[attr_name];
			kwargs['output'] = os.path.join(output, "%s.png" % attr_name);
			if(attr):
				try:
					ncounter = dict([ (float(x[0]), x[1]) for x in attr.items() ])
					flag = "hist"
				except:
					ncounter = attr;
					flag = "pie"
					
				if(flag == "hist"):
					pyplot_extension.histogram(ncounter, **kwargs);
				elif(flag == "pie"):
					pyplot_extension.pie(ncounter, **kwargs);
					
			else:
				pass
					
			
			
			
			
			
			
			
			
			
			
			
			
			
		
