# /usr/bin/python
'''Library supporting clustering of mapping hits'''
import sys
import os

from nrlbio.pybedtools_extension import construct_gff_interval

class Cluster(object):
	'''Represents cluster(peak) of mapping hits
	
	Attributes:
		name str: name of the cluster
		chrom str: chromosome cluster lies on 
		start int: start position of the cluster, 0-based, inclusive
		end int: end position of the cluster, 0-based, exclusive
		segments list of pysam.AlignedSegment: mapping hits supporting a cluster
		support int: number of mapping hits supporting the cluster;
	'''
	def __init__(self, name, chrom, start, end, strand="."):
		self.name = name
		self.chrom = chrom
		self.start = start
		self.end = end
		self.strand = strand
		self.support = 0;
		
	def add_segments(self, segments):
		self.segments = segments;
		self.support = len(segments);
		
	def __len__(self):
		return self.end-self.start
	
	def gff(self):
		return construct_gff_interval(self.chrom, self.start, self.end, 'mc', score=str(self.support), strand=self.strand, source='mc', frame='.', attrs=[('ID', self.name)])

