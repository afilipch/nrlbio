'''Library contains functionality to operate with RNA:RNA interactions'''

import sys
#from collections import namedtuple

import pybedtools

from nrlbio.itertools_extension import cmp_attributes
from nrlbio.pybedtools_extension import  generate_overlaping_intervals




class Interaction(object):
	def __init__(self, name, intervals, read_names, attrs = {}):
		self.name = name
		self.intervals = intervals
		self.read_names = read_names
		self.attrs = attrs
	
	@classmethod
	def from_chimeras(cls, name, chimeras):
		chrom1 = chimeras[0][6]
		start1 = min([int(x[7]) for x in chimeras])
		stop1 = min([int(x[8]) for x in chimeras])
		name1 = name #"_".join((name, '1'))
		score1 = str(max([float(x[10]) for x in chimeras]))
		strand1 = chimeras[0][11]
		interval1 = pybedtools.Interval(chrom1, start1, stop1, name1, score1, strand1) 
		
		
		chrom2 = chimeras[0][0]
		start2 = min([int(x[1]) for x in chimeras])
		stop2 = min([int(x[2]) for x in chimeras])
		name2 = name #"_".join((name, '2'))
		score2 = str(max([float(x[4]) for x in chimeras]))
		strand2 = chimeras[0][5]	
		interval2 = pybedtools.Interval(chrom2, start2, stop2, name2, score2, strand2) 
		
		read_names = [x[3] for x in chimeras];
		
		return cls(name, (interval1, interval2), read_names)
		
		
	@classmethod
	def from_intervals(cls, name, interacting_intervals):
		r = [];
		for c, intervals in enumerate(interacting_intervals):
			chrom = intervals[0].chrom
			start = min([int(x.start) for x in intervals])
			stop = min([int(x.stop) for x in intervals])
			n = "|".join((name, str(c)))
			score = str(max([float(x[4]) for x in intervals]))
			strand = intervals[0].strand
			gaps = [int(x[10]) for x in intervals]
			gap = min(gaps, key = lambda x: x**2)
			interval = pybedtools.Interval(chrom, start, stop, n, score, strand, otherfields = [str(gap)]) 
			r.append(interval)

		read_names = [x[3].split("|")[0] for x in interacting_intervals[0]];		
		return cls(name, r, read_names);
		
		
	def doublebed(self):
		'''converts interaction to a doublebed file entry'''
		l = [];
		for i in self.intervals:
			l.append("%s\t%s\t%s\t%s\t%s\t%s\t%s" % tuple(list(i)));
		return "\n".join(l);	
			
		
		
		
	def __str__(self):
		l = [];
		for i in self.intervals:
			l.append("%s\t%s\t%s\t%s\t%s\t%s" % tuple(i));
		return "\t".join(l);	
	
	
	
def bed2interactions(bed, distance, name='interaction'):
	c = 0;
	for m1 in generate_overlaping_intervals(bed, distance[0]):
		tbed = [pybedtools.Interval(i[6], int(i[7]), int(i[8]), i[9], i[10], i[11], otherfields = i[:6] + i[12:]) for i in m1];
		
		for m2 in generate_overlaping_intervals(sorted(tbed, cmp=lambda x,y: cmp_attributes(x,y,['chrom', 'strand', 'start'])), distance[1]):
			c+=1
			yield Interaction.from_chimeras("%s_%d" % (name, c), m2)

		
		
		


		


#testing section
if(__name__ == "__main__"):

	bed = pybedtools.BedTool(sys.argv[1])
	bed2interactions(bed, [12, 12]);