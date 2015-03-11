'''Library contains functionality to operate with RNA:RNA interactions'''

import sys
from copy import copy
#from collections import namedtuple

import pybedtools

from nrlbio.itertools_extension import cmp_attributes
from nrlbio.pybedtools_extension import  interval2seq, construct_gff_interval
from nrlbio.statistics.sam import get_alignment;
from nrlbio.html import get_link


	

class Interaction(object):
	def __init__(self, name, intervals, aligned_reads = [], read_names = []):
		self.name = name;
		self.intervals = intervals;
		self.aligned_reads = aligned_reads;
		self.read_names = read_names;
		
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
		return cls(name, r, read_names = read_names);
		
		
	def set_extended_intervals(self, reference=None, extensions=None):
		'''Sets extended intervals and sequences for the interaction
		
			interaction Interaction: interaction between intervals
			reference dict: key - chrmosome(feature name), value - Bio.SeqRecord.SeqRecord object. Normally reference odject is created via 'Bio.SeqIO.to_dict(Bio.SeqIO.parse([path_to_fasta_file], "fasta"))' call
			extensions iterable of iterables:Controls how far interacting intervals' sequences will be extended. For example extensions=[(1,4), (0,7)] will extend first interval 1nt upstream and 4nts downstream, while the second interaction will be extended 0nt upstream and 7nts downstream
		'''
		if(extensions):
			if(not reference):
				raise ValueError('\'reference\' argument cannot be None, while \'extensions\' is set to be not None\n')
			
			self.extended_intervals = [];
			for interval, (eu, ed) in zip(self.intervals, extensions):
				if(interval.strand == "+"):
					estart = max(0, interval.start - eu);
					estop = interval.stop+ed;
				else:
					estart = max(0, interval.start - ed);
					estop = interval.stop+eu;
				newint = construct_gff_interval(interval.chrom, estart, estop, "extinterval", score='0', strand=interval.strand, source='.', frame='.', attrs=[])
				newint.attrs['seq'] = interval2seq(newint, reference)
				self.extended_intervals.append(newint)
				
		elif(reference):
			for interval in self.intervals:
				interval.attrs['seq'] = interval2seq(interval, reference)
			self.extended_intervals = self.intervals;
			
		else:
			self.extended_intervals = self.intervals;
		#for interval in self.extended_intervals:
			#print interval.file_type
		#print "_"*120	
			
		
		
	def	set_html_attributes(self, system):
		if (not self.aligned_reads):
			raise AttributeError("HTML attributes cannot be set without aligned_reads. Please add them to Interaction");
		if(not getattr(self, 'extended_intervals', None)):
			raise AttributeError("Interaction.set_extended_intervals() has to be called in advance");
			
		self.icolors = ['green', 'lightblue', 'purple', 'yellow'];
		self.ilinks = [get_link(interval, system, internal = False) for interval in self.extended_intervals]
			
		def _set_part(arw):
			pass;
			
		
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
	


		
		
		


		


#testing section
if(__name__ == "__main__"):

	bed = pybedtools.BedTool(sys.argv[1])
	bed2interactions(bed, [12, 12]);