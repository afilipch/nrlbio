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
		self.read_names = read_names;
		# we have to bin aligned_reads according to which interval do they fall
		self.aligned_reads = [[] for i in range(len(intervals))];
		for ar in aligned_reads:
			for c, interval in enumerate(intervals):
				if((interval.chrom == ar.chrom) and (interval.strand == ar.strand) and (ar.start >= interval.start) and (ar.stop <= interval.stop)):
					self.aligned_reads[c].append(ar);
			else:
				sys.stderr.write("%s\nThere is no interval for aligned_read: %s %d %d %s\n\n" % (str("_"*100), ar.chrom, ar.start, ar.stop, ar.qname))
				for interval in intervals:
					sys.stderr.write("interval coordinates: %s %d %d\n" % (interval.chrom, interval.start, interval.stop))
					 
		
	@classmethod
	def from_intervals(cls, name, interacting_intervals):
		r = [];
		for c, intervals in enumerate(interacting_intervals):
			chrom = intervals[0].chrom
			start = min([int(x.start) for x in intervals])
			stop = max([int(x.stop) for x in intervals])
			n = "|".join((name, str(c)))
			score = str(max([float(x[4]) for x in intervals]))
			strand = intervals[0].strand
			gaps = [int(x[10]) for x in intervals]
			gap = min(gaps, key = lambda x: x**2)
			unique_reads = len(intervals)
			interval = construct_gff_interval(chrom, start, stop, 'chimera', score=score, strand=strand, source='chiflex', frame='.', attrs=[("ID", n), ('gap', gap), ('n_uniq', unique_reads)])
			#interval = pybedtools.Interval(chrom, start, stop, n, score, strand, otherfields = [str(gap), str(unique_reads)]) 
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
		#for arw in self.aligned_reads:
			#for a in arw:
				#sys.stderr.write("%s\n" % str(a.aligned_read.query))
			#sys.stderr.write("%s\n" %  str("_"*100))
		
		#set constants
		if (not self.aligned_reads):
			raise AttributeError("HTML attributes cannot be set without aligned_reads. Please add them to Interaction");
		if(not getattr(self, 'extended_intervals', None)):
			raise AttributeError("Interaction.set_extended_intervals() has to be called in advance");
			
		self.icolors = ['green', 'lightblue', 'purple', 'yellow'];
		self.ilinks = [get_link(interval, system, internal = False) for interval in self.intervals]
		##################################################################################################
		
		#set alignments	
		def _set_part(arw, interval):
			if(interval.strand == '+'):
				ladj = arw.start - interval.start; 
			elif(interval.strand == '-'):
				ladj = interval.stop - arw.stop
				
			match_mismatch = [[" "*ladj, ''], ];
			for nref, nread in get_alignment(arw.aligned_read):
				if(not nref):
					continue;
					
				if(nref == nread):
					if(not match_mismatch[-1][1]):
						match_mismatch[-1][0]+=nread;
					else:
						match_mismatch.append([nread, ''])
				else:
					if(nread):
						match_mismatch[-1][1] += nread;
					else:
						match_mismatch[-1][1] += '-';
			#sys.stderr.write("%d\t%s\n" % (ladj, str(match_mismatch)))
			
			#sys.stderr.write("%s\n" % str("*"*100))
			return match_mismatch;			
		
		self.match_mismatch = [[] for i in range(len(self.intervals))]
		for c, (interval, aligned_reads) in enumerate(zip(self.extended_intervals, self.aligned_reads)):
			for arw in aligned_reads:
				self.match_mismatch[c].append(_set_part(arw, interval))
		##################################################################################################
		
		#set read sequences
		def _set_seq(aligned_reads):
			seq_fragments = [];
			seq = aligned_reads[0].aligned_read.seq
			
			borders = [(x[0].aligned_read.qstart, x[0].aligned_read.qend, x[1]) for x in zip(aligned_reads, self.icolors)];
			borders.sort(key = lambda x: x[0]);
			borders.insert(0, (0,0));
			borders.append((len(seq), len(seq)));
			
			for i in range(1, len(borders)-1):
				overlap = seq[borders[i][0]: borders[i-1][1]]
				unmapped = seq[borders[i-1][1]: borders[i][0]]
				mapped = seq[max(borders[i-1][1], borders[i][0]): min(borders[i+1][0], borders[i][1])]
				color = borders[i][2]
				seq_fragments.append((overlap, unmapped, mapped, color));
			else:	
				seq_fragments.append(('', seq[borders[-2][1]: borders[-1][0]], '', 'white'));	
				
			return seq_fragments
						
		self.seq_fragments = [];
		#try:
		for nrow in range(len(self.aligned_reads[0])):
			self.seq_fragments.append(_set_seq([x[nrow] for x in self.aligned_reads]))
		#except:
			#sys.stderr.write("_"*110 + "\n")
			#for ars in self.aligned_reads:
				#for ar in ars:
					#sys.stderr.write("aligned_read: %s %d %d %s\n\n" % (ar.chrom, ar.start, ar.stop, ar.qname))
			#for interval in self.intervals:
				#sys.stderr.write("interval coordinates: %s %d %d\n" % (interval.chrom, interval.start, interval.stop))
				#sys.exit()
			
		##################################################################################################
		
		#set reads' ids
		self.read_ids = [x.qname for x in self.aligned_reads[0]]
			
		
	def doublebed(self):
		'''converts interaction to a doublebed file entry'''
		l = [];
		for i in self.intervals:
			l.append(str(i));
		return "".join(l);	
			

	def __str__(self):
		'''converts interaction to a doublebed file entry'''
		l = [];
		for i in self.intervals:
			l.append(str(i));
		return "".join(l);		
	


		
		
		


		


#testing section
if(__name__ == "__main__"):

	bed = pybedtools.BedTool(sys.argv[1])
	bed2interactions(bed, [12, 12]);