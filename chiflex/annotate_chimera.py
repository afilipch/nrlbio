#! /usr/lib/python
'''assignes the type(linera/circular splice junctions, intra/inter-molecular interaction) to each chimera/interaction.''' 
import argparse
import sys;
from collections import defaultdict, OrderedDict

from pybedtools import BedTool, Interval

from nrlbio.pybedtools_extension import list2interval, bed2gff

parser = argparse.ArgumentParser(description='assignes the type(linear/circular splice junctions, intra/inter-molecular interaction) to each chimera/interaction');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to bed file");
parser.add_argument('-e', '--exons', nargs = '?', required = True, type = str, help = "Path to the bed file of transcripts\' exons");
parser.add_argument('-d', '--distance', nargs = '?', default = 2, type = int, help = "Maximum distance between region and exons edges allowed");
parser.add_argument('-s', '--stranded', nargs = '?', const=True, default = False, type = bool, help = "Are sequenced reads stranded");
args = parser.parse_args();

#sys.stderr.write(str(args.stranded))



class AnnotatedChimera(object):
	def __init__(self, intervals, intersections):
		self.intervals = intervals
		self.intersections = intersections
		
		
	def annotate_intramolecular_chimera(self, i1, i2, distance):
		distance_start_first = abs(self.intervals[0].start - i1.start)
		distance_stop_first = abs(self.intervals[0].stop - i1.stop)
		distance_start_second = abs(self.intervals[1].start - i2.start)
		distance_stop_second = abs(self.intervals[1].stop - i2.stop)
		
		exon_number_difference = int(i2.score) - int(i1.score)
		
		sense_condition = distance_stop_first<=distance and distance_start_second<=distance
		antisense_condition = distance_start_first<=distance and distance_stop_second<=distance
		
		if(sense_condition and exon_number_difference>0):
			return 'lsj'; 
		elif(antisense_condition and exon_number_difference < 0):
			return 'lsj'
		
		elif(sense_condition and exon_number_difference<=0):
			return 'csj'; 
		elif(antisense_condition and exon_number_difference>=0):
			return 'csj'
		
		else:
			return "intra";
					
				
	def annotate(self, distance):
		types = [];
		for transcript_name, d in self.intersections.items():
			if(len(d) == 2 and all(d.values())):
				chimera_type = self.annotate_intramolecular_chimera(d[0], d[1], distance)
				if(chimera_type=="lsj"):
					return chimera_type;
				else:	
					types.append(chimera_type);
				  
		for t in ['csj', 'intra']:
			if(t in types):
				return t;
				
		return 'inter';

	
	def yield_annotated(self, distance, offset):
		for k, interval in self.intervals.items():
			if(interval.file_type=='bed'):
				gi = bed2gff(interval, feature='ch')
				gi.attrs['qstart'] = interval[6]
				gi.attrs['qend'] = interval[7]
				gi.attrs['chscore'] = interval[9]
				gi.attrs['gap'] = interval[10]
			else:
				gi = interval;
				
			gi.attrs['interaction'] = self.annotate(distance)
			yield gi[:offset];
			
			
	def __str__(self):
		l = []
		for k, i in self.intervals.items():
			l.append("%s:\t%s" % (k, str(i).strip()));
			
		l.append("\n")	
		
		for k1, d in self.intersections.items():
			for k2, i in d.items():
				l.append("%s:\t%s:\t%s" % (k1, k2, str(i).strip()));
			
		return "\n".join(l)	





chimeras = BedTool(args.path)
exons = BedTool(args.exons)
offset  = chimeras.field_count();
chimeras_vs_exons = chimeras.intersect(exons, s=args.stranded, wao=True)


first = chimeras_vs_exons[0]

curname, cur_pair_number = first.name.split("|");
intersections = defaultdict(dict);
intervals = OrderedDict();
intervals[int(cur_pair_number)] = first;
intersections[first[offset+3]][int(cur_pair_number)] = list2interval(first[offset:])

for i in chimeras_vs_exons[1:]:
	name, pair_number = i.name.split("|");
	if(name == curname):
		if(pair_number == cur_pair_number):
			intersections[i[offset+3]][int(pair_number)] = list2interval(i[offset:])
		else:
			cur_pair_number = pair_number;
			intervals[int(pair_number)] = i;
			intersections[i[offset+3]][int(pair_number)] = list2interval(i[offset:])
	else:
		ac = AnnotatedChimera(intervals, intersections);
		for ani in ac.yield_annotated(args.distance, offset):
			print "\t".join(ani)
		curname = name
		cur_pair_number = pair_number;
		intersections = defaultdict(dict);
		intervals = OrderedDict();
		intervals[int(pair_number)] = i;
		intersections[i[offset+3]][int(pair_number)] = list2interval(i[offset:])
		
else:
	ac = AnnotatedChimera(intervals, intersections);
	for ani in ac.yield_annotated(args.distance, offset):
		print "\t".join(ani)



























