#! /usr/lib/python
'''assignes the type(linera/circular splice junctions, intra/inter-molecular interaction) to each chimera/interaction.''' 
import argparse
import sys;
from collections import defaultdict, OrderedDict

from pybedtools import BedTool, Interval;

from nrlbio.pybedtools_extension import list2interval

parser = argparse.ArgumentParser(description='assignes the type(linear/circular splice junctions, intra/inter-molecular interaction) to each chimera/interaction');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to bed file");
parser.add_argument('-e', '--exons', nargs = '?', required = True, type = str, help = "path to the bed file of transcripts\' exons");
parser.add_argument('-d', '--distance', nargs = '?', default = 2, type = int, help = "maximum distance between region and exons edges allowed");
args = parser.parse_args();




class AnnotatedChimera(object):
	def __init__(self, intervals, intersections):
		self.intervals = intervals
		self.intersections = intersections
		
		
	def _intra(self, i1, i2, distance):
		distances1 = abs(self.intervals[0].start - i1.start), abs(self.intervals[0].stop - i1.stop)
		distances2 = abs(self.intervals[1].start - i2.start), abs(self.intervals[1].stop - i2.stop)
		exon_number_difference = int(i2.score) - int(i1.score)
		#test, if it is a splice junction;
		if( (distances1[1]<=distance and distances2[0]<=distance) or (distances1[0]<=distance and distances2[1]<=distance) ):
			
			#test, if it is a linear splice junction
			#forward case
			if(exon_number_difference == 1 and distances1[1]<=distance and distances2[0]<=distance):
				return 'lsj'
			#backward case
			elif(exon_number_difference == -1 and distances1[0]<=distance and distances2[1]<=distance):
				return "lsj"
			
			#test, if it is a circular splice junction
			#one exon circle
			if(exon_number_difference == 0):
				return "csj"
			#forward case
			elif(exon_number_difference < 0 and distances1[1]<=distance and distances2[0]<=distance):
				return 'csj'
			#backward case
			elif(exon_number_difference > 0 and distances1[0]<=distance and distances2[1]<=distance):
				return "csj"			
		else:
			return "intra";
					
				
	def annotate(self, distance):
		order = ['csj', 'intra'];
		types = [];
		for exon_name, d in self.intersections.items():
			if(len(d) == 2 and all(d.values())):
				t = self._intra(d[0], d[1], distance)
				if(t=="lsj"):
					return t;
				else:	
					types.append(self._intra(d[0], d[1], distance));
		for t in order:
			if(t in types):
				return t;
				
		return 'inter'		

	
	def yield_annotated(self, distance, offset):
		for k, interval in self.intervals.items():
			interval.attrs['interaction'] = self.annotate(distance)
			yield interval[:offset];
			
			
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
chimeras_vs_exons = chimeras.intersect(exons, s=True, wao=True)


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
		#if (True or intervals[0].name == 'sg_01_1|0'):
			#for nk, nd in intersections.items():
				#sys.stderr.write(str(intersections.keys()) + "\n");
				#for nkk, ni in nd.items():
					#sys.stderr.write(str(nkk) + "\n");
					#sys.stderr.write(str(ni) + "\n");
					
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



























