'''converts bed-like file of chimeric reads into interactions. That is merging chimeras with intersecting regions'''

import sys
import argparse;
from collections import defaultdict
from operator import itemgetter, attrgetter

from pybedtools import BedTool

from nrlbio.pybedtools_extension import doublebed2dict, generate_overlaping_intervals
from nrlbio.interaction import Interaction

parser = argparse.ArgumentParser(description='converts bed-like file of chimeric reads into interactions. That is merging chimeras with intersecting regions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to chimeras bed-like file");
parser.add_argument('-d', '--distance', nargs = '?', default = -12, type = int, help = "minimum overlap(negative number)/maximum distance(positive number) in nucleotides to merge intervals");
parser.add_argument('-n', '--name', nargs = '?', required = True, type = str, help = "name for interactions")
parser.add_argument('-oi', '--interactions', nargs = '?', required = True, type = str, help = "path to output interactions file")
parser.add_argument('-od', '--dictionary', nargs = '?', required = True, type = str, help = "path to output \"interaction to read id\" file")
parser.add_argument('--order', nargs = '?', default = False, const=True, type = int, help = "keeps order of left to right parts in interactions");
args = parser.parse_args();

od = open(args.dictionary, 'w');
oi = open(args.interactions, 'w');


def intervals2interaction(intervals, distance, number, order=False):
	if(order):
		merged_regions = [[], []];
		for interval in intervals:
			#sys.stderr.write("%d" % int(interval.name.split("|")[1]))
			merged_regions[int(interval.name.split("|")[1])].append(interval)			
	else:
		intervals.sort(key = attrgetter('chrom','strand','start'));
		merged_regions = list(generate_overlaping_intervals(intervals, distance));

	if(len(merged_regions)==2):
		interaction = Interaction.from_intervals("%s_%d" % (args.name, number), merged_regions)
		oi.write("%s\n" % interaction.doublebed())
		od.write("%s\t%s\n" % (interaction.name, ",".join(interaction.read_names)))
	else:
		sys.stderr.write("Warning: interacting regions may be further split or merged")
		for i in intervals:
			sys.stderr.write(str(i))
		sys.stderr.write("_"*140	+ "\n")	
			
	


bed = BedTool(args.path);
	
name2intervals = doublebed2dict(bed);
name2interaction = defaultdict(list);
interaction2name = defaultdict(list);

for c, i in enumerate(bed.merge(s=True, d=args.distance, c='4,6', o='distinct', delim=';')):
	 for name in i.name.split(";"):
		 n = name.split("|")[0]
		 name2interaction[n].append(c);
		 
for k, v in	name2interaction.iteritems():
	key = tuple(sorted(v));
	interaction2name[key].append(k)
		 
number_chimeras = len(name2intervals);
number_interactions = 0
number_self_intersection = 0;
		
for number, (k, v) in enumerate(interaction2name.items()):
	if(k[0]!=k[1]):
		intervals = []
		for name in v:
			for i in name2intervals[name]:
				intervals.append(i)
		intervals2interaction(intervals, args.distance, number+1, order=args.order);
		number_interactions += 1;
	else:
		number_self_intersection += 1;
		
		
sys.stderr.write("number of chimeras\t%d\ninteractions generated\t%d\nself interacting removed\t%d\n" % (number_chimeras, number_interactions, number_self_intersection));

