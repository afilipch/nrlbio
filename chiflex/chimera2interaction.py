'''converts bed-like file of chimeric reads into interactions. That is merging chimeras with intersecting regions'''

import sys
import argparse;
from collections import defaultdict
from operator import itemgetter, attrgetter

from pybedtools import BedTool

from nrlbio.pybedtools_extension import doublebed2dict, generate_overlaping_intervals
from nrlbio.interaction import bed2interactions, Interaction

parser = argparse.ArgumentParser(description='converts bed-like file of chimeric reads into interactions. That is merging chimeras with intersecting regions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to chimeras bed-like file");
parser.add_argument('-d', '--distance', nargs = '?', default = -12, type = int, help = "minimum overlap(negative number)/maximum distance(positive number) in nucleotides to merge intervals");
parser.add_argument('-n', '--name', nargs = '?', required = True, type = str, help = "name for interactions")
parser.add_argument('-oi', '--interactions', nargs = '?', required = True, type = str, help = "path to output interactions file")
parser.add_argument('-od', '--dictionary', nargs = '?', required = True, type = str, help = "path to output \"interaction to read id\" file")
args = parser.parse_args();

od = open(args.dictionary, 'w');
oi = open(args.interactions, 'w');


def intervals2interaction(intervals, distance, number):
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
		print "_"*140		
			
	


bed = BedTool(args.path);
	
name2intervals = doublebed2dict(bed);
name2interaction = defaultdict(list);
interaction2name = defaultdict(list);

for c, i in enumerate(bed.merge(s=True, nms=True, d=args.distance)):
	 for name in i.name.split(";"):
		 name2interaction[name].append(c);
		 
for k, v in	name2interaction.iteritems():
	key = tuple(sorted(v));
	interaction2name[key].append(k)
		 
		
		
for number, (k, v) in enumerate(interaction2name.items()):
	if(k[0]!=k[1]):
		intervals = []
		for name in v:
			for i in name2intervals[name]:
				intervals.append(i)
		intervals2interaction(intervals, args.distance, number+1)
	

