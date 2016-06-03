#! /usr/bin/python
'''Assignes to a circle number of reads supporting outer splice junctions'''


import sys;
import argparse
from collections import defaultdict

from pybedtools import BedTool, Interval

#from nrlbio.pybedtools_extension import construct_gff_interval(chrom, start, stop, feature, score='0', strand='.', source='un', frame='.', attrs=[])


parser = argparse.ArgumentParser(description='Assignes to a circle number of reads supporting outer splice junctions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the circles, gff format");
parser.add_argument('--lsj', nargs = '?', required = True, type = str, help = "Path to the linear splice junction expression file (collapsed.lsj.gff), gff format");
args = parser.parse_args();

left_support = defaultdict(float)
right_support = defaultdict(float)

for interval in BedTool(args.lsj):
	left_support[interval.end] += float(interval.attrs['n_uniq'])
	right_support[interval.start] += float(interval.attrs['n_uniq'])
	
for interval in BedTool(args.path):
	interval.attrs['left_lsj'] = str(left_support[interval.start])
	interval.attrs['right_lsj'] = str(right_support[interval.end])
	sys.stdout.write(str(interval))
	