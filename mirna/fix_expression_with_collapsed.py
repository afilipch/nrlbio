#! /usr/bin/python
'''Fixes small RNA expression using the information on nonunique mappings'''
import argparse
import sys;
from collections import Counter, defaultdict

from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Fixes small RNA expression using the information on nonunique mappings');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the expression file, gff format");
parser.add_argument('--collapsed', nargs = '?', required = True, type = str, help = "Path to the auxillary file with collapsed nonunique mappings");
args = parser.parse_args();

identicals = defaultdict(list)

for interval in BedTool(args.collapsed):
	identicals[interval.name].append(interval.chrom)
	
adjustments = defaultdict(float)

for k, mirids in identicals.iteritems():
	adjustment = 1.0/len(mirids)
	adjustments[mirids[0]] += (adjustment-1)
	for mirid in mirids[1:]:
		adjustments[mirid] += adjustment
		
		
totalreads = sum([float(x.attrs['raw_expr']) for x in BedTool(args.path)])/1000000.0

for interval in BedTool(args.path):
	interval.attrs['raw_expr'] = str(float(interval.attrs['raw_expr']) + adjustments[interval.name])
	interval.attrs['norm_expr'] = '%1.2f' % (float(interval.attrs['raw_expr'])/totalreads)
	sys.stdout.write(str(interval))
	
		

