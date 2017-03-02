#! /usr/bin/python
'''Intersects interactions with other interactions. Reports the consensus for the first'''


import sys;
import argparse

from pybedtools import BedTool

from nrlbio.generators import generator_doublebed


parser = argparse.ArgumentParser(description='Intersects interactions with a given bed/gff file and outputs only those interactions which have or not have an overlap with provided genomic regions');
parser.add_argument('path', metavar = 'N', nargs = 2, type = str, help = "Path to the interactions file1 interactions file2, single gff file (targets.gff)");
#parser.add_argument('--intervals', nargs = '?', type = str, help = "path to the intervals to intersect with, bed/gff file");
#parser.add_argument('--invert', '-v', nargs = '?', const=True, default=False, type = bool, help = "If set interactions with NO overlap will be output");
args = parser.parse_args();

OFFSET = 9

interactions1 = BedTool(args.path[0])
interactions2 = BedTool(args.path[1])

def get_mirna(interval, offset):
	for l in interval[offset].split(';'):
		if(l.startswith('mirid')):
			return l.split('=')[1];
	
	
passed = []
for interval in interactions1.intersect(interactions2, wo=True, s=True, f=0.95):
	samemirna = get_mirna(interval, OFFSET-1) == get_mirna(interval, OFFSET*2-1)
	if(samemirna):
		passed.append(interval.name)

print passed
print len(passed)
print len(list(set(passed)))

#if(args.invert):
	#for i1, i2 in generator_doublebed(args.path):
		#cid = i1.attrs['ID'].split('|')[0]
		#if(cid not in passed):
			#sys.stdout.write("%s%s" % (str(i1), str(i2)))
#else:
	#for i1, i2 in generator_doublebed(args.path):
		#cid = i1.attrs['ID'].split('|')[0]
		#if(cid in passed):
			#sys.stdout.write("%s%s" % (str(i1), str(i2)))	




