#! /usr/bin/python
'''Intersects interactions with a given bed/gff file and outputs only those interactions which have or an overlap with provided genomic regions'''


import sys;
import argparse

from pybedtools import BedTool

from nrlbio.generators import generator_doublebed


parser = argparse.ArgumentParser(description='Intersects interactions with a given bed/gff file and outputs only those interactions which have or not have an overlap with provided genomic regions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the interactions, double bed/gff file");
parser.add_argument('--intervals', nargs = '?', type = str, help = "path to the intervals to intersect with, bed/gff file");
parser.add_argument('--invert', '-v', nargs = '?', const=True, default=False, type = bool, help = "If set interactions with NO overlap will be output");
args = parser.parse_args();


intervals = BedTool(args.intervals)
interactions = BedTool(args.path)
	
	
passed = set();
for interval in interactions.intersect(intervals, u=True, s=True):
	passed.add(interval.attrs['ID'].split('|')[0]);

if(args.invert):
	for i1, i2 in generator_doublebed(args.path):
		cid = i1.attrs['ID'].split('|')[0]
		if(cid not in passed):
			sys.stdout.write("%s%s" % (str(i1), str(i2)))
else:
	for i1, i2 in generator_doublebed(args.path):
		cid = i1.attrs['ID'].split('|')[0]
		if(cid in passed):
			sys.stdout.write("%s%s" % (str(i1), str(i2)))	




#if(args.invert):
	#intersection = interactions.intersect(intervals, v=True, s=True);
#else:
	#intersection = interactions.intersect(intervals, u=True, s=True);
	
#previous = None;
#pid = ''
#for interval in intersection:
	#cid = interval.attrs['ID'].split('|')[0]
	#if(cid==pid):
		#sys.stdout.write("%s%s" % (str(previous), str(interval)))
	#previous = interval
	#pid = cid
	
