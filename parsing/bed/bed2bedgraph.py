#! /usr/lib/python
'''Converts bed file with scores into bedgraph files''' 
import argparse
import sys;
from collections import defaultdict

from pybedtools import BedTool, Interval;

from nrlbio.itertools_extension import smooth_values

parser = argparse.ArgumentParser(description='Converts bed file with score into bedgraph files');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the bed file");
parser.add_argument('-s', '--smooth', nargs = '?', default = 0, type = int, help = "If set, scores will be aggregated in [--smooth]-length windows");
args = parser.parse_args()

def smooth(interval, window_size):
	return smooth_values((float(x) for x in interval[-1].split(",")), window_size)

def no_smooth(interval, window_size):
	return (float(x) for x in interval[-1].split(","))

if(args.smooth):
	fun = smooth;
else:
	fun = no_smooth


curstart = 0;
curscore = ''
for interval in BedTool(args.path):
	for pos, score in enumerate(fun(interval, args.smooth)):
		if(pos and score!=curscore):
			print "%s\t%d\t%d\t%1.4f" % (interval.name, curstart, pos, curscore);
			curstart = pos;
		curscore=score;
	else:
		print "%s\t%d\t%d\t%1.4f" % (interval.name, curstart, pos, curscore);