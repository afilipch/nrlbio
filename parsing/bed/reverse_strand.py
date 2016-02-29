#! /usr/lib/python
'''Script just changes strand to the opposite''' 
import argparse
import sys;

from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Script just changes strand to the opposite');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the bed file");
args = parser.parse_args()

d = {'-': '+', '+': '-'}

for interval in BedTool(args.path):
	interval.strand = d[interval.strand]
	sys.stdout.write(str(interval));