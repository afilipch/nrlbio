#! /usr/lib/python
'''Changes name for each bed/gff interval according to the provided translational table'''
import argparse
import sys;

from pybedtools import BedTool, Interval



parser = argparse.ArgumentParser(description='Changes name for each bed/gff interval according to the provided translational table');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the bed/gff file");
parser.add_argument('--table', nargs = '?', default = True, type = str, help = "Path to the translational table, tsv file. First collumn - old names, second column - new names");
args = parser.parse_args();

old2new = {};
with open(args.table) as f:
	for l in f:
		a = l.strip().split("\t")
		old2new[a[0]] = a[1]

for interval in BedTool(args.path):
	newname = old2new.get(interval.name, interval.name)
	#if(newname, interval.name):
		#print newname, interval.name
	interval.name = newname
	sys.stdout.write(str(interval))