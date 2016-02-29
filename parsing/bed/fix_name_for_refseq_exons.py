#! /usr/lib/python
'''Changes the names of refseq exons to the appropriate (for splicing annotation) ones''' 
import argparse
import sys;

from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Changes the names of refseq exons to the appropriate (for splicing annotation) ones');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the bed file");
args = parser.parse_args()


for interval in BedTool(args.path):
	a = interval.name.split("_")
	interval.score = a[3];
	interval.name = "_".join(a[:2]);
	sys.stdout.write(str(interval));