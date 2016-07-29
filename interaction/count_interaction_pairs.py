#! /usr/bin/python
'''Counts number of miRNA:gene pairs''' 
import argparse
import sys;
from collections import defaultdict

from pybedtools import BedTool, Interval

parser = argparse.ArgumentParser(description='Counts number of miRNA:gene pairs');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated interactions, gff format");
args = parser.parse_args();

pairs_targets = defaultdict(float);
pairs_reads = defaultdict(float);

for interval in BedTool(args.path):
	tg = interval.attrs['gene_symbols']
	if(tg):
		nuniq = float(interval.attrs['n_uniq'])
		mirid = interval.attrs['mirid']
		gsymbols = tg.split(":")
		norm = float(len(gsymbols))
		for g in gsymbols:
			pairs_targets[(mirid,g)] += 1/norm
			pairs_reads[(mirid,g)] += nuniq/norm
			
			
pairs = sorted(pairs_reads.items(), key = lambda x: x[1], reverse=False)
for pair, value in pairs:
	print "%s\t%s\t%1.2f\t%1.2f" % (pair[0], pair[1], value, pairs_targets[pair]) 