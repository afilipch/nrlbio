#! /usr/bin/python
'''Compares relative expression of miRNAs and corresponding destructive sites for two different compartments'''


import sys;
import argparse
from collections import defaultdict
from itertools import izip
from math import log

from pybedtools import BedTool
import pysam
from scipy.stats.stats import pearsonr 


parser = argparse.ArgumentParser(description='Compares relative expression of miRNAs and corresponding destructive sites for two different compartments');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the mirna expression, tsv custom format");
parser.add_argument('-c1', '--condition1', nargs = '?', required = True, type = str, help = "Path to the RNA expression related to the first experiment (condition), gff file");
parser.add_argument('-c2', '--condition2', nargs = '?', required = True, type = str, help = "Path to the RNA expression related to the second experiment (condition), gff file");
args = parser.parse_args();

mir2expr = {};

with open(args.path) as f:
	f.next()
	for l in f:
		a = l.strip().split('\t')
		mir2expr[a[0]] = float(a[-1])
		
		

pairs = []
addition = 0.00001
for i1, i2 in izip(BedTool(args.condition1), BedTool(args.condition2)):
	score = float(i1.attrs['destructive_score'])
	if(score > 35):
		mirid = i1.attrs['mirid']
		e1 = float(i1.attrs['norm_coverage'])
		e2 = float(i2.attrs['norm_coverage'])
		if(e1+e2>0.005):
			lfc = log((e1+addition)/(e2+addition), 2)
			pairs.append((i1.chrom, i1.start, i1.end, i1.chrom, mirid, (e1+e2)/2, score, lfc, mir2expr[mirid]))
	
	
#mirlfc = [x[-1] for x in pairs]
#targetlfc = [x[-2] for x in pairs]
#print pearsonr(mirlfc, targetlfc)

pairs.sort(key = lambda x: x[-2]*x[-1]);
for p in pairs:
	print "\t".join([str(x) for x in p])





	
	
	