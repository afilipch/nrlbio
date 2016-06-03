#! /usr/bin/python
'''Compares expression of small RNAs between two conditions (replicates are allowed)'''
import argparse
import sys;
from collections import defaultdict
from math import log

import numpy
import pysam;
from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Compares expression of small RNAs between two conditions (replicates are allowed)');
parser.add_argument('-c1', '--condition1', nargs = '+', required = True, type = str, help = "Path to the small RNA expression files related to the replicates of the first experiment (condition), gff file");
parser.add_argument('-c2', '--condition2', nargs = '+', required = True, type = str, help = "Path to the small RNA expression files related to the replicates of the second experiment (condition), gff file");
parser.add_argument('--names', nargs = 2, default = ('condition_1', 'condition_2'), type = str, help = "Descriptive names for the first and second conditions");
parser.add_argument('--minexpr', nargs = '?', default = 0.0, type = float, help = "Min normalized expression threshold for sRNA to be analysed");
args = parser.parse_args();

size1 = len(args.condition1)
size2 = len(args.condition2)


def compexpr(c1, c2, varcoeff):
	mean1 = numpy.mean(c1);
	mean2 = numpy.mean(c2);
	diff =  abs(mean1-mean2)/numpy.mean(numpy.concatenate((c1,c2)))
	vdiff = diff/varcoeff
	inner_diff = max([(numpy.std(c1)/mean1)/varcoeff, (numpy.std(c2)/mean2)/varcoeff])
	lfc = log((mean1+0.1)/(mean2+0.1), 2) 
	return list(c1) + list(c2) + [inner_diff, vdiff, lfc]


#Get actual normalized expression for two conditions(with replicates)
mirids = set()
condition1 = defaultdict(lambda: defaultdict(float))
for en, cpath in enumerate(args.condition1):
	for interval in BedTool(cpath):
		expr = float(interval.attrs['norm_expr'])
		mirids.add(interval.attrs['Name'])
		condition1[interval.attrs['Name']][en] = expr;
		
condition2 = defaultdict(lambda: defaultdict(float))
for en, cpath in enumerate(args.condition2):
	for interval in BedTool(cpath):
		expr = float(interval.attrs['norm_expr'])
		mirids.add(interval.attrs['Name'])
		condition2[interval.attrs['Name']][en] = expr;
		
		
#Get global normalized standard deviation (coefficient of variance) on basis of replicates
cvs = [];
for mirid, edict in condition1.items():
	if(len(edict)>1):
		a = numpy.array(edict.values())
		if(max(a)>=args.minexpr):
			cvs.append(numpy.std(a)/numpy.mean(a))
		
for mirid, edict in condition2.items():
	if(len(edict)>1):
		a = numpy.array(edict.values())
		if(max(a)>=args.minexpr):
			cvs.append(numpy.std(a)/numpy.mean(a))
		
		
varcoeff = numpy.mean(numpy.array(cvs));
sys.stderr.write("Global variation coefficient is equal to %1.4f\n" % varcoeff)
		
		
		
#Calculate log fold change and zscore between the two means(between conditions)
expression = []
for mirid in mirids:
	cd1 = condition1[mirid]
	cd2 = condition2[mirid]
	c1 = numpy.array([cd1[x] for x in range(size1)])
	c2 = numpy.array([cd2[x] for x in range(size2)])
	if(max(c1)>args.minexpr or max(c2)>args.minexpr):
		expression.append([mirid] + compexpr(c1, c2, varcoeff))


#Output the results
header = ['id'] + ["%s_%d" % (args.names[0], x+1) for x in range(size1)] + ["%s_%d" % (args.names[1], x+1) for x in range(size2)] + ['zscore_inner', 'zscore', 'lfc2']
sys.stdout.write("\t".join(header) + "\n") 

expression.sort(key = lambda x: x[-1], reverse=True)
for l in expression:
	sys.stdout.write("\t".join([str(x) for x in l]) + "\n")
	






