# /usr/bin/python
'''Sorts cells according to the expression for each gene'''

import sys
import os
import copy
import argparse;
from collections import defaultdict
from math import log;

import numpy as np



parser = argparse.ArgumentParser(description='assignes hybridization energy to interactions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the digital expression matrix(rows are genes, collumns are cells), tsv format");
args = parser.parse_args();



def get_centers(expression):
	res = {};
	curexpr = expression[0][1];
	curlist = [expression[0][0]]
	curbase = -1;
	for base, (num, expr) in enumerate(expression[1:]):
		if(expr==curexpr):
			curlist.append(num);
		else:
			center = curbase + (len(curlist)+1.0)/2;
			for el in curlist:
				res[el] = center;
			curbase = base;
			curlist = [num];
			curexpr = expr
	else:
		center = curbase + (len(curlist)+1.0)/2;
		for el in curlist:
			res[el] = center;
			
	return res;


genes = {};

with open(args.path) as f:
	cellnames = f.next().strip().split("\t")
	for l in f:
		a = l.strip().split("\t");
		genename = a[0];
		expression = list(enumerate([float(x) for x in a[1:]]))
		expression.sort(key = lambda x: x[1])
		#cell2expression = dict(expression)
		cell2center = get_centers(expression);
		#print cell2center
		print "\t".join([str(cell2center[x]) for x in range(len(cell2center))])

		
		
		
		
		
		

