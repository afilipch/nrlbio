# /usr/bin/python
'''Sorts cells according to the expression for each gene'''

import sys
import os
import copy
import argparse;
from collections import defaultdict
from itertools import combinations
from math import log;

import numpy as np



parser = argparse.ArgumentParser(description='assignes hybridization energy to interactions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the digital expression matrix(rows are genes, collumns are cells, entries are centers of the expressions), tsv format");
args = parser.parse_args();

def get_distance(g1, g2, size):
	return (sum((g1 - g2)**2)/size)**0.5



def get_genes(path):
	genes = [];
	with open(path) as f:
		for l in f:
			genes.append(np.array( [float(x) for x in l.strip().split("\t")] ));
	return genes



if __name__ == "__main__":
	genes = get_genes(args.path)
	size = len(genes[0]);
	for c, ((n1, g1), (n2, g2)) in enumerate(combinations(enumerate(genes), 2)):
		print n1, n2, get_distance(g1, g2, size);
		if(c and c % 10000 == 0):
			sys.stderr.write("%d processed\n" % c) 
	