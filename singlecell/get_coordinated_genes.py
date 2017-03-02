# /usr/bin/python
'''Finds groups of coordinated genes'''

import sys
import os
import copy
import argparse;
from collections import defaultdict
from math import log;

import numpy as np



parser = argparse.ArgumentParser(description='Finds groups of coordinated genes');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the ddistances between the genes");
args = parser.parse_args();

distances = defaultdict(int)

with open(args.path) as f:
	for l in f:
		a = l.strip().split(" ")
		d = float(a[2])
		distances[int(d)] += 1;
		
		
for k in sorted(distances.keys()):
	print "%d\t%d" %  (k, distances[k])
		
#import matplotlib.pyplot as plt
#boundaries = range(2000)
#weights = [distances[x] for x in boundaries]
#ax.hist(boundaries, weights=weights, bins=bins, label=('0', '1', '2'), rwidth=1, color=colors)