# /usr/bin/python
'''reads already constructed 2d grid and tests LRGFDR'''

import argparse
from collections import defaultdict
from random import randint

from nrlbio import LRGFDR;
from visualize_LRGFDR import  grid2heatmap, animate_clustering

parser = argparse.ArgumentParser(description='reads already constructed 2d grid and tests LRGFDR');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the grid");
parser.add_argument('-v', '--video', nargs = '?', type = str, default = False, help = "path to the video output");
args = parser.parse_args();

grid = LRGFDR.Grid.from_csv(args.path, delimiter="\t");
clusters, nclusters = LRGFDR.generate_clusters(grid, support = 0.01, maxiter = 10, nciter=1,  fdr=0.05, lookforward=10, fit_function=LRGFDR.ff_balanced)



def ff_balanced(x):
	maxfdr = 0.05;
	if(x.real):
		return x.real*(maxfdr - x.imag/(x.imag+x.real+0.01))
	#else:
		#return maxfdr-1
#grid2heatmap(grid, ff_balanced)
if(args.video):
	animate_clustering(grid, ff_balanced, clusters, nclusters, args.video)