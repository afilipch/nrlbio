'''investigates shared sites between different miRNA'''
import interaction_lib;
import mir_lib;
import sys;
import os;
import argparse;
import random;
from collections import *;


parser = argparse.ArgumentParser(description='investigates shared sites between different miRNA');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to interactions (output/interactions.bed)");
parser.add_argument('-s', '--system', nargs = '?', required = True, type = str, help = "genome system")
parser.add_argument('--overlap', nargs = '?', default = 16, type = int, help = "minimum overlap of targets");
parser.add_argument('-o', '--output', nargs = '?', default = "rstatistics/family.tsv", type = str, help = "name of the output")
args = parser.parse_args();

def run(loci):
	seeds = defaultdict(int);
	for l in loci:
		l.multitargeting()
		models = set()
		for mset in l.mirids:
			for mirid in mset:
				models.add(mirid.split("-")[0])
		if(len(models) > 1):
			print l.represent();
					
			

interactions = interaction_lib.get(args.path[0]);
loci = interaction_lib.interactions2loci(interactions, overlap = args.overlap)
run(loci);

