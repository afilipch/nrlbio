#! /usr/lib/python
'''Investigates some peculiarities of miRNA:miRNA interactions'''

import argparse
import os
import sys
from collections import defaultdict

from pybedtools import BedTool, Interval



parser = argparse.ArgumentParser(description='Investigates some peculiarities of miRNA:miRNA interactions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the miRNA precursors intersted with miRNA binding sites (mirmir_targets_vs_interactions.py output), gff format");
#parser.add_argument("--mirexpr", nargs = '?', type=str, required=True, help="Path to the miRNA expression, gff format")
parser.add_argument("--minfraction", nargs = '?', type=float, default=0.5, help="Minimal required overlap between a binding site and a mature miRNA on a precursor as a fraction of the binding site length")
args = parser.parse_args();

def readsite(s):
	if(s):
		a = s.split(':');
		return a[0], int(a[1]), int(a[2])
	else:
		return None
	
	
def check_overlap(bsite, arm, minfraction):
	length = float(bsite[2]-bsite[1]);
	overlap = min(bsite[2], arm[2]) - max(bsite[1], arm[1]);
	return overlap/length > minfraction

def check_dicer_duplex(name1, name2):
	a1 = name1.split('-');
	a2 = name2.split('-');
	if(tuple(a1[:3]) == tuple(a2[:3])):
		if(len(set([a1[-1], a2[-1]])) == 2):
			return True;
		else:
			return False;
	else:
		return False;
	
	
interactors = defaultdict(list)

for interval in BedTool(args.path):
	n_uniq = float(interval.attrs['n_uniq'])
	bsite = readsite(interval.attrs['bsite'])
	arms = [];
	for s in ['arm1', 'arm2']:
		arm = readsite(interval.attrs[s])
		if(arm):
			arms.append(arm);
			if(check_overlap(bsite, arm, args.minfraction) and not check_dicer_duplex(bsite[0], arm[0])):
				print '%s\t%s\t%d\t%s' % (bsite[0], arm[0], n_uniq, check_dicer_duplex(bsite[0], arm[0]))
				#print interval
		
