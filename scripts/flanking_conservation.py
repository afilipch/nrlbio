#! /usr/lib/python
'''Compares conservation of selected regions to the conservation of flanking regions. Assignes differential conservation score to the intervals corresponding to the conservation blocks'''

import argparse
import sys
from collections import defaultdict

from pybedtools import BedTool

from nrlbio.generators import targets_generator




parser = argparse.ArgumentParser(description='Compares conservation of selected regions to the conservation of flanking regions. Assignes differential conservation score to the intervals corresponding to the conservation blocks');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the intervals, gff format");
parser.add_argument('--mafasta', nargs = '?', required = True, type = str, help = "Path to the conservation blocks corresponding to the intervals, fasta file");
parser.add_argument('-l', '--left', nargs = '?', required = True, type = int, help = "Length of the left flank");
parser.add_argument('-r', '--right', nargs = '?', required = True, type = int, help = "Length of the right flank");
parser.add_argument('--refspecie', nargs = '?', required=True, type = str, help = "Name of the reference specie genome [hg19, mm10,..]"); 
args = parser.parse_args();

def get_cscore(reference, alseqs, positions):
	score = 0;
	for pos, (rn, block) in enumerate(zip(reference, zip(*alseqs))):
		if(pos in positions):
			for an in block:
				if(an == rn):
					score += 1;
	return float(score)/len(positions);


name2score = {}

for name, stub, targets in targets_generator(args.mafasta):
	target = targets.pop(args.refspecie, None)
	tl = len(target.replace("-", ''))
	rightbound = tl - args.right;
	
	flanks = []
	region = []
	onwogaps = 0
	for p, n in enumerate(target):
		if(n != '-'):
			if(onwogaps<args.left or onwogaps>=rightbound):
				flanks.append(p);
			else:
				region.append(p);
			onwogaps += 1;
	
	rscore = get_cscore(target, targets.values(), set(region))
	fscore = get_cscore(target, targets.values(), set(flanks))
	if(fscore):
		diffscore = rscore/fscore 
	else:
		diffscore = rscore
		
	name2score[name] = diffscore;
		
	#if(diffscore> 1.5):
		#print target
		#for v in targets.values():
			#print v
		#print
		#print
	
	
for interval in BedTool(args.path):
	interval.attrs['flankcons'] = "%1.2f" % name2score[interval.name]
	sys.stdout.write(str(interval))
	
	
	
			
