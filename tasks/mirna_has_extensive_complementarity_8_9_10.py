#! /usr/bin/python
'''Script answers to the question: If miRNA:target pair bound via 2-7(with or without 1mm) has an extensive complementarity at positions 8,9,10'''
import sys;
import argparse
from collections import defaultdict

from nrlbio.chimirna import interaction_lib;
from nrlbio.mirna import Mirna;
from nrlbio.sequencetools import RevComplDict;
from nrlbio.formatting import  feature_dict_total


parser = argparse.ArgumentParser(description='Script answers to the question: If miRNA:target pair bound via 2-7(with or without 1mm) has an extensive complementarity at positions 8,9,10');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to interactions (output/interactions.bed)");
parser.add_argument('--seed_start', nargs = '?', default = 1, type = int, help = "start postion of mirna seed, 0-based and inclusive")
parser.add_argument('--seed_stop', nargs = '?', default = 7, type = int, help = "stop postion of mirna seed, 0-based and exclusive")
parser.add_argument('--steps', nargs = '?', default = 4, type = int, help = "controls how far seed pairing will be extended downstream miRNA")
parser.add_argument('--mismatch', nargs = '?', default = False, const = True, type = bool, help = "check for 1mm seed matches as well")
args = parser.parse_args();

def extend_match(mirseq, tseq, pos, steps):
	for step in range(steps):
		npos = pos-step-1
		if(npos < 0):
			return -1;
		elif(tseq[npos] != RevComplDict[mirseq[args.seed_stop+step]]):
			return step;
		else:
			pass;
	return steps;
	
	
def check(mirna, tseq, steps = 4, mismatch=False):
	pos = tseq.find(mirna.match);
	if(pos == -1):
		if(mismatch):
			mirna.set_1mm_matches()
			for match in mirna.matches_1mm:
				pos = tseq.find(match)
				if(pos != -1):
					return extend_match(mirna.seq, tseq, pos, steps);
			else:
				return -1;
		else:
			return -1;
	else:
		return extend_match(mirna.seq, tseq, pos, steps);
		

interactions = interaction_lib.get(args.path, undef = False);
extensions = defaultdict(int)
control_extension = defaultdict(int)

for interaction in interactions:
	mirna = Mirna(interaction.mirid, interaction.mirseq, seed_start = args.seed_start, seed_stop = args.seed_stop)
	extensions[check(mirna, interaction.tseq, steps = args.steps, mismatch=args.mismatch)] += 1;
	#control_extensions[check(mirna, interaction.tseq, steps = args.steps, mismatch=args.mismatch)] += 1;
	
print feature_dict_total(extensions, top = 0, key_names = None)


