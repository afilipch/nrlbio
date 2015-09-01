'''explore sequence complementarity'''
import interaction_lib;
import sys;
import os;
import argparse;
import copy
from collections import *;
from Bio.Seq import reverse_complement;


parser = argparse.ArgumentParser(description='explore sequence complementarity');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to interactions (output/interactions.bed)");
parser.add_argument('--minl', nargs = '?', default = 4, type = int, help = "min length of window")
parser.add_argument('--maxl', nargs = '?', default = 11, type = int, help = "max length of window")
args = parser.parse_args();

def run(interactions, minl, maxl):
	r = defaultdict(int)
	for inter in interactions:
		rs = reverse_complement(inter.tseq)
		ms = inter.mirseq
		for l in range(minl, maxl):
			for i in range(len(ms)-l+1):
				m = ms[i:i+l]
				if(m in rs):
					r[(i, i+l)] += 1;
	return r;
	
	
def run_control(interactions, minl, maxl, trials = 10):
	r = defaultdict(int)
	for stub in range(10):	
		cinteractions = interaction_lib.shufseq(interactions);
		for inter in cinteractions:
			rs = reverse_complement(inter.tseq)
			ms = inter.mirseq
			for l in range(minl, maxl):
				for i in range(len(ms)-l+1):
					m = ms[i:i+l]
					if(m in rs):
						r[(i, i+l)] += 1;
	for k, v in r.iteritems():
		r[k] = v/trials	
	return r;		
		
interactions = interaction_lib.get(args.path[0]);

signal = run(interactions, args.minl, args.maxl)
control = run_control(interactions, args.minl, args.maxl, trials = 10)


for l in range(args.minl, args.maxl):
	for i in range(max([x[0] for x in signal.keys()])):
		print "length %d\tstart %d\tend %d\tsignal %d\tcontrol %d" % (l, i, l+i, signal[(i, i+l)], control[(i, i+l)])
	print	


