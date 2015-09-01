#! /usr/bin/python	
'''Script outputs presence of crosslink relative to seed-matches with 1 mismatch'''
import interaction_lib;
import mir_lib;
import parsing;
import sys;
import os;
import argparse;
from collections import *;
from Bio.Seq import reverse_complement;


parser = argparse.ArgumentParser(description='Script outputs presence of certain binding modes in the interactions');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to interactions (final/interactions.bed)");
parser.add_argument('-mr', '--mir', nargs = '?', type = str, required = True, help = "path mirna in fasta format");
parser.add_argument('-r', '--reads', nargs = '?', type = str, required = True, help = "path to read2int file");
args = parser.parse_args();



def execute(interactions, seeds, int2crosslink):
	r = defaultdict(int);
	control = defaultdict(int);
	for inter in interactions:	
		seed = seeds[inter.mirid];
		if(seed.match in inter.tseq):
			continue;
		
		for mm in seed.mismatched[0]:
			if mm.seq in inter.tseq:
				pos = inter.tseq.find(mm.seq);
				break;
		else:
			continue
		
		
		
		for cr in int2crosslink[inter.iid]:
			if(cr == -1):
				continue;
			else:	
				r[cr-pos] += 1;
				for i, nucl in enumerate(inter.tseq):
					if (nucl == "T"):
						control[i-pos]+=1;
				
	return r, control;
	
int2crosslink = defaultdict(list);	
f = open(args.reads)
for l in f:
	a = l.strip().split();
	int2crosslink[a[1]].append(int(a[2]))
f.close()	
	
mirdict = parsing.fasta2dict(args.mir);
seeds = mir_lib.get_seeds(mirdict, start=1, end=7);		
interactions = interaction_lib.get(args.path[0]);

r, control = execute(interactions, seeds, int2crosslink);

#>>>>>>>>>>>>>> normalization

t1, t2 = float(sum(r.values())), float(sum(control.values()));

for i in range(-15, 11):
	print "%d\t%1.4f" % (i, (r[i]/t1)/((control[i]+0.01)/t2))



a = r.items()
a.sort(key = lambda x: x[1]);
