#! /usr/bin/python	
'''Script outputs presence of certain binding modes in the interactions'''
import interaction_lib;
import mir_lib;
import sys;
import os;
import argparse;
from collections import *;
from Bio.Seq import reverse_complement;


parser = argparse.ArgumentParser(description='Script outputs presence of certain binding modes in the interactions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to interactions (output/interactions.bed)");
parser.add_argument('-o', '--output', nargs = '?', default = "rstatistics/modes.tsv", type = str, help = "name of the output")
parser.add_argument('-a', '--add', nargs = '?', default = False, const = True, type = bool, help = "add modes annotation to the interacions.bed")
args = parser.parse_args();

order = ["seed27", "seed38", "wobble", "mismatch28", "insertion28", "deletion28", "no mode"]

def execute(interactions):
	r = {};
	for inter in interactions:	
		m = mir_lib.modes(inter.tseq, inter.mirseq);
		r[inter.iid] = m
	return r;	
	
def first(interactions):
	a = defaultdict(int);
	for inter in interactions:
		s27 = reverse_complement(inter.mirseq[1:7]);
		pos = inter.tseq.find(s27);
		if(pos > -1 and pos+6 < len(inter.tseq)):
			nucl = inter.tseq[pos+6]
			a[nucl] += 1
	return a;		
	
	
interactions = interaction_lib.get(args.path, undef = True);
control = interaction_lib.shufseq(interactions);
t = float(len(interactions));


il = execute(interactions);
cl = execute(control)

ic, cc = Counter(il.values()), Counter(cl.values());

o = open(args.output, 'w');
for i in [0,1,3,2,4,5,-1]:
	o.write("%s\t%1.4f\t%1.4f\n" % (order[i], ic[i]/t, cc[i]/t));
o.close()	
	
if(args.add):
	name = ".".join(args.path.split(".")[:-1] + ["modes", args.path.split(".")[-1]]);
	with open(args.path, 'r') as inp, open(name, 'w') as out:
		for l in inp:
			a = l.strip().split("\t")
			m = str(il.get(a[3], -2)+1)
			out.write("\t".join(a[:13] + [m]) + "\n")
	
#>>>>> first nt bias
nf = first(interactions);
norma = float(sum(nf.values()))
f = open("rstatistics/first_nucleotide.tsv", 'w')
for nucl in "ACTG":
	f.write("%s\t%d\t%1.3f\n"  % (nucl, nf[nucl], nf[nucl]/norma));
f.close();	
	
	
