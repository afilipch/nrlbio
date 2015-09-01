#! /usr/bin/python	
'''Script outputs presence of certain binding modes in the interactions'''
import interaction_lib;
import mir_lib;
import parsing;
import sys;
import os;
import argparse;
from collections import *;
from Bio.Seq import reverse_complement;


parser = argparse.ArgumentParser(description='Script outputs presence of certain binding modes in the interactions');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to interactions (output/interactions.bed)");
parser.add_argument('-o', '--output', nargs = '?', default = "rstatistics/modes2.tsv", type = str, help = "name of the output")
parser.add_argument('-mr', '--mir', nargs = '?', type = str, required = True, help = "path mirna in fasta format");
parser.add_argument('-a', '--add', nargs = '?', default = False, const = True, type = bool, help = "add modes annotation to the interacions.bed")
args = parser.parse_args();


def execute(interactions, seeds):
	r = {};
	for inter in interactions:	
		seed = seeds[inter.mirid.split(",")[0]];
		mode, variants = seed.find_once(inter.tseq[:])
		r[inter.iid] = mode;
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
	
def get_weighted(signal, control, order, weighted = True):
	ans = [];
	weight = 1;
	st = float(sum(signal.values()));
	ct = float(sum(control.values()));
	sa = 0;
	ca = 0;
	for i, k in enumerate(order):
		ans.append((k, signal[k]/st, control[k]*weight/ct))
		sa += signal[k];
		ca += control[k];
		if(weighted and i < len(order)-1):
			weight = ((st-sa)/st)/((ct-ca)/ct)
	return ans;
	
mirdict = parsing.fasta2dict(args.mir);

seeds = mir_lib.get_seeds(mirdict, start=1, end=7);	

	
interactions = interaction_lib.get(args.path[0], undef = True);



il = execute(interactions, seeds);
cl = execute(control, seeds)

ic, cc = Counter(il.values()), Counter(cl.values());

order = ['match', 'mismatch', 'insertion', 'mismatch2', 'deletion', 'no_match']
ans = get_weighted(ic, cc, order, weighted = True);



	
	

o = open(args.output, 'w');
for el in ans:
	o.write("%s\t%1.4f\t%1.4f\n" % el);
o.close()	


	
if(args.add):
	f = open(args.path[0], 'r');
	a = [];
	for l in f:
		iid = l.strip().split("\t")[3]
		m = il.get(iid, "not_assessed")
		a.append(l.strip() + "\t" + m + "\n");
	f.close()
	
	f = open(args.path[0], 'w');
	for s in a:
		f.write(s)
	f.close()
	
#>>>>> first nt bias
nf = first(interactions);
norma = float(sum(nf.values()))
f = open("rstatistics/first_nucleotide.tsv", 'w')
for nucl in "ACTG":
	f.write("%s\t%d\t%1.3f\n"  % (nucl, nf[nucl], nf[nucl]/norma));
f.close();	
	
	
