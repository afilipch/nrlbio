#! /usr/bin/python	
'''Script outputs binding energy of given interactions'''
import interaction_lib;
import rnahybrid;
import sys;
import os;
import argparse;
import miscellaneous;
import parsing
from collections import *;


parser = argparse.ArgumentParser(description='Script outputs binding energy of given interactions');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to interactions (output/interactions.bed)");
parser.add_argument('-a', '--add', nargs = '?', default = False, const = True, type = bool, help = "add modes annotation to the interacions.bed")
args = parser.parse_args();

order = ["seed27", "seed38", "wobble", "mismatch28", "insertion28", "deletion28", "no mode"]

def execute(interactions):
	pattern = defaultdict(int);
	energies = {}
	for inter in interactions:		
		p, e = rnahybrid.run(inter.tseq, inter.mirseq[:]);
		energies[inter.iid] = e;
		for i in range(len(p)):
			pattern[i] += p[i];
	return pattern, energies;	
	
def execute_c(interactions):
	pattern = defaultdict(int);
	energies = defaultdict(int);
	for j in range(20):
		for inter in interactions:
			u = inter.mirseq[:];
			t = parsing.shuffle_string(inter.tseq)
			p, e = rnahybrid.run(t, u);
			energies[inter.iid] += e;
			for i in range(len(p)):
				pattern[i] += p[i];
	for k, v in pattern.iteritems():
		pattern[k] = v/20;
	for k, v in energies.iteritems():
		energies[k] = v/20;			
	return pattern, energies;	

interactions = interaction_lib.get(args.path[0]);
shufseq = interaction_lib.shufseq(interactions);
shufpairs = interaction_lib.shufpairs(interactions);
t = float(len(interactions));


ip, iel = execute(interactions);
sp, sel = execute(shufseq)
pp, pel = execute(shufpairs);

ie, se , pe = Counter(iel.values()), Counter(sel.values()), Counter(pel.values()),


o = open(os.path.join("rstatistics", "energy.tsv"), 'w');
for e in sorted(list(set(ie.keys() + se.keys() + pe.keys()))):
	o.write("%1.1f\t%1.6f\t%1.6f\t%1.6f\n" % (e, ie[e]/t, se[e]/t, pe[e]/t));
o.close()	
	
o = open(os.path.join("rstatistics", "pattern.tsv"), 'w');
for i in sorted(list(set(ip.keys() + sp.keys() + pp.keys()))):
	o.write("%d\t%1.3f\t%1.3f\t%1.3f\n" % (i, ip[i]/t, sp[i]/t, pp[i]/t));
o.close()


if(args.add):
	f = open(args.path[0], 'r');
	a = [];
	for l in f:
		iid = l.strip().split("\t")[3]
		m = iel.get(iid, 0)
		a.append(l.strip() + "\t" + str(m) + "\n");
	f.close()
	
	f = open(args.path[0], 'w');
	for s in a:
		f.write(s)
	f.close()





