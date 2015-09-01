#! /usr/bin/python	
'''Script outputs preference of miRNA family members target the locus'''
import interaction_lib;
import mir_lib;
import sys;
import os;
import argparse;
import random;
from collections import *;


parser = argparse.ArgumentParser(description='Script outputs preference of miRNA family members target the locus');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to interactions (output/interactions.bed)");
parser.add_argument('-s', '--system', nargs = '?', required = True, type = str, help = "genome system")
parser.add_argument('--overlap', nargs = '?', default = 16, type = int, help = "minimum overlap of targets");
parser.add_argument('-o', '--output', nargs = '?', default = "rstatistics/family.tsv", type = str, help = "name of the output")
args = parser.parse_args();

def run(loci):
	seeds = defaultdict(int);
	for l in loci:
		l.multitargeting()
		if(l.onefam):
			for k, v in l.s2m.iteritems():
				if(len(v) > 1):
					seeds[k] += 1;
					
	return seeds;				

interactions = interaction_lib.get(args.path[0], undef = False);
print len(interactions)
loci = interaction_lib.interactions2loci(interactions, overlap = args.overlap)
true = run(loci);

#>>> control;
control = defaultdict(list);
for i in range(100):
	r = [];
	ms = [(x.mirid, x.mirseq) for x in interactions];
	for inter in interactions:
		index = random.randint(0,len(ms) - 1);
		r.append(ms.pop(index));
	for l in loci:
		new = [];
		for inter in l.interactions:
			m = r.pop()
			chromosome, start, end, iid, score, strand, mirid, mirseq, tseq, map_quality, totalreads, indreads, indseqs = list(inter);
			new.append(interaction_lib.Interaction(chromosome, start, end, iid, score, strand, m[0], m[1], tseq, map_quality, totalreads, indreads, indseqs))
			l.interactions = new;
	for k, v in run(loci).iteritems():
		control[k].append(v);
		
for k, v in control.iteritems():
	v.sort();

	

s2n = mir_lib.seeds2name(args.system)
print [len(x) for x in control.values()];
#print [len(x) for x in true.values()];

o = open(args.output, 'w');
for k,v in Counter(true).most_common(5):
	print len(control[k])
	o.write("%s\t%d\t%d\t%d\t%d\n"  % (s2n[k], v, int(sum(control[k])/float(len(control[k]))), control[k][int(len(control[k])*0.05)], control[k][int(len(control[k])*0.95)]));
o.close()	

#print run(loci)


#mdict = {};
#for inter in interactions:
	#mdict[inter.mirid] = inter.mirseq;