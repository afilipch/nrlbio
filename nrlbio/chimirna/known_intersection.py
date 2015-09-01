import random;
import sys;
import re;
import os;
import argparse;
import multiprocessing;
from collections import defaultdict;




parser = argparse.ArgumentParser(description='explore enrichment of seeds in targets');
parser.add_argument('path', metavar = 'N', nargs = 1, type = str, help = "path to mir2genes.file");
parser.add_argument('-k', '--known', nargs = '?', type = str, help = "path to known interactions from miRTarBase")
parser.add_argument('-d', '--dataset', nargs = '+', type = str, help = "name of dataset")
args = parser.parse_args();

dname = " ".join(args.dataset);

def execute(mg, known):
	r = defaultdict(set);
	for mirid, gene, pid in known:
		if(gene in mg[mirid]):
			r[mirid].add((gene,pid));	
	return r;		

real = defaultdict(set);

f = open(args.path[0]);
for l in f:
	a = l.strip().split("\t")
	for k in a[0].split(","):
		m = k.split("-");
		if(len(m) <= 3):			
			real["-".join(m)].add(a[1]);
		else:
			real["-".join(m[:-1])].add(a[1]);
f.close();


known = [];
f = open(args.known);
f.readline()
for l in f:
	a = l.strip().split("\t")
	m = a[1].split("-");
	if(len(m) <= 3):			
		mirid = "-".join(m)
	else:
		mirid = "-".join(m[:-1])	
	known.append((mirid, a[3], a[8]));
f.close();

r = execute(real, known);
for mirid, gset in r.iteritems():
	for gene, pid in gset:
		print "\t".join((mirid, gene, dname, pid));
	
t = sum(len(x) for x in r.values())	
sys.stderr.write("%d\n" % t) 

#>>>>>>>> control

for i in range(100):
	control = defaultdict();
	mirs = [x for x in real.keys()];
	for mirid, gset in real.iteritems():
		index = random.randint(0,len(mirs) - 1);
		m = mirs.pop(index);
		control[m] = gset
	r = execute(control, known);	
	tc = 	sum(len(x) for x in r.values())
	if(tc>=t):
		sys.stderr.write("%d\n" % tc) 
	