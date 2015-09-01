import parsing;
import mir_lib;
import gconst;
import sys;
import re;
import os;
#import time;
#import csv;
#import copy;
import argparse;
import multiprocessing;
import numpy;
#import random;
from Bio.Seq import reverse_complement;
#from Bio import pairwise2;
from collections import defaultdict, Counter;
from itertools import  permutations;
from scipy.stats import ks_2samp
from scipy.stats import hypergeom


parser = argparse.ArgumentParser(description='explore conservation of seed containing targets');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to custom maf file (interactions.bed.maf file)");
parser.add_argument('--interactions', nargs = '?', type = str, required = True, help = "path to initial interactions_inutr.bed file");
parser.add_argument('--control', nargs = '+', type = str, required = True, help = "path to control maf files");
parser.add_argument('-s', '--system', nargs = '?', required = True, type = str, help = "genome system")
parser.add_argument('-m', '--matches', nargs = '?', default = 0, type = int, help = "site is considered to be conserved if it is conserved in at least \'matches\' species")
args = parser.parse_args();

aligned_species = gconst.conservation[args.system]; 
my_specie = args.system;
if(args.matches < 1):
	matches = len(aligned_species) - 1;
else:
	matches = args.matches;

def get_control(path):
	maf_control = defaultdict(dict);
	f = open(path)
	for l in f:
		l = l.strip()  
		if(l and l[0] == ">"):
			a = l.split(":");
			if(len(a) > 2):
				specie_raw, stub = a[:2]
				gene = ":".join((a[2:]))
				specie = specie_raw[1:].split(".")[0];
			else:
				specie = l[1:];      
		elif(l and specie in aligned_species):
			maf_control[gene][specie] = l.upper();
		else:
			pass
	f.close();
	return maf_control;



def conservation(seed, dictionary, pos, matches):
	m = 0;
	for key, seq in dictionary.iteritems():
		if(seed == seq[pos:pos+6]):
			m += 1;
	if(m <= matches):
		return 0;
	else:
		return 1;
	
  
def count_conservation(maf_dict, iid2seed, matches):
	counts = defaultdict(int);# list of conservation scores
	total = defaultdict(int);
	for iid, d in maf_dict.iteritems():
		seed = iid2seed[iid]
		pos = d[my_specie].find(seed)
		if(pos > -1):
			counts[seed] += conservation(seed, d, pos, matches)
			total[seed] += 1;
	return counts, total 
  
  
def count_control(maf_control, seeds, matches):
	counts = defaultdict(int);# list of conservation scores
	total = defaultdict(int);
	for g,d in maf_control.iteritems():

		for seed in seeds:      
			positons = parsing.findsubstring(d[my_specie], seed)
			for pos in positons:
				counts[seed] += conservation(seed, d, pos, matches)
				total[seed] += 1;
	return counts, total     
  
#def represent(counts_true, counts_control, total_true, total_control, seed_dict, path):
  #order = sorted(total_true.keys(), key = lambda x: total_true[x], reverse = True)  
  #simple = open(path, 'w');
  #for seed in order:
    #name = sorted(seed_dict[seed], key = lambda x: int(re.search(r"\d+", x).group()))[0]
    ##print name
    ##name = re.search("\w+-\d+", name).group()
    #simple.write("%s\t%d\t%d\t%1.3f\t%d\t%d\t%1.3f\n" % (name, counts_true[seed], total_true[seed], float(counts_true[seed])/total_true[seed], counts_control[seed], total_control[seed], float(counts_control[seed])/total_control[seed]))
  #simple.close()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> get seeds for miRNA which are present in interactions 
  
     
#>>>>>>>>>>>>>>>we have to know what miRNA(seed) is part of interaction & get seeds for miRNA which are present in interactions
iid2seed = {};
mirids = set()
f = open(args.interactions)
for l in f:
	a = l.strip().split("\t");
	iid2seed[a[3]] = reverse_complement(a[7][1:7]);
	mirids.update(a[6].split(","));
f.close();	  

seeds = set()
mirdict = parsing.fasta2dict(gconst.system2mir[args.system])
for k,v in mir_lib.mir2seed(mirdict).iteritems():
	if(k in mirids):
		seeds.add(v)
		
seed2mirs = mir_lib.seed2mirs(mirdict)	


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  

#>>>>>>>>>>>>>>>>>>>>>>> read maf file of signal
signal = defaultdict(dict)# keys are names of targets (hsa-let-7c_12). Values are dicts: keys are specie names, values are sequences
f = open(args.path[0])
for l in f:
	l = l.strip();  
	if(l and l[0] == ">"):
		a = l.split(":");
		if(len(a) == 3):
			specie_raw, stub, iid = a
			specie = specie_raw[1:].split(".")[0];
		else:
			specie = l[1:];
	elif(l and specie in aligned_species):
		signal[iid][specie] = l.upper();
	else:
		pass
f.close();
sys.stderr.write("%d blocks were read from\t%s\n\n" % (len(signal), args.path[0]))
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>read maf files of control
control = []# is a list of:   maf_control = defaultdict(dict)#keys are 3utr gene names. Values are dicts: keys are specie names, values are sequences
for path in args.control:
	control.append(get_control(path))
	sys.stderr.write("%d blocks were read from\t%s\n" % (len(control[-1]), path))

counts, total = count_conservation(signal, iid2seed, matches);	

r_control = []
for maf_control in control:
	cc,tc = count_control(maf_control, seeds, matches)
	r_control.append((cc,tc))
	

for k,v in Counter(total).most_common(10):
	name = mir_lib.shortname(seed2mirs[k])
	a = [name, counts[k], v, counts[k]/float(v)]
	for cc,tc in r_control:
		rawp = (1 - hypergeom.cdf(counts[k], tc[k] ,cc[k] , v)) 
		p_value = "%.2e" % rawp;
		a += [cc[k], tc[k], cc[k]/(tc[k] + 0.1), p_value]
	print "\t".join([str(x) for x in a])

#micros = parsing.getDictFromFa(args.micro);
#micros = parsing.to_upper_case(micros);

#for key, value in micros.iteritems():
  #seed_dict[reverse_complement(value[1:7])].append(key);


#counts_true, total_true = count_conservation(maf_dict, micros)
#seeds = set([reverse_complement(x[1:7]) for x in micros.values()])
#counts_control, total_control  = count_control(maf_control, seeds)


#represent(counts_true, counts_control, total_true, total_control, seed_dict, "micro_conservation.tsv")


