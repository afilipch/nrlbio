import parsing;
import mir_lib;
import gconst;
import conservation_lib;

import sys;
import re;
import os;

import argparse;
import multiprocessing;
import numpy;

from Bio.Seq import reverse_complement;
from collections import defaultdict, Counter;
from itertools import  permutations;
from scipy.stats import ks_2samp
from scipy.stats import hypergeom


parser = argparse.ArgumentParser(description='explore conservation of seed containing targets');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to custom maf file (interactions.bed.maf file)");
parser.add_argument('--interactions', nargs = '?', type = str, required = True, help = "path to initial interactions_inutr.bed file");
parser.add_argument('--control', nargs = '+', type = str, required = True, help = "path to control maf files");
parser.add_argument('-s', '--system', nargs = '?', required = True, type = str, help = "genome system")
parser.add_argument('-cm', '--count_mismatch', nargs = '?', default = False, const = True, help = "count conservation of mismatched seeds")
#parser.add_argument('-m', '--matches', nargs = '?', default = 0, type = int, help = "site is considered to be conserved if it is conserved in at least \'matches\' species")
args = parser.parse_args();


if(args.system == "hg19"):
	primates = ['rheMac2', 'ponAbe2', 'calJac1', 'gorGor1', 'panTro2', 'otoGar1', 'papHam1'];
	distant = ['xenTro2', 'tetNig2', 'fr2', 'gasAcu1', 'oryLat2', 'danRer6', 'petMar1']
	species = ["hg19", "anoCar1", "bosTau4", "calJac1", "canFam2", "cavPor3", "choHof1", "danRer6", "dasNov2", "dipOrd1", "echTel1","equCab2","eriEur1","felCat3","fr2","galGal3","gasAcu1","gorGor1","loxAfr3","macEug1","micMur1","mm9","monDom5","myoLuc1","ochPri2","ornAna1","oryCun2","oryLat2","otoGar1","panTro2","papHam1","petMar1","ponAbe2","proCap1","pteVam1","rheMac2","rn4","sorAra1","speTri1","taeGut1","tarSyr1","tetNig2","tupBel1","turTru1","vicPac1","xenTro2"];
	aligned_species = set(species) - set(distant) - set(primates);
elif(args.system == "ce6"):	
	species = ["ce6", "caeJap1", "caePb2", "caeRem3", "cb3", "priPac1"];
	distant = ["caeJap1", "priPac1"]

	aligned_species = set(species) - set(distant)
	
#print list(aligned_species);

#def get_mismatch2conservation(seeds):
	#r = defaultdict(list)	
	#for k, seed in seeds.iteritems():
		#for i, l in enumerate(seed.mismatched):
			#for j, mm in enumerate(l):
				#for m in mm.mismatches:
					#r[m.fr, m.to, m.position] += seed.ass_mismatched[i][j];	
	##print len(r)				
	#return r;	
	
def output_resolution_mismatch(seeds, path):
	f = open(path, 'w')
	for k, seed in seeds.iteritems():
		for i, l in enumerate(seed.mismatched):
			for j, mm in enumerate(l):
				for m in mm.mismatches:
					if(seed.ass_mismatched[i][j]):
						f.write("%s\t%d\t%s\t%s\t%s\n" % (seed.match, m.position, m.fr, m.to, ",".join([str(x) for x in seed.ass_mismatched[i][j]])));	
	f.close()
	
def output_resolution_perfect(seeds, path):
	f = open(path, 'w')
	for k, seed in seeds.iteritems():
		if(seed.ass_match):
			f.write("%s\t%d\t%s\t%s\t%s\n" % (seed.match, 0, "N", "N", ",".join([str(x) for x in seed.ass_match])));	
	f.close()	
	
	
def process_signal_perfect(path, seeds, iid2seed, aligned_species=None):
	for k, seed in seeds.iteritems():
		seed.associated_values(list);
		seed.generate_mismatches(upto=1);
	
	maf_dict = conservation_lib.read_maf(args.path[0], aligned_species=aligned_species);
	
	for iid, seqs in maf_dict.iteritems():
		k = iid2seed[iid];
		seed = seeds[k];	
		conservation_lib.conservation_match(seed, seqs, args.system);

	output_resolution_perfect(seeds, os.path.basename(args.path[0]).split(".")[0] + ".perf.cons.tsv");
	
def process_control_perfect(path, seeds, aligned_species=None):	
	for k, seed in seeds.iteritems():
		seed.associated_values(list);
		seed.generate_mismatches(upto=1);
		
	maf_dict_control = conservation_lib.read_maf_control(path, aligned_species = aligned_species)

	for iid, seqs in maf_dict_control.iteritems():
		for seed in seeds.values():	
			conservation_lib.conservation_match(seed, seqs, args.system);

	output_resolution_perfect(seeds, os.path.basename(path.split(".")[0] + ".perf.cons.tsv"));	
	
	
def process_signal_mismatches(path, seeds, iid2seed, aligned_species=None):
	for k, seed in seeds.iteritems():
		seed.generate_mismatches(upto=1);
		seed.associated_values(list);
	
	maf_dict = conservation_lib.read_maf(args.path[0], aligned_species=aligned_species);
	
	for iid, seqs in maf_dict.iteritems():
		k = iid2seed[iid];
		seed = seeds[k];	
		conservation_lib.conservation_mismatch(seed, seqs, args.system);

	output_resolution_mismatch(seeds, os.path.basename(args.path[0]).split(".")[0] + ".mm.cons.tsv");

def process_control_mismatches(path, seeds, aligned_species=None):	
	for k, seed in seeds.iteritems():
		seed.associated_values(list);
		
	maf_dict_control = conservation_lib.read_maf_control(path, aligned_species = aligned_species)

	for iid, seqs in maf_dict_control.iteritems():
		for seed in seeds.values():	
			conservation_lib.conservation_mismatch(seed, seqs, args.system);

	output_resolution_mismatch(seeds, os.path.basename(path.split(".")[0] + ".mm.cons.tsv"));

	
	
##preparing required seeds data
iid2seed = {};
seeds = {};
f = open(args.interactions)
for l in f:
	a = l.strip().split("\t");
	k = reverse_complement(a[7][1:7]);
	iid2seed[a[3]] = k;
	if(k not in seeds):
		seeds[k] = mir_lib.Seed(a[7], a[6], start=1, end = 7);
f.close();	
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

if(args.count_mismatch):
	process_signal_mismatches(args.path[0], seeds, iid2seed, aligned_species=aligned_species)
	for path in args.control:
		process_control_mismatches(path, seeds, aligned_species=aligned_species)
else:
	process_signal_perfect(args.path[0], seeds, iid2seed, aligned_species=aligned_species)
	for path in args.control:
		process_control_perfect(path, seeds, aligned_species=aligned_species)	


print len(seeds)

