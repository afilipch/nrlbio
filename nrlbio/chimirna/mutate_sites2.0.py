'''fast, but suboptimal version of in silico mutation'''
#! /usr/bin/python
import sys;
import parsing;
import re;
import os;
import time;
import csv;
import copy;
import argparse;
import random
import multiprocessing;
import itertools
import commands;
from collections import *;
from Bio import SeqIO;
from generators import *;
from Bio import pairwise2;
from Bio.Seq import reverse_complement


parser = argparse.ArgumentParser(description='tries to mutate binding sites for microRNAs to increase binding energy and not create seed matches for other top expressed miRNAs STDOUT');
parser.add_argument('-m', '--manual', nargs = '?', type = argparse.FileType('r'), help = "path to manually curated file in tsv format: miRNA id, miRNA seq, binding site id, binding site seq, mutable letters in lowercase");
parser.add_argument('-t', '--top_expressed', nargs = '?', type = argparse.FileType('r'), required = True, help = "path to top expressed miRNA");
parser.add_argument('--mir', nargs = '?', type = argparse.FileType('r'), required = True, help = "path to fasta file of microRNAs");
#parser.add_argument('--max_mutations', nargs = '?', type = int, default = 1, help = "maximum number of mutations allowed");
#parser.add_argument('--show', nargs = '?', type = int, default = 5, help = "number of mutants to output");
parser.add_argument('--modes', nargs = '*', default = ['nocreate', 'nodisrupt'], choices=['nocreate', 'nodisrupt'], help = "can be nocreate and/or nodisrupt seed matches");
args = parser.parse_args();

if("nodisrupt" in args.modes and "nocreate" in args.modes):
	mode = r'site.seeds == tuple(s)'
elif("nocreate" in args.modes):
	mode = r'set(site.seeds) >= set(tuple(s))'
elif("nodisrupt" in args.modes):
	mode = r'set(site.seeds) <= set(tuple(s))'
else:
	mode = 'True';

#print mode	
	
class Site(object):
	"""positions - mutable positions in site, mutated - list of new mutated sequences with the highest energy; in form sequence energy"""
	def __init__(self, mirna_id, mirna_seq, site_id, site_seq, positions, start, end):
		self.mirna_id, self.mirna_seq, self.site_id, self.site_seq, self.positions = mirna_id, mirna_seq, site_id, site_seq, positions
		self.positions.sort();
		self.seed = reverse_complement(self.mirna_seq[start:end]);
		self.current = site_seq; 
		self.mutated = [];
	def get_seeds(self, seeds):	
		s = [];
		for seed in seeds:
			if(seed in self.site_seq):
				s.append(seed)				
		#s.append(self.seed)		
		self.seeds = tuple(s);

	
def mutate(mirna, target):
	#mirna, target = arr;
	call = "RNAhybrid -s 3utr_worm \"" + target + "\"  \"" + mirna + "\" -b 1 -c ";
	status, output = commands.getstatusoutput(call)
	if(not output):
		return -100; 
	energy = float(output.split(":")[4]);
	return target, energy
	
def test_target(site, target, seeds, mode):
	s = [];
	for seed in seeds:
		if(seed in target):
			s.append(seed)
	if(eval(mode)):# and site.seed not in target):
		return True
	else:
		return False
		
def mutate_piece(site, seeds, mode, length = 4):
	targets = set();
	lseq = list(site.current)
	mnum = min(length, len(site.positions))
	mutable_letters = list(itertools.product("ACTG", repeat = mnum))
	for ml in mutable_letters:
		target = copy.copy(lseq);
		for i in range(mnum):
			if(False and ml[i] == site.site_seq[site.positions[i]]):
				break;
			target[site.positions[i]] = ml[i];
		else:	
			target = "".join(target);
			if(test_target(site, target, seeds, mode)):
				targets.add(target);
	
	r = [];
	for target in targets:
		r.append(mutate(site.mirna_seq, target));
		
	#>>>> mutate all of them ??	
	#print site.current
	if(mnum <= len(site.positions)):
		site.mutated = list(sorted(r, key = lambda x: x[1], reverse = True))[:10]
		
	site.current = max(r, key = lambda x: x[1])[0];
	site.positions = site.positions[mnum:];
	
def run(site, seeds, mode, length = 4):
	while site.positions:
		mutate_piece(site, seeds, mode, length)
		#represent(site)
			
	
		

	


	
def represent(site):
	#stub, energy = mutate((site.mirna_seq, site.site_seq))
	print "-"*100
	print "miRNA %s %s\tintact binding site %s %s" % (site.mirna_id, site.mirna_seq, site.site_id, site.site_seq)
	call = "RNAhybrid -s 3utr_worm \"" + site.site_seq + "\"  \"" + site.mirna_seq + "\" -b 1 ";
	status, output = commands.getstatusoutput(call)
	print "\n".join(output.split("\n")[4:13]);
	print "\nMutants:\n";
	
	mutant = ""
	for j in range(len(site.site_seq)):
		char = site.site_seq[j]
		if(char == site.current[j]):
			mutant += site.current[j];
		else:
			mutant += site.current[j].lower()
	
	print "miRNA %s %s\tmutated binding site number: %s" % (site.mirna_id, site.mirna_seq, mutant)
	call = "RNAhybrid -s 3utr_worm \"" + site.current + "\"  \"" + site.mirna_seq + "\" -b 1";
	status, output = commands.getstatusoutput(call)
	print "\n".join(output.split("\n")[4:13]);
	print "\n";
	
def represent_variants(site):
	#stub, energy = mutate((site.mirna_seq, site.site_seq))
	print "-"*100
	print "miRNA %s %s\tintact binding site %s %s" % (site.mirna_id, site.mirna_seq, site.site_id, site.site_seq)
	call = "RNAhybrid -s 3utr_worm \"" + site.site_seq + "\"  \"" + site.mirna_seq + "\" -b 1 ";
	status, output = commands.getstatusoutput(call)
	print "\n".join(output.split("\n")[4:13]);
	print "\nMutants:\n";
	
	for i in range(len(site.mutated)):
		m = site.mutated[i][0]
		mutant = ""
		for j in range(len(site.site_seq)):
			char = site.site_seq[j]
			if(char == m[j]):
				mutant += m[j];
			else:
				mutant += m[j].lower()
		
		print "miRNA %s %s\tmutated binding site number:%d %s" % (site.mirna_id, site.mirna_seq, i+1, mutant)
		call = "RNAhybrid -s 3utr_worm \"" + m + "\"  \"" + site.mirna_seq + "\" -b 1";
		status, output = commands.getstatusoutput(call)
		print "\n".join(output.split("\n")[4:13]);
		print "\n";	
	
		
		
	

def get_top_mirna():
	top = set();
	ans = {};
	for line in args.top_expressed:
		top.add(line.strip());
	temp = SeqIO.to_dict(SeqIO.parse(args.mir,"fasta"))
	for key, value in temp.iteritems():
		if(key in top):
			ans[key] = str(value.seq).upper()
	args.mir.close()
	args.top_expressed.close()
	return ans;
	
def get_manual(seeds, start, end):
	sites = [];
	for line in args.manual:
		arr = line.strip().split("\t")
		positions = [];
		for i in range(len(arr[3])):
			char = arr[3][i];
			if(char.islower()):
				positions.append(i)
		site = 	Site(arr[0], arr[1], arr[2], arr[3].upper(), positions, start, end)
		site.get_seeds(seeds)
		sites.append(site)		
	args.manual.close()	
	return sites;

	
	

top = get_top_mirna()
seeds = set()
for k,v in top.iteritems():
	seeds.add(reverse_complement(v[1:7]));

if(args.manual):
	sites = get_manual(seeds, 1, 6)

#print sites[0].seeds	
	
for site in sites[13:14]:
	#print site.positions
	run(site, seeds, mode, 4)
	represent_variants(site);
	#mutate_piece(site, seeds, mode)
	#represent(site);	
	
#site1 = sites[8];	
#site2 = sites[9];
#site2.mutated = site1.mutated;
#represent_variants(site1)
#represent_variants(site2);
	
#set1 = set([x[0] for x in sites[8].mutated]);
#set2 = set([x[0] for x in sites[9].mutated]);
#for el in set1 & set2:
	#sys.stderr.write("%s\n" % el)


	
