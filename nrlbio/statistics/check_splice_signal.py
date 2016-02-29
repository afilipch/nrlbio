#! usr/bin/python
'''Script outputs statistics on edge nucleotides of provided sequences'''
import os;
import sys;
import argparse;
from collections import defaultdict
from itertools import product

from Bio import SeqIO




parser = argparse.ArgumentParser(description='Script outputs statistics on edge nucleotides of provided sequences');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to fasta file"); 
args = parser.parse_args();

start = defaultdict(int);
end = defaultdict(int);
total = 0.0;

for seqrecord in SeqIO.parse(args.path, 'fasta'):
	start[str(seqrecord.seq[:2]).upper()] += 1;
	end[str(seqrecord.seq[-2:]).upper()] += 1;
	total += 1;
	
	
	
	
	
print "dinucleotide\ttotal start\tfraction start\ttotal end\tfraction end"	
for p in product('ACTG', repeat=2):
	dn = "".join(p);
	print "%s\t%d\t%1.3f\t%d\t%1.3f" % (dn, start[dn], start[dn]/total, end[dn], end[dn]/total)
	
	

