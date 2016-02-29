#! usr/bin/python
'''Checks for identical sequences in provided fasta file. Can be used only for small fasta file, will be improved further'''
import os;
import sys;
import argparse;
from collections import defaultdict, Counter
#from itertools import product

from Bio import SeqIO




parser = argparse.ArgumentParser(description='Checks for identical sequences in provided fasta file. Can be used only for small fasta file, will be improved further');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to fasta file"); 
args = parser.parse_args();

sequences = defaultdict(int);


for seqrecord in SeqIO.parse(args.path, 'fasta'):
	sequences[str(seqrecord.seq.upper())] += 1;
	
	
print Counter(sequences.values());
	
	
	
