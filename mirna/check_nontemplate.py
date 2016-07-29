#! /usr/bin/python
'''Checks if the chimeric read carry  information of nontemplate addition of nucleotides to miRNAs''' 
import argparse
import sys;
from collections import defaultdict

import pysam;
from Bio import SeqIO



parser = argparse.ArgumentParser(description='Checks if the chimeric read carry  information of nontemplate addition of nucleotides to miRNAs');
parser.add_argument('path', metavar = 'N', nargs = 2, type = str, help = "path to the mapping hits for the first and second run");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "path to the miRNA sequences, fasta file");
args = parser.parse_args();

mirid2len = {};
for seqrecord in SeqIO.parse(args.mir, 'fasta'):
	mirid2len[seqrecord.name] = len(seqrecord)


first = pysam.Samfile(args.path[0]);
second = pysam.Samfile(args.path[1]);

cutdict = defaultdict(float)
addnucl = defaultdict(float)
addlen = defaultdict(float)

for s1, s2 in zip(first.fetch(until_eof=True), second.fetch(until_eof=True)):
	rname = first.getrname(s1.tid)
	if(s1.query_name != s2.query_name):
		sys.stderr.write('WARNING: Hits in first and second files are not ordered. Please, order themn according to the query names\n');
	else:
		cut = mirid2len[rname] - s1.reference_end
		cutdict[cut] += 1
		
		nucl = s2.query_sequence[:s2.query_alignment_start]
		addnucl[nucl] += 1
		addlen[len(nucl)] += 1
		
		
print addnucl;
print cutdict;
		
first.close();
second.close();