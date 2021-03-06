#! /usr/bin/python
'''Extra''' 
import argparse
import sys;

from Bio import SeqIO;


parser = argparse.ArgumentParser(description='Creates a table which connects MirBase names to the UCSC genomes names');
#parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the MirBase miRNA mature sequences, fasta format");
parser.add_argument('--mirbase', nargs = '?', required = True, type = str, help = "Path to the MirBase miRNA mature sequences, fasta format")
parser.add_argument('--ucsc', nargs = '?', required = True, type = str, help = "Path to the UCSC maf alignment description, txt file")
args = parser.parse_args();

mirbase_names = {}

for seqrecord in SeqIO.parse(args.mirbase, 'fasta'):
	mname = seqrecord.name.split("-")[0]
	specie_name = tuple(seqrecord.description.split(' ')[2:4])
	mirbase_names[specie_name] = mname
	
#for key in mirbase_names.keys():
	#print key
	
	
#m = set()


	
with open(args.ucsc) as f:
	for l in f:
	a = [x for x in l.strip().split(" ") if x];
	for r in a:
		ga = r.split('/') 
		if(len(ga) == 2 and ga[1][-1].isdigit()):
			print ga[1]