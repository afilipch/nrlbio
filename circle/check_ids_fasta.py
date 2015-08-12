#! /usr/bin/python
'''Checks if reads which get rise to circles are actually in fasta file''' 
import argparse
import sys;

from Bio import SeqIO



parser = argparse.ArgumentParser(description='Checks if reads which get rise to circles are actually in fasta file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the reads, fasta file");
parser.add_argument('--circle2read', nargs = '?', required = True, type = str, help = "path to a dictionary file, which connects circular ids with read ids, it is based on. fasta format");
args = parser.parse_args();


circids = set();
for seqrecord in SeqIO.parse(args.circle2read, 'fasta'):
	circids.add(seqrecord.description.split(" ")[1])
total = len(circids);	

found = 0;
for seqrecord in SeqIO.parse(args.path, 'fasta'):
	if(seqrecord.name in circids):
		found+=1;
		circids = circids - set((seqrecord.name,))
		
		
		


print ("total number of circular reads:\t%d\nnumber of circular reads detected:\t%d\nratio:\t%1.5f") % (total, found, float(found)/total);
