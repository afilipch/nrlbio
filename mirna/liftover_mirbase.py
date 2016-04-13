#! /usr/bin/python
'''Creates a table which connects ids from one miRBase version with another one based on sequence''' 
import argparse
import sys;

from Bio import SeqIO;


parser = argparse.ArgumentParser(description='Creates a table which connects ids from one miRBase version with another one based on sequence');
parser.add_argument('path', metavar = 'N', nargs = 2, type = str, help = "Path to the MirBase miRNA mature sequences. Id's form the first file will be mapped to the second");
args = parser.parse_args();



forward = {}
reverse = {}

for seqrecord in SeqIO.parse(args.path[1], 'fasta'):
	mirid = seqrecord.id
	mseq = str(seqrecord.seq.upper()).replace('U', 'T');
	forward[mirid] = mseq
	reverse[mseq] = mirid
	
	
for seqrecord in SeqIO.parse(args.path[0], 'fasta'):
	mirid = seqrecord.id
	mseq = str(seqrecord.seq.upper()).replace('U', 'T');
	if(mirid in forward):
		print "%s\t%s" % (mirid, mirid)
	elif(mseq in reverse):
		print "%s\t%s" % (mirid, reverse[mseq])
	else:
		print "%s\t%s" % (mirid, mirid)