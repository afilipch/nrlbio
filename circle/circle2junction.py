#! /usr/bin/python
'''For given circles (sequences) creates fasta files with sequencing corresponding to splice junctions''' 
import sys;
import argparse
from collections import defaultdict, Counter;


from Bio import SeqIO



parser = argparse.ArgumentParser(description='For given circles (sequences) creates fasta files with sequencing corresponding to splice junctions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the file with circular splice junctions, fasta file");
parser.add_argument('--length', nargs = '?', default = 60, type = int, help = "length of splice junction")
args = parser.parse_args();

piece = args.length/2;


for seqrecord in SeqIO.parse(args.path, 'fasta'):
	seq = str(seqrecord.seq.upper())
	apiece = min(piece, len(seq)/2)
	print ">%s|csj|%d" % (seqrecord.name, apiece*2)
	print "".join(( seq[-apiece:], seq[:apiece] ))
	