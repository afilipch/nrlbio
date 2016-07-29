#! /usr/bin/python
'''Selects miRNAs with an expression higher than provided cutoff''' 
import argparse
import sys;

from Bio import SeqIO;
from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Selects miRNAs with an expression higher than provided cutoff');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the MirBase miRNA mature sequences, fasta format");
parser.add_argument('--expression', nargs = '?', required = True, type = str, help = "Path to the expression of miRNAs, gff format")
parser.add_argument('--minexpr', nargs = '?', default = 100, type = float, help = "Min expression cutoff, in counts per million")
args = parser.parse_args();

mirnas = set()

for interval in BedTool(args.expression):
	if(float(interval.attrs['norm_expr'])>=args.minexpr):
		mirnas.add(interval.attrs['Name']);
		


for seqrecord in SeqIO.parse(args.path, 'fasta'):
	if(seqrecord.name in mirnas):
		print ">%s\n%s" % (seqrecord.id, str(seqrecord.seq.upper()).replace('U', 'T'));
