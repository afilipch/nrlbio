#! /usr/bin/python
'''Selects for targets with a destructive potential''' 
import argparse
import sys;

from Bio import SeqIO;


parser = argparse.ArgumentParser(description='Selects miRNA sequences from MirBase mature.fa file for a particular specie');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the MirBase miRNA mature sequences, fasta format");
parser.add_argument('--specie', nargs = '?', required = True, type = str, help = "Name of the specie to select (MirBase coding: cel, hsa, mmu and so on)")
args = parser.parse_args();


for seqrecord in SeqIO.parse(args.path, 'fasta'):
	if(seqrecord.name.split("-")[0] == args.specie):
		print ">%s" % seqrecord.id
		print str(seqrecord.seq.upper()).replace('U', 'T');