#! /usr/bin/python
'''Selects mature miRNAs for specie of interest'''
import argparse;

from Bio import SeqIO

parser = argparse.ArgumentParser(description='Selects mature miRNAs for specie of interest');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the fasta file");
parser.add_argument('--specie', nargs = '?', type = str, required = True, help = "Specie name, like it is written in MirBase file (hsa, mmu, cel)")
args = parser.parse_args();

for seqrecord in SeqIO.parse(args.path, "fasta"):
	if(seqrecord.id.split('-')[0] == args.specie):
		print ">%s" % seqrecord.id
		print str(seqrecord.seq.upper()).replace('U', 'T')