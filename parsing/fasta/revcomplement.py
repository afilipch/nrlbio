#! /usr/bin/python
'''Outputs reverse complement for given fasta file'''
import argparse;

from Bio.Seq import reverse_complement;
from Bio import SeqIO;

parser = argparse.ArgumentParser(description='Outputs reverse complement for given fasta file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to fasta file");
parser.add_argument('-ll', '--line_length', nargs = '?', default = 50, type = int, help = "Max length of line in output fasta file. Makes it more readable");
args = parser.parse_args();

def reverse_fast(seqrecord, line_length):
	line = [];
	for s in seqrecord.seq.reverse_complement():
		line.append(s);
		if(len(line) == line_length):
			yield "".join(line);
			line[:] = [];
	if(line):
		yield "".join(line);


for seqrecord in SeqIO.parse(args.path, 'fasta'):
	print ">revcompl_%s" % seqrecord.name
	for ss in reverse_fast(seqrecord, args.line_length):
		print ss;