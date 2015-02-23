#! /usr/lib/python
'''assigns sequence to each gff entry'''

import argparse
import sys

from Bio import SeqIO;
from pybedtools import BedTool

parser = argparse.ArgumentParser(description='assigns sequence to each gff entry');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to gff file to be annotated");
parser.add_argument('-f', '--fasta', nargs = '?', required = True, type = str, help = "path to reference(genome) fasta file");
args = parser.parse_args();

def interval2seq(interval, reference):
	if(interval.strand == '+'):
		return str(reference[interval.chrom][interval.start:interval.stop].seq.upper())
	elif(interval.strand == '-'):
		return str(reference[interval.chrom][interval.start:interval.stop].seq.reverse_complement().upper())
	else:
		sys.stderr.write("Strand is not defined. Plus strand sequence will be returned\n")
		return str(reference[interval.chrom][interval.start:interval.stop].seq.upper())

reference = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

for interval in BedTool(args.path):
	seq = interval2seq(interval, reference);
	interval.attrs['seq'] = seq;
	sys.stdout.write(str(interval))
	