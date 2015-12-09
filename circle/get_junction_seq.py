#! /usr/bin/python
'''For given circles creates fasta files with sequencing corresponding to splice junctions''' 
import sys;
import argparse
from collections import defaultdict, Counter;

from pybedtools import BedTool
from Bio import SeqIO
from Bio.Seq import reverse_complement


parser = argparse.ArgumentParser(description='For given circles creates fasta files with sequencing corresponding to splice junctions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the file with circular splice junctions, bed/gff file");
parser.add_argument('-g', '--genome', nargs = '?', required = True, type = str, help = "path to the reference(genome), fasta format")
parser.add_argument('--length', nargs = '?', default = 60, type = int, help = "length of splice junction")
args = parser.parse_args();

def get_junction(circle, seq, piece_length):
	js = "".join( [str(x) for x in  (seq[circle.end-piece_length:circle.end], seq[circle.start:circle.start+piece_length])] )
	if(circle.strand == '+'):
		return js.upper()
	else:
		return reverse_complement(js.upper())


genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))

for circle in BedTool(args.path):
	piece_length = min(len(circle), args.length)/2;
	print ">%s" % circle.name
	print get_junction(circle, genome[circle.chrom].seq, piece_length)