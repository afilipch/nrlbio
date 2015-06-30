#! /usr/lib/python
'''Outputs reverse complement for STDIN sequences'''
import argparse;
from Bio.Seq import reverse_complement;

parser = argparse.ArgumentParser(description='Outputs reverse complement for STDIN sequences');
parser.add_argument('sequences', metavar = 'N', nargs = '+', type = str, help = "sequences to get reverse complement");
args = parser.parse_args();

for seq in args.sequences:
	print reverse_complement(seq);