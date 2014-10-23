#! /usr/lib/python
'''Converts expected nucleotides in the end of matches into 1, otherwise into 0'''
import argparse;
import sys;

parser = argparse.ArgumentParser(description='Converts expected nucleotides in the end of matches into 1, otherwise into 0');
parser.add_argument('-b', '--bias', nargs = '+', required = True, choices = ['G', 'T', 'C', 'A'], type = str, help = "expected nucleotides");
parser.add_argument('-i', '--indices', nargs = '+', required = True, type = int, help = "indices of the nucleotide fields in chimera file");
args = parser.parse_args();

def reassign(n):
	if(n in args.bias):
		return '1';
	else:
		return '0'
	


for l in sys.stdin:
	a = l.strip().split("\t");
	for i in args.indices:
		a[i] = reassign(a[i]);
	print "\t".join(a);	

