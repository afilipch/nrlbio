#! /usr/lib/python
'''Converts refseq gene prediction track into exons bed file'''
import argparse
import sys;

from nrlbio.genome_system import refseq2interval



parser = argparse.ArgumentParser(description='Converts refseq gene prediction track into exons bed file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to bed file");
args = parser.parse_args();


f = open(args.path);
f.readline();
for i in refseq2interval(f):
	sys.stdout.write(str(i))