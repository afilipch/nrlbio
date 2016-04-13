#! /usr/bin/python
'''Converts chromosome name from UCSC name convention to the ENSEMBL one'''
import sys;
import argparse

from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Converts chromosome name from UCSC name convention to the ENSEMBL one');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to a file in gff format");
parser.add_argument('--mname', nargs = '?', default = 'MT', type = str, help = "Name of the mitochondrial DNA in ENSEMBL");
args = parser.parse_args();

mname = args.mname


for interval in BedTool(args.path):
	if(interval.chrom == 'chrM'):
		interval.chrom = mname
	else:
		interval.chrom =  interval.chrom[3:]
	sys.stdout.write(str(interval))