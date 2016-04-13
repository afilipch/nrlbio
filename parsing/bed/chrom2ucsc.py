#! /usr/bin/python
'''Converts chromosome name from ENSEMBL name convention to the UCSC one'''
import sys;
import argparse

from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Converts chromosome name from ENSEMBL name convention to the UCSC one');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to a file in gff format");
parser.add_argument('--mname', nargs = '?', default = 'MT', type = str, help = "Name of the mitochondrial DNA  in ENSEMBL");
args = parser.parse_args();

mname = args.mname


for interval in BedTool(args.path):
	if(interval.chrom == mname):
		interval.chrom = 'chrM'
	else:
		interval.chrom =  "chr%s" % interval.chrom
	sys.stdout.write(str(interval))