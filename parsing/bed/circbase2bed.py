#! /usr/bin/python
'''Converts circbase scv file into bed format'''
import sys;
import argparse

from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Converts circbase scv file into bed format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to a circbase file, csv format");
#parser.add_argument('-a', '--attributes', nargs = '+', required = False, type = str, help = "Names of gff attributes to be written as additional bed fields in an order they are provided");
args = parser.parse_args();


with open(args.path) as f:
	for l in f:
		a = l.strip().split(",")
		chrom, positions = a[0].split(":");
		start, end = [int(x) for x in positions.split("-")];
		strand = a[1];
		name = a[2];
		score = '0';
		print "%s\t%d\t%d\t%s\t%s\t%s" % (chrom, start, end, name, score, strand);
