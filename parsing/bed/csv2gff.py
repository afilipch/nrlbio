#! /usr/bin/python
'''Converts csv/tsv file into gff format''' 
import argparse
import sys;
from collections import defaultdict, namedtuple
from itertools import groupby 
from operator import itemgetter, attrgetter


from nrlbio.pybedtools_extension import construct_gff_interval


parser = argparse.ArgumentParser(description='Converts csv/tsv file into gff format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the csv/tsv file. First line has to be a headers for collumns. First six collumns have to be: chromosome, start, end, name, score, strand");
parser.add_argument('-d', '--delimiter', nargs = '?', default = '\t', type = str, help = "delimeter used in csv/tsv file");
args = parser.parse_args();

with open(args.path) as f:
	header = f.next().strip().split(args.delimiter);
	attribute_names = ["_".join(x.split()) for x in header[6:]]
	for l in f:
		a = l.strip().split()
		chrom, start, end, name, score, strand = a[:6];
		attribute_values = a[6:]
		attrs = [('ID', name)] + zip(attribute_names, attribute_values)
		sys.stdout.write(str(construct_gff_interval(chrom, int(start), int(end), 'un', score=score, strand=strand, source='un', frame='.', attrs=attrs)))