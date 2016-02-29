#! /usr/bin/python
'''Converts headers of artificialy generated reads into doublebed format'''
import argparse;
import sys

from Bio import SeqIO

from nrlbio.pybedtools_extension import construct_gff_interval


parser = argparse.ArgumentParser(description='Converts headers of artificialy generated reads into doublebed format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the fasta file");
args = parser.parse_args();


def regions2intervals(part):
	for r in regions.split("&"):
		chrom, strand, start, end = r.split(":")
		start, end = int(start), int(end)
		yield chrom, strand, start, end
		

	
for seqrecord in SeqIO.parse(args.path, "fasta"):
	number, read_type, regions, genomic_type, left_right, conversion = seqrecord.id.split("|")
	left, right, gap  = left_right.split(":")
	for (chrom, strand, start, end), conv in zip(regions2intervals(regions), conversion.split(':')):
		sys.stdout.write(str(construct_gff_interval(chrom, start, end, read_type, score='0', strand=strand, source='ag', frame='.', attrs=[('conv', conv), ('ID', number), ('left', left), ('right', right), ('gap', gap), ('gtype', genomic_type), ('rtype', read_type)])))
	