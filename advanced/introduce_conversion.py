#! /usr/lib/python
'''introduces conversion into fastq reads. Generates variants of the same read with conversions introduced at every possible spot'''
import argparse;

from nrlbio.generators import generator_fastq;

parser = argparse.ArgumentParser(description='introduces conversion into fastq reads. Generates variants of the same read with conversions introduced at every possible spot');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to fasta file");
parser.add_argument('--From', nargs = '?', required = True, type = str, help = "conversion from");
parser.add_argument('--To',   nargs = '?', required = True, type = str, help = "conversion to");
args = parser.parse_args();

def get_variants(seq, From, To):
	variants = [(seq, -1)];
	for p, n in enumerate(seq):
		if(n == From):
			variants.append(("".join([seq[:p], To, seq[p+1:]]), p))
	return variants		


for fastq in generator_fastq(args.path):
	for v in  get_variants(fastq.seq, args.From, args.To):
		print "\n".join(["%s_%s:%s:%d" % (fastq.id, args.From, args.To, v[1]), v[0], fastq.sign, fastq.qual])