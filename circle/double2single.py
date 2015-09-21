#! /usr/bin/python
'''Converts doublebed formatted circles into single bed lines''' 
import sys;
import argparse

from nrlbio.generators import generator_doublebed;
from nrlbio.pybedtools_extension import construct_gff_interval


parser = argparse.ArgumentParser(description='Converts doublebed formatted circles into single bed lines');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the chimeras, double bed/gff file");
args = parser.parse_args()


def double2single(i1, i2):
	start = min(i1.start, i2.start)
	stop = max(i1.stop, i2.stop)
	score = str(int(i1.score) + int(i2.score));
	circ_attrs = [('gap', i1.attrs['gap']), ('score1', i1.score), ('score2', i2.score), ('ID', i1.name.split("|")[0])]
	
	return construct_gff_interval(chrom=i1.chrom, start=start, stop=stop, feature='circle', score=score, strand=i1.strand, source='.', frame='.', attrs=circ_attrs)


for i1, i2 in generator_doublebed(args.path):
	sys.stdout.write(str(double2single(i1, i2)));