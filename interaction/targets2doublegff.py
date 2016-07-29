#! /usr/bin/python
'''Converts single gff to double gff for miRNA:target interactions''' 
import argparse
import sys;
from collections import defaultdict

from pybedtools import BedTool

from nrlbio.pybedtools_extension import construct_gff_interval

parser = argparse.ArgumentParser(description='Converts single gff to double gff for miRNA:target interactions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the interactions, gff format");
args = parser.parse_args();


for interval in BedTool(args.path):
	attrs1 = [('n_uniq', interval.attrs['n_uniq']), ('chscore', interval.attrs['chscore']), ('gap', interval.attrs['gap']), ('ID', interval.attrs['ID']+"|0")]
	attrs2 = [('n_uniq', interval.attrs['n_uniq']), ('chscore', interval.attrs['chscore']), ('gap', interval.attrs['gap']), ('ID', interval.attrs['ID']+"|1")]
	i1 = construct_gff_interval(interval.attrs['mirid'], 0, 22, 'ch', score=interval.score, strand='+', source='ch', frame='.', attrs=attrs1)
	i2 = construct_gff_interval(interval.chrom, interval.start, interval.end, 'ch', score=interval.score, strand=interval.strand, source='ch', frame='.', attrs=attrs2)
	sys.stdout.write("%s%s" % (i1, i2))
