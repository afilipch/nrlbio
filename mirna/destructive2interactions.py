#! /usr/lib/python
'''Converts potential candidates to be destructive sites into interactions double gff for the sake of downstream compatibility'''

import argparse
import os
import sys
from collections import defaultdict

from pybedtools import BedTool

from nrlbio.pybedtools_extension import construct_gff_interval



parser = argparse.ArgumentParser(description='Tool for destructive miRNA sites prediction');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the potential destructive sites, bed format");
args = parser.parse_args();

for c, interval in enumerate(BedTool(args.path)):
	for mirid in interval.name.split(','):
		i1 = construct_gff_interval(mirid, 0, 30, 'ds', score=interval.score, strand='+', source='un', frame='.', attrs=[('ID', "cid%d|0" % (c+1)), ('n_uniq', 1)])
		i2 = construct_gff_interval(interval.chrom, interval.start, interval.stop, 'ds', score=interval.score, strand=interval.strand, source='un', frame='.', attrs=[('ID', "cid%d|1" % (c+1)), ('n_uniq', 1)])
		sys.stdout.write(str(i1));
		sys.stdout.write(str(i2));