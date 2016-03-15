#! /usr/lib/python
'''Converts old formatted interactions into ChiFlex format'''

import argparse
import sys

from pybedtools import BedTool

from nrlbio.pybedtools_extension import construct_gff_interval


parser = argparse.ArgumentParser(description='Converts old formatted interactions into ChiFlex format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the olf-fashion interactions, bed file");
args = parser.parse_args();

for interval in BedTool(args.path):
	for mirid in interval[6].split(","):
		sys.stdout.write(str(construct_gff_interval(mirid, 0, 22, 'newint', score=interval.score, strand='+', source='un', frame='.', attrs=[('ID', "%s|0" % interval.name), ('seq', interval[7]), ('n_uniq', interval[11])])));
		sys.stdout.write(str(construct_gff_interval(interval.chrom, interval.start, interval.stop, 'newint', score=interval.score, strand=interval.strand, source='un', frame='.', attrs=[('ID', "%s|1" % interval.name), ('seq', interval[8]), ('n_uniq', interval[11])])));