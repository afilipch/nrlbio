#! /usr/bin/python
'''Script splits RNA-RNA into certain subtypes. What are these subtypes is not really clear'''
import argparse
import sys;
import os;

from pybedtools import BedTool
from collections import defaultdict

from nrlbio.generators import generator_doublebed
from nrlbio.pybedtools_extension import get_distance;

parser = argparse.ArgumentParser(description='Script splits RNA-RNA into certain subtypes. What are these subtypes is not really clear');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to interactions, double bed/gff file");
parser.add_argument('--distance', nargs = '?', default = 100, type = int, help = "Maximum allowed distance between intervals on different strands to be annotated as \'strand_switch\'");
parser.add_argument('--overlap', nargs = '?', default = 1, type = int, help = "Mimimum allowed overlap between intervals to be annotated as \'overlapped\'");
args = parser.parse_args();


def annotate(i1, i2, overlap, distance):
	if(i1.chrom == i2.chrom):
		if(i1.strand == i2.strand):
			if(-1*get_distance(i1, i2) >= overlap):
				return 'overlapped';
			else:
				return "true"
		else:
			if(get_distance(i1, i2) <= distance):
				return 'strand_switch'
			else:
				return 'true'
	else:
		return "true"	



for i1, i2 in generator_doublebed(args.path):
	itype = annotate(i1, i2, args.overlap, args.distance);
	i1.attrs['itype'] = itype
	i2.attrs['itype'] = itype
	sys.stdout.write("%s%s\n" % (i1, i2));

	
	

