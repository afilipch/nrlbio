#! /usr/bin/python
'''Convert Chiflex doublebed file into Circular Pipeline format''' 
import argparse
import sys;

from pybedtools import Interval

from nrlbio.generators import generator_doublebed;

parser = argparse.ArgumentParser(description='Convert Chiflex doublebed file into Circular Pipeline format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to doublebed file");
args = parser.parse_args();



for i1, i2 in generator_doublebed(args.path):
	i = Interval(i1.chrom, min(i1.start, i2.start), max(i1.end, i2.end), name=i1.name.split("|")[0], score=i1.attrs['n_uniq'], strand=i1.strand, otherfields=None)
	sys.stdout.write(str(i))