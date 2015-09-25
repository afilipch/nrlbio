#! /usr/bin/python
'''Converts single formatted linear/circular splice junctions into double bed/gff format format''' 
import sys;
import argparse

from pybedtools import BedTool, Interval

#from nrlbio.generators import generator_doublebed;
#from nrlbio.pybedtools_extension import construct_gff_interval


parser = argparse.ArgumentParser(description='Converts single formatted linear/circular splice junctions into double bed/gff format format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the chimeras, single bed/gff file");
parser.add_argument('--jtype', nargs = '?', required = True, choices=['csj', 'lsj'], type = str, help = "type of junction [csj|lsj]")
parser.add_argument('--length', nargs = '?', default=50, type = int, help = "length of new intervals produced from a single record")
args = parser.parse_args()



def linear(interval, length):
	start1 = interval.start - length;
	stop1 = interval.start;	
	start2 = interval.stop;
	stop2 = interval.stop+length;
	
	if(interval.strand == '-'):
		start1, start2 = start2, start1
		stop1, stop2 = stop2, stop1
		
	i1 = Interval(interval.chrom, start1, stop1, name="|".join((interval.name, "0")), score=interval.score, strand=interval.strand, otherfields=None)
	i2 = Interval(interval.chrom, start2, stop2, name="|".join((interval.name, "1")), score=interval.score, strand=interval.strand, otherfields=None)
	
	return i1, i2
	


def circular(interval, length):
	start1 = interval.stop - length;
	stop1 = interval.stop;
	start2 = interval.start;
	stop2 = interval.start+length;
	
	if(interval.strand == '-'):
		start1, start2 = start2, start1
		stop1, stop2 = stop2, stop1
		
	i1 = Interval(interval.chrom, start1, stop1, name="|".join((interval.name, "0")), score=interval.score, strand=interval.strand, otherfields=None)
	i2 = Interval(interval.chrom, start2, stop2, name="|".join((interval.name, "1")), score=interval.score, strand=interval.strand, otherfields=None)
	
	return i1, i2



if(args.jtype=='lsj'):
	func = linear;
else:
	func = circular;


for interval in BedTool(args.path):
	sys.stdout.write("%s%s" % func(interval, args.length));

	
