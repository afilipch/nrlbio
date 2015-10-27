#! /usr/bin/python
'''Converts gff file into csv/tsv format''' 
import argparse
import sys;

from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Converts gff file into csv/tsv format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the gff file");
parser.add_argument('-d', '--delimiter', nargs = '?', default = '\t', type = str, help = "delimeter to use in csv/tsv file");
args = parser.parse_args();

delim = args.delimiter;

bed = BedTool(args.path)
akeys = bed[0].attrs.keys()
header = ['chrom', 'start', 'end', 'name', 'score', 'strand'] + akeys;
print delim.join(header)
for interval in bed:
	print delim.join([interval.chrom, str(interval.start), str(interval.end), interval.name, interval.score, interval.strand] + [interval.attrs.get(x) for x in akeys]);
	#print "%s\t%d\t%d\t%s\t%s\t%s\t%s" % (interval.chrom, interval.start, interval.end, interval.name, interval.score, interval.strand, "\t".join([interval.attrs.get(x) for x in akeys]));
