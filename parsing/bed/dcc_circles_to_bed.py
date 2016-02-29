#! /usr/bin/python
'''Converts DCC suplementary data tsv file into decent bed file''' 
import argparse
import sys;



parser = argparse.ArgumentParser(description='Converts DCC suplementary data tsv file into decent bed file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to DCC suplementary data tsv file");
args = parser.parse_args();


names = set();
skipped = 0;

with open(args.path) as f:
	f.next()
	for l in f:
		a = l.strip().split("\t")
		if(a[0] not in names):
			names.add(a[0]);
			start = str(int(a[2])-1)
			print "chr%s\t%s\t%s\t%s\t%s\t." % (a[1], start, a[3], a[0], a[5]);
		else:
			skipped += 1;
			
			
sys.stderr.write("%d skipped regions with identical IDs\n" % skipped)