#! /usr/bin/python
'''Removes interactions coming from the same gene'''
import sys;
import argparse;
import os

from nrlbio.generators import generator_doublebed




parser = argparse.ArgumentParser(description='Removes interactions coming from the same gene');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the annotated interactions, double gff format");
args = parser.parse_args();




passed, removed = 0, 0
for i1, i2 in generator_doublebed(args.path):
	g1 = set([x for x in i1.attrs['gene_symbols'].split(':') if x])
	g2 = set([x for x in i2.attrs['gene_symbols'].split(':') if x])
	if(not (g1 & g2)):
		sys.stdout.write("%s%s" % (str(i1), str(i2)))
		passed += 1;
	else:
		removed += 1;
		
		
sys.stderr.write("\npassed: %d\nremoved: %d\n\n" % (passed, removed))
		
