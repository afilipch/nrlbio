#! /usr/bin/python
'''Selects interactions only for with ids found in provided bed/gff file''' 
import argparse
import sys;

from pybedtools import BedTool

from nrlbio.generators import generator_doublebed



parser = argparse.ArgumentParser(description='Selects interactions only for with ids found in provided bed/gff file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the interactions to select from, double gff");
parser.add_argument('--selection', nargs = '?', type = str, help = "path to the intervals to select ids for filtering, bed/gff file");
args = parser.parse_args();

ids = [x.name.split("|")[0] for x in BedTool(args.selection)]
ids = set(ids);

for i1, i2 in generator_doublebed(args.path):
	if(i1.name.split("|")[0] in ids):
		sys.stdout.write(str(i1));
		sys.stdout.write(str(i2));