#! usr/bin/python
'''this script generates statistics in html format for provided sam/bam files'''
import os;
import sys;
import argparse;

import pysam

from sam import Stat;




parser = argparse.ArgumentParser(description='script generates clash pipeline makefile and creates required folder tree');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "paths to sam/bam files"); 
parser.add_argument('-o', '--output', nargs = '+', default = [], type = str, help = "paths to the output reports. If not set, reports will be generated in current folder, deriving file names from sam/bam files");
parser.add_argument('--sparse_coefficient', nargs = '?', default = 1, type = int, help = "If sparse_coefficient set to 3, then only each 3rd record will be analyzed, saves time in a case of large files")
parser.add_argument('--short_reference', nargs = '?', const = True, default = False, type = str, help = "If set, provides additional statistics. Useful for sam files derived from mapping to short reference(miRNAs, piRNAs, enhancers, etc.)")
args = parser.parse_args();

#parse options
if(args.output):
	if(len(args.output) == len(args.path)):
		report_names = args.output;
	else:
		sys.exit('number of output files have to be equal to the number of input ones, or \'--output\' should not be set\n')
else:
	report_names = [os.path.basename(x) for x in args.path]
############################################################################################################################

for path, name in zip(args.path, report_names):
	samfile = pysam.Samfile(path)
	stat = Stat(name=name);	
	stat.fill_stat_sheet(samfile.fetch(until_eof=True), short_reference = args.short_reference, detailed = True, sparse_coefficient = args.sparse_coefficient)
	stat.tohtml(output = name, template = "statistic_tables.html", top_entries = 20);