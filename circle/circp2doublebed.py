#! /usr/lib/python
'''Convert Circular Pipeline ouput file into Chiflex doublebed format''' 
import argparse
import sys;

from pybedtools import BedTool

from nrlbio.pybedtools_extension import construct_gff_interval;

parser = argparse.ArgumentParser(description='Convert Circular Pipeline ouput file into Chiflex doublebed format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to Circular Pipeline ouput file ");
parser.add_argument('--length', nargs = '?', default = 50, type = int, help = "assumed length of linked intervals");
args = parser.parse_args();



for i in BedTool(args.path):
	i1 = construct_gff_interval(i.chrom, i.start, i.start+args.length, 'chimera', score=i[7], strand=i.strand, source='cp', frame='.', 
	attrs=[ ('n_uniq', i[6]), ("ID", "|".join((i.name, "0"))) ])
	i2 = construct_gff_interval(i.chrom, i.end-args.length, i.end, 'chimera', score=i[8], strand=i.strand, source='cp', frame='.', 
	attrs=[('n_uniq', i[6]), ("ID", "|".join((i.name, "1"))) ])
	sys.stdout.write(str(i1))
	sys.stdout.write(str(i2))
	
	
	
#first columns of the Circular Pipeline ouput file: chrom start   end     name    n_reads strand  n_uniq  best_qual_A     best_qual_B 