#! /usr/lib/python
'''converts SAM/BAM file into fastq'''
import argparse

import pysam;






parser = argparse.ArgumentParser(description='converts SAM/BAM file into fastq');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
args = parser.parse_args();

# open output and input sam/bam files
samfile = pysam.Samfile(args.path)
for ar in samfile.fetch(until_eof=True):
	print "\n".join([ar.qname, ar.seq, "+", ar.qual])