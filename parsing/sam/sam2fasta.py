#! /usr/lib/python
'''creates fasta/fastq file from sam/bam file'''
import argparse
import sys;

import pysam;

parser = argparse.ArgumentParser(description='creates fasta/fastq file from sam/bam file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the sam/bam file");
parser.add_argument('--format', nargs = '?', choices = ['fasta', 'fastq'], default='fasta', type = str, help = "format of output file");
args = parser.parse_args();

samfile = pysam.Samfile(args.path)

if(args.format == 'fasta'):
	for aligned_read in samfile.fetch(until_eof=True):
		sys.stdout.write(">%s\n%s\n" % (aligned_read.qname, aligned_read.seq));
elif(args.format == 'fastq'):
	for aligned_read in samfile.fetch(until_eof=True):
		sys.stdout.write("@%s\n%s\n+\n%s\n" % (aligned_read.qname, aligned_read.seq, aligned_read.qual));