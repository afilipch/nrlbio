#! /usr/lib/python
'''converts SAM/BAM file into fastq masking mapped portion of the read and everything upstream'''
import argparse

import pysam;






parser = argparse.ArgumentParser(description='converts SAM/BAM file into fastq masking mapped portion of the read and everything upstream');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
parser.add_argument('-b', '--bias', nargs = '+', default = [], type = str, help = "accounts for RNAse bias. Masked portion of nucleotides has be enclosed by one of provided nucleotides");
parser.add_argument('--pad', nargs = '?', default = 'N', type = str, help = "Symbol to pad mapped part with");
args = parser.parse_args();


pad = args.pad;
# open output and input sam/bam files
samfile = pysam.Samfile(args.path)
#print len(samfile.lengths)
for ar in samfile.fetch(until_eof=True):
	pos = ar.qend

	if(args.bias):
		
		if(samfile.lengths[ar.tid] != len(ar.query)):
			for i, n in enumerate(ar.query[::-1]):
				if(n in args.bias):
					pos = pos - i;
					break;
		else:
			pass
	else:	
		pass;
	
	padding = pad*pos
	seq = "".join((padding, ar.seq[pos:]));
	qual = ar.qual[pos-len(padding):]
	print "\n".join(['@'+ar.qname, seq, "+", qual])