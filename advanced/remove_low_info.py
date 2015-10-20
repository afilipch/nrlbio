#! /usr/lib/python
'''Removes reads with chunks of low information content'''
import argparse;
import sys;

from Bio import SeqIO;

from nrlbio.generators import generator_fastq;
from nrlbio.sequencetools import chunk_entropy

parser = argparse.ArgumentParser(description='removes reads with chunks of low information content');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to fastq file");
parser.add_argument('-l', '--length', nargs = '?', default = 14, type = int, help = "length of chunks");
parser.add_argument('-me', '--min_entropy', nargs = '?', default=1.35, type = float, help = "minimum entropy requiremt for the worst(in sense of entropy) chunk in the read sequence");
parser.add_argument('-t', '--stype', nargs = '?', choices = ['fasta', 'fastq'], default='fastq', type = str, help = "type of sequencence file");
args = parser.parse_args();

passed, removed = 0, 0;

if(args.stype == 'fastq'):
	for fastq in generator_fastq(args.path):
		me = chunk_entropy(fastq.seq, args.length, step = 1, order = 1);
		if(me>args.min_entropy):
			print "\n".join([fastq.id, fastq.seq, fastq.sign, fastq.qual])
			passed+=1;
		else:
			removed+=1;
			
elif(args.stype == 'fasta'):
	for seqrecord in SeqIO.parse(args.path, 'fasta'):
		seq = str(seqrecord.seq.upper());
		me = chunk_entropy(seq, args.length, step = 1, order = 1);
		if(me>args.min_entropy):
			print "\n".join([">%s" % seqrecord.name, seq])
			passed+=1;
		else:
			removed+=1;		
			
			
sys.stderr.write("total reads: %d\npassed reads: %d\nremoved reads: %d\n" % (passed + removed, passed, removed))