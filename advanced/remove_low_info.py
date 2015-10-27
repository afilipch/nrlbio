#! /usr/lib/python
'''Removes reads with chunks of low information content'''
import argparse;
import sys;

from Bio import SeqIO;

from nrlbio.generators import generator_fastq;
from nrlbio.sequencetools import chunk_entropy, entropy

parser = argparse.ArgumentParser(description='Removes reads with chunks of low information content');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to fastq file");
parser.add_argument('-l', '--length', nargs = '?', default = 16, type = int, help = "Length of chunks. Default 16");
parser.add_argument('-w', '--whole', nargs = '?', default=False, const=True, type = bool, help = "If set, entropy is calculated for the whole sequence, to speed up the algorithm");
parser.add_argument('-me', '--min_entropy', nargs = '?', default=1.4, type = float, help = "minimum entropy requiremt for the worst(in sense of entropy) chunk in the read sequence. Default 1.4");
parser.add_argument('--order', nargs = '?', default = 1, type = int, help = "Size of nucleotide blocks to calculate entropy. Default 16");
parser.add_argument('-t', '--stype', nargs = '?', choices = ['fasta', 'fastq'], default='fastq', type = str, help = "type of sequencence file");
args = parser.parse_args();

passed, removed = 0, 0;

if(args.whole):
	def centropy(seq):
		return entropy(seq, order=args.order);
else:
	def centropy(seq):
		return chunk_entropy(seq, args.length, step = 1, order = args.order);


if(args.stype == 'fastq'):
	for fastq in generator_fastq(args.path):
		me = centropy(fastq.seq)
		if(me>args.min_entropy):
			print "\n".join([fastq.id, fastq.seq, fastq.sign, fastq.qual])
			passed+=1;
		else:
			removed+=1;
			
elif(args.stype == 'fasta'):
	for seqrecord in SeqIO.parse(args.path, 'fasta'):
		seq = str(seqrecord.seq.upper());
		#e1 = entropy(seq, order=1)
		#if(e1<2):
			#print seq 
			#print
			#print entropy(seq, order=1)
			#print entropy(seq, order=2)
			#print entropy(seq, order=3)
			#print
			#print "_"*140
		me = centropy(seq);
		if(me>args.min_entropy):
			print "\n".join([">%s" % seqrecord.name, seq])
			passed+=1;
		else:
			removed+=1;		
			
			
sys.stderr.write("total reads: %d\npassed reads: %d\nremoved reads: %d\n" % (passed + removed, passed, removed))