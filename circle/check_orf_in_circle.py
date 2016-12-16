#! /usr/bin/python
'''Checks for an ORFs troughout the given circles'''


import sys;
import argparse
from collections import defaultdict

from pybedtools import BedTool, Interval
from Bio import SeqIO

from nrlbio.sequencetools import multifind


parser = argparse.ArgumentParser(description='Checks for an ORFs troughout the given circles');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the circles, bed/gff format");
parser.add_argument('--minlength', nargs = '?', required = True, type = int, help = "Minimum length of a potential ORF");
args = parser.parse_args();

START = 'ATG'
STOP = {'TAG', 'TAA', 'TGA'}

def seq2frames(seq, minlength):
	starts = multifind(seq, START, overlap = True);
	#print len(seq)
	#print seq.find('TGA')
	frames = []
	for start in starts:
		localseq = seq[start+3:] + seq[:start]
		localstop = 0;
		for cs in range(0,len(localseq), 3):
			if(localseq[cs:cs+3] in STOP):
				localstop = cs+start+3;
				#break;
				if(localstop - start>minlength):
					if(localstop>=len(seq)):
						frames.append((start, localstop-len(seq)))
					else:
						frames.append((start, localstop))
	return frames;		
			

#test = "ATGAAAAGGGATGTGTGTGCCCCCAAATGCATGGCCCGG"
#print seq2frames(test, args.minlength)

for seqrecord in SeqIO.parse(args.path, 'fasta'):
	seq = str(seqrecord.seq.upper()).replace('U', 'T');
	frames = seq2frames(seq, args.minlength);
	for c, frame in enumerate(frames):
		print "%s\t%d\t%d\t%s:orf%d" % (seqrecord.name, frame[0], frame[1], seqrecord.name, c);