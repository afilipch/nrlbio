#! /usr/bin/python
'''Checks for a splice signal around the edges of the provided junctions. Filters if asked''' 
import argparse
import sys;
from collections import defaultdict, Counter

from Bio import SeqIO
from pybedtools import BedTool


from nrlbio.genome_system import seqrecord2seq



parser = argparse.ArgumentParser(description='Checks for a splice signal around the edges of the provided junctions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the linear splice junctions, bed/gff file");
parser.add_argument('--reference', nargs = '?', required = True, type = str, help = "path to the reference(genome) to extract sequences from");
parser.add_argument('--allowed_signal', nargs = '+', default = [], type = str, help = "if provided, only the junctions with given flanks will pass the filter. Pass the arguments in a format: GT,AG GC,AG");
args = parser.parse_args();

allowed_signal = [tuple(x.split(",")) for x in args.allowed_signal]


def check_splice_site(interval, seqrecord):
	if(interval.strand == '+'):
		flank5 = seqrecord2seq(seqrecord, interval.start, interval.start+2, strand='+')
		flank3 = seqrecord2seq(seqrecord, interval.end-2, interval.end, strand='+') 
	elif(interval.strand == '-'):
		flank5 = seqrecord2seq(seqrecord, interval.end-2, interval.end, strand='-') 
		flank3 = seqrecord2seq(seqrecord, interval.start, interval.start+2, strand='-')

	return flank5, flank3
		

flank_count = defaultdict(int);
reference = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))
total = 0;
passed = 0;
for interval in BedTool(args.path):
	seqrecord = reference[interval.chrom];
	flanks = check_splice_site(interval, seqrecord)
	flank_count[flanks]+=1;
	total+=1
	if(not allowed_signal or flanks in allowed_signal):
		sys.stdout.write(str(interval));
		passed += 1;
	else:
		pass;
	
	
sys.stderr.write("total junctions: %d\npassed junctions: %d\nfraction passed %1.5f\n\n" % (total, passed, float(passed)/total));
sys.stderr.write("type of signal\tnum of junctions\n")	
	
for (flank5, flank3), count in flank_count.items():
	sys.stderr.write("%s<->%s\t%d\n" % (flank5, flank3, count))