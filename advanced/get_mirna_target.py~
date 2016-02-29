#! /usr/lib/python
'''produces old-fashion interaction.bed file of mirna and their targets'''
import argparse;
import sys;
import re;

from pybedtools import BedTool, Interval
from Bio import SeqIO
from Bio.Seq import reverse_complement

from nrlbio.generators import generator_doublebed
from nrlbio.pybedtools_extension import interval2seq


#sys.exit()

#import pysam;
#from nrlbio.filters_for_sam import *
#from nrlbio.chimera import arlist2chimera
#from nrlbio import chimera

parser = argparse.ArgumentParser(description='produces old-fashion interaction.bed file of mirna and their targets');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the new-style(double bed/gff file) interaction file");
parser.add_argument('-m', '--mirna', nargs = '?', required = True, type = str, help = "path to a miRNA fasta file");
parser.add_argument('-f', '--fasta', nargs = '?', required = True, type = str, help = "path to a reference(genome) fasta file");
args = parser.parse_args();

#exec("from sequence_data.systems import %s as gsys" % args.system);

def reassign_coordinates(interval):
	chrom, strand, start, stop = interval.chrom.split("|")[:4]
	start = int(start)
	stop = int(stop)
	
	if(strand == '+'):
		stop = start + interval.stop
		start = start + interval.start
	else:
		start = stop - interval.stop
		stop = stop - interval.start
	
	return chrom, start, stop, interval.name, interval.score, strand


#interactions = BedTool(args.path)
mirna = SeqIO.to_dict(SeqIO.parse(args.mirna, "fasta"))
reference = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

s = 0;	
n = 0;
for i1, i2 in generator_doublebed(args.path):
	mirseq = str(mirna[i1.chrom].seq.upper())
	
	
	chrom, start, stop, name, score, strand = reassign_coordinates(i2)
	interval = Interval(chrom, start, stop, name=name, score=score, strand=strand)
	tseq = interval2seq(interval, reference)
	
	print "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, start, stop, i2.name[:-2], i2.score, strand, i1.chrom, mirseq, tseq)
	
	if(reverse_complement(mirseq[1:7]) in tseq):
		s+=1;
	else:
		n+=1;

sys.stderr.write("total interactions: %d\n2-7 seed interactions: %d\nnon 2-7 seed interactions: %d\nfraction of 2-7 seeds: %1.5f\n" % (s + n, s, n, s/float(s+n)));		