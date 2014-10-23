#! /usr/lib/python
'''produces old-fashion interaction.bed file of mirna and their targets'''
import argparse;
import sys;
import re;

from pybedtools import BedTool
from Bio import SeqIO
from Bio.Seq import reverse_complement




#import pysam;
#from nrlbio.filters_for_sam import *
#from nrlbio.chimera import arlist2chimera
#from nrlbio import chimera

parser = argparse.ArgumentParser(description='produces old-fashion interaction.bed file of mirna and their targets');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the new-style interaction file");
parser.add_argument('-m', '--mirna', nargs = '?', required = True, type = str, help = "path to miRNA fasta file");
parser.add_argument('-s', '--system', nargs = '?', required = True, choices = ['ce6', 'mm9', 'hg19'], type = str, help = "genome ce6|mm9|hg19");
args = parser.parse_args();

def reassign_coordinates(a):
	chrom, strand, start, stop = a[0].split("|")
	start = int(start)
	stop = int(stop)
	
	if(strand == '+'):
		stop = start + int(a[2])
		start = start + int(a[1])		
	else:
		start = stop - int(a[2])
		stop = stop - int(a[1])
	
	return chrom, start, stop, a[3], a[4], strand


interactions = BedTool(args.path)
mirna = SeqIO.to_dict(SeqIO.parse(args.mirna, "fasta"))

s = 0;	
n = 0;
for i in interactions:
	mirid = i[0]
	mirseq = str(mirna[mirid].seq.upper())
	
	
	chrom, start, stop, name, score, strand = reassign_coordinates(i[6:12])
	
	tseq = gsys.genome.get_oriented(chrom, start, end, strand).upper()
	
	print "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, start, end, i[3], i[10], strand, mirid, mirseq, tseq)
	
	if(reverse_complement(mirseq[1:7]) in tseq):
		s+=1;
	else:
		n+=1;

sys.stderr.write("total interactions: %d\n2-7 seed interactions: %d\nnon 2-7 seed interactions: %d\nfraction of 2-7 seeds: %1.5f\n" % (s + n, s, n, s/float(s+n)));		