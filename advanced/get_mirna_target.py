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

exec("from sequence_data.systems import %s as gsys" % args.system);
pat = re.compile('[:)(]')

interactions = BedTool(args.path)
mirna = SeqIO.to_dict(SeqIO.parse(args.mirna, "fasta"))
#for k, v in mirna.iteritems():
	#print k
	#print v.seq
	#print

	
s = 0;	
n = 0;
for i in interactions:
	mirid = i[0]
	mirseq = str(mirna[mirid].seq.upper())
	
	
	a = re.split(pat, i[6])
	start, end = [int(x) for x in a[1].split("-")]
	chrom = a[0];
	strand = a[2];

	
	if(strand == '+'):
		end = start + int(i[8])
		start = start + int(i[7])		
	else:
		start = end - int(i[8])
		end  = end - int(i[7])	
	
	tseq = gsys.genome.get_oriented(chrom, start, end, strand).upper()
	
	print "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, start, end, i[3], i[10], strand, mirid, mirseq, tseq)
	
	if(reverse_complement(mirseq[1:7]) in tseq):
		s+=1;
	else:
		n+=1;

sys.stderr.write("total interactions: %d\n2-7 seed interactions: %d\nnon 2-7 seed interactions: %d\nfraction of 2-7 seeds: %1.5f\n" % (s + n, s, n, s/float(s+n)));		