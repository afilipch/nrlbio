#! /usr/bin/python
'''Produces fasta file of spliced transcripts for a given gff3 file'''
import sys;
import argparse
from collections import defaultdict

from pybedtools import BedTool
from Bio import SeqIO
from Bio.Seq import reverse_complement

from nrlbio.sequencetools import splitstring


parser = argparse.ArgumentParser(description='Produces fasta file of spliced transcripts for a given gff3 file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genome system annotation, gff3 format");
parser.add_argument('-g', '--genome', nargs = '?', required = True, type = str, help = "Path to the genome, fasta file");
args = parser.parse_args();


transcripts = defaultdict(list)
for interval in BedTool(args.path):
	if(interval[2] == 'exon'):
		transcripts[interval.attrs['Parent'].split(':')[1]].append(interval);
		

genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
for k in genome.keys():
	if(k == 'chrM'):
		genome['MT'] = genome.pop(k)
	else:
		a = k.split("_")
		if(len(a) == 1):
			genome[k[3:]] = genome.pop(k)
		else:
			genome.pop(k);
			
print genome.keys()
		


for name, exons in transcripts.iteritems():
	splice = [];
	exons.sort(key = lambda x: x.start);
	seqrec = genome.get(exons[0].chrom)
	
	if(seqrec):
		for exon in exons:
			splice.append( str(seqrec[exon.start:exon.end].seq.upper()) )
		
		seq = ''.join(splice)
		if(exons[0].strand == '-'):
			seq = reverse_complement(seq)
			
		print ">%s\n%s" % (name, "\n".join(splitstring(seq, 60)));
			
			
			
			