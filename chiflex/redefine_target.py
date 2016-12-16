#! /usr/lib/python
'''Redefines genomic regions of miRNA targets, based on seed-anchored hybridization'''

import argparse
import sys

from Bio import SeqIO

from nrlbio.generators import generator_doublebed
from nrlbio.rnahybrid import get_rnahybrid, get_shuffled_rnahybrid
from nrlbio.pybedtools_extension import interval2seq


parser = argparse.ArgumentParser(description='Redefines genomic regions of miRNA targets, based on seed-anchored hybridization');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to doublebed/gff file");
#parser.add_argument('-f', '--fasta', nargs = '+', required = True, type = str, help = "Path to reference(genome, mirBase, transcriptome) fasta files");
#parser.add_argument('--genome', nargs = '?', required = True, type = str, help = "Path to the reference for targets (genome), fasta format");
#parser.add_argument('--mirna', nargs = '?', required = True, type = str, help = "Path to the reference for miRNA sequences, fasta format");
parser.add_argument('--pval_cutoff', nargs = '?', default = 0.1, type = float, help = "Max p value required for additional(split) interactions to be output");
parser.add_argument('--upstream', nargs = '?', default = 10, type = int, help = "Maximum upstream extension of the target site allowed");
parser.add_argument('--downstream', nargs = '?', default = 10, type = int, help = "Maximum downstream extension of the target site allowed");
args = parser.parse_args();



#def seqrecord_iterator(paths):
	#for path in paths:
		#for seqrecord in SeqIO.parse(path, "fasta"):
			#yield seqrecord;


#reference = SeqIO.to_dict(seqrecord_iterator(args.fasta))


#for i1, i2 in generator_doublebed(args.path):
	#if(i2.strand == '+'):
		#i2.start = i2.start - args.upstream
		#i2.end = i2.end + args.downstream
	#else:
		#i2.start = i2.start - args.downstream
		#i2.end = i2.end + args.upstream	
	
	#tseq = interval2seq(i2, reference);
	#mirseq = reference[i1.chrom];
	
	#i2.attrs['seq'] = tseq
	#i1.attrs['seq'] = str(mirseq.seq)
	
	#sys.stdout.write("%s%s" % (i1, i2))
	


for i1, i2 in generator_doublebed(args.path):
	stub, stub, stub, pval, hybstart = get_rnahybrid(i2.attrs['seq'], i1.attrs['seq'], extended=True);