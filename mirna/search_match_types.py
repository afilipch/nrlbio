#! /usr/lib/python
'''Searches for seed-matches types in given fasta sequences'''

import argparse
import os
import sys

from Bio import SeqIO 
from collections import defaultdict, Counter

from nrlbio.mirna import fasta2mirnas, assign_expression, mirnas2families, find_family
from nrlbio.pyplot_extension import histogram

parser = argparse.ArgumentParser(description='Finds seed matches for provided miRNAs. If ezpression of miRNAS is provided, than only seeds for top expressed miRNA families may be detected');	
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sequence to look for seed matches, fasta format");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "path to the miRNAs, fasta format");
args = parser.parse_args();


order = ('m27', 'm28', 'm29', 'm27a', 'm28a', 'm29a')



#get dictionary of mirna.Mirna objects
mirdict = fasta2mirnas(args.mir);


print "%s\t%s\t%s" % ('seq_id', 'mirna_id', "\t".join(order));
for seqrecord in SeqIO.parse(args.path, "fasta"):
	tseq = str(seqrecord.seq.upper())
	for mirid, mirna in mirdict.items():
		types_ = mirna.find_fixed_types(tseq);
		if(any(types_.values())):
			print "%s\t%s\t%s" % (seqrecord.name, mirid, "\t".join([str(types_[x]) for x in order]));
		



	
	
	


	
