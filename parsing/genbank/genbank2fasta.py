#! /usr/bin/python
'''script convert pre-selected features(transcripts, exons, CDs) in provided genbank files into fasta file. Fasta file is printed to stdout'''

import sys, os
from Bio import SeqIO
import argparse
from nrlbio.ncbi import seq_record2fasta

parser = argparse.ArgumentParser(description='script convert pre-selected features(transcripts, exons, CDs) in provided genbank files into fasta file. Fasta file is printed to stdout');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to genbank files");
parser.add_argument('-ft', '--feature_type', nargs = '+', default = ['mRNA', 'misc_RNA'], type = str, help = "select only the transcripts of provided types, might be \'exon\', \'gene\', \'mRNA\', \'CDS\', \'misc_RNA\', \'misc_feature\'");
parser.add_argument('-mt', '--miscRNA_type', nargs = '+', default = None, type = str, help = "select only the non-coding transcripts of provided types, might be \'ncRNA\', \'pseudogene\', \'tRNA\', \'snoRNA\', \'miRNA\', \'lincRNA\', \'rRNA\', \'snRNA\'");
args = parser.parse_args();
	
				
for myfile in args.path:
	seq_record = SeqIO.read(myfile, "genbank")
	seq_record2fasta(seq_record, feature_types = args.feature_type, miscRNA_type = args.miscRNA_type);
	