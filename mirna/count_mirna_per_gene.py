#! /usr/bin/python
'''Counts how many miRNA target sites (found via chimeras) are on each gene'''
import sys;
import argparse
from collections import defaultdict

from pybedtools import BedTool


parser = argparse.ArgumentParser(description='Counts how many miRNA target sites (found via chimeras) are on each gene');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the annotated mirna targets, gff format");
parser.add_argument('--utr3', nargs = '?', default = False, const = True, type = bool, help = "If set, only 3'utr residents are counted");

args = parser.parse_args();


#attrs = args.attributes

#if(args.description):
	#base = ["#chrom", "start", "end", "ID", "score", "strand"]
	#print "\t".join(base + attrs)

#if(args.addchr):
	#def fix_chrom(interval):
		#if(interval.chrom == 'MT'):
			#return 'chrM'
		#else:
			#return "chr%s" % interval.chrom
#else:
	#fix_chrom = lambda x: x.chrom
	
	
genes_mirna_counts = defaultdict(lambda: defaultdict(float))
genes_mirna_mcounts = defaultdict(lambda: defaultdict(float))

for interval in BedTool(args.path):
	gsymbols = interval.attrs['gene_symbols'].split(':')
	biotypes = interval.attrs['biotypes'].split(':')
	transcription = interval.attrs['transcription'].split(':')
	mirid = interval.attrs['mirid'];
	mode = interval.attrs['mode'];
	
	filtered_symbols = []
	
	for gsymbol, biotype, tr in zip(gsymbols, biotypes, transcription):
		if(gsymbol and 'exon' in tr):
			if((not args.utr3) or 'utr3' in biotype):
				filtered_symbols.append(gsymbol);
	
	if(filtered_symbols):
		score = float(interval.attrs['n_uniq'])/len(filtered_symbols);
		if(mode not in ['none', 'mm28', 'm38']):
			mscore = score;
		else:
			mscore = 0.0;
			
		for gsymbol in filtered_symbols:
			genes_mirna_counts[gsymbol][mirid] += score;
			genes_mirna_mcounts[gsymbol][mirid] += mscore;
	#if(len(gsymbols)>1):
		#print biotypes, gsymbols, transcription
		
for gene, d in genes_mirna_counts.iteritems():
	for mirid, count in d.items():
		print "%s\t%s\t%1.1f\t%1.1f" % (gene, mirid, count, genes_mirna_mcounts[gene][mirid])
		
		
