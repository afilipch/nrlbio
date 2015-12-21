#! /usr/lib/python
'''Looks for interactions and interactors with the highest number of uniq read support'''

import argparse
import os
import sys
from collections import defaultdict

from pybedtools import BedTool

from nrlbio.generators import generator_doublebed
from nrlbio.mirna import mirnas2families, fasta2mirnas



parser = argparse.ArgumentParser(description='Looks for interactions and interactors with the highest number of uniq read support');	
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the interactions, doublegff format");
parser.add_argument('--collapsed', nargs = '+', type = str, help = "Path to the file(s) with collpased nonunique mappings. If set, unique reads count for each interaction will be split between the regions which are putative origin for one of the interactors. For the case of doubleproject 2 files should be provided in an order mapping was done (that is first miRNA, then targets)");
parser.add_argument('--nosplit', nargs = '?', default=False, const=True, type = bool, help = "If set, read counts will be added to all of them. Otherwise read counts will be split between collapsed interactors.");
parser.add_argument('--circles', nargs = '?', type = str, help = "Path to the circbase annotation(if mapping was done to circles)");
parser.add_argument('--genes', nargs = '?', default=False, const=True, type = bool, help = "Merge interactors on gene level");
parser.add_argument('--mirfam', nargs = '?', type = str, help = "Path to miRNA sequences, fasta format. If set, interactors will be merged on miRNA family level");
args = parser.parse_args();


interactions_counts = defaultdict(lambda: defaultdict(int))

for i1, i2 in generator_doublebed(args.path):
	counts = int(i1.attrs['n_uniq']);
	interactions_counts[(i1.chrom, i2.chrom)]['n_uniq'] += counts;
	interactions_counts[(i1.chrom, i2.chrom)]['n_spots'] += 1;
	
if(args.circles):
	circles = {}
	for interval in BedTool(args.circles):
		circles[interval.name] = interval;
		
if(args.mirfam):
	fam_dict = {};
	for family in mirnas2families(fasta2mirnas(args.mirfam).values()):
		for name in family.names:
			fam_dict[name] = family.name
	
	fam_counts = defaultdict(lambda: defaultdict(int))
	for (chrom1, chrom2), dcount in interactions_counts.items():
		fam_counts[(fam_dict[chrom1], chrom2)]['n_uniq'] += dcount['n_uniq'];
		fam_counts[(fam_dict[chrom1], chrom2)]['n_spots'] += dcount['n_spots'];
	interactions_counts = fam_counts;
	
	
print "miRNA\tCircle\tread_support_interactions\tnumber_of_target_sites\tgene_name\tread_support_circular_jucntion\tlength_of_circle"
	
	
if(args.genes):
	gene_counts = defaultdict(lambda: defaultdict(int))
	for (chrom1, chrom2), dcount in interactions_counts.items():
		circle = circles[chrom2];
		gene = circle[8]
		gene_counts[(chrom1, gene)]['n_uniq'] += dcount['n_uniq'];
		gene_counts[(chrom1, gene)]['n_spots'] += dcount['n_spots'];
		gene_counts[(chrom1, gene)]['length'] = circle[6];
		gene_counts[(chrom1, gene)]['rcount'] += int(circle.score);
	
	for (chrom1, chrom2), dcount in sorted(gene_counts.items(), key = lambda x: x[1]['n_uniq'], reverse=True):
		print "%s\t%s\t%d\t%d\t%s\t%s\t%s" % (chrom1, chrom2, dcount['n_uniq'], dcount['n_spots'], chrom2, dcount['rcount'], dcount['length'])	

else:
	for (chrom1, chrom2), dcount in sorted(interactions_counts.items(), key = lambda x: x[1]['n_uniq'], reverse=True):
		circle = circles[chrom2] 
		print "%s\t%s\t%d\t%d\t%s\t%s\t%s" % (chrom1, chrom2, dcount['n_uniq'], dcount['n_spots'], circle[8], circle.score, circle[6])