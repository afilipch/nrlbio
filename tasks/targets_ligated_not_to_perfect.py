#! /usr/lib/python
'''Script answers to the question: How many clusters have a perfect seed match for one of the top N expressed miRNAs families, but were actually found ligated to another miRNA?'''

import argparse
import os
import sys

from pybedtools import BedTool
from collections import defaultdict

from nrlbio.mirna import fasta2mirnas, assign_expression, mirnas2families, find_family
from nrlbio.pyplot_extension import histogram

parser = argparse.ArgumentParser(description='Script answers to the question: How many clusters have a perfect seed match for one of the top N expressed miRNAs families, but were actually found ligated to another miRNA?');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to interaction.bed file");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "path to the file with miRNAs in fasta format");
parser.add_argument('--expr', nargs = '?', required=True, type = str, help = "path to the expression file (mirid expression tsv format)");
args = parser.parse_args();

seed_start = 1;
seed_stop = 7


def bound2other(mirid, families, tseq):
	lfam=find_family(mirid, families)
		
		
	for fam in families:
		if(fam!=lfam and fam.match in tseq):
			return 1;
	else:
		return 0;

mirdict = fasta2mirnas(args.mir, seed_start, seed_stop);
assign_expression(args.expr, mirdict, sep="\t");
families = mirnas2families(mirdict.values())
families.sort(key=lambda x: x.expression, reverse = True)

result = defaultdict(int);

total_int=0;
for interval in BedTool(args.path):
	total_int += 1;
	mirid, tseq = interval[6].split(",")[0], interval[8];
	for i in range(1,11):
		result[i] += bound2other(mirid, families[:i], tseq)
		
		
histogram(result, title='interactions that have perfect seed match to another miRNA from top expressed families', ylabel='number of interactions(total %d)' % total_int, xlabel='number of top expressed families', xticks=range(1, 11), xticklabels=None, xticksrotation = 0, output='targets_ligated_not_to_perfect.pdf', color='skyblue', align=u'left', rwidth=0.5)
	
	
	#if(not find_family(mirid, families)):
		#print mirid;
	
