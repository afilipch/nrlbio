#! /usr/lib/python
'''Script answers to the questions: How many clusters have perfect seed matches for >1 of the expressed miRNAs (families); How many clusters do not have any seed match for the top 10(?) expressed miRNAs (families)?'''

import argparse
import os
import sys

from Bio import SeqIO 
from collections import defaultdict

from nrlbio.mirna import fasta2mirnas, assign_expression, mirnas2families, find_family
from nrlbio.pyplot_extension import histogram

parser = argparse.ArgumentParser(description='Script answers to the questions: How many clusters have perfect seed matches for >1 of the expressed miRNAs (families)?; How many clusters do not have any seed match for the top 10(?) expressed miRNAs (families)?');	
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to clusters.fa file");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "path to the file with miRNAs in fasta format");
parser.add_argument('--expr', nargs = '?', required=True, type = str, help = "path to the expression file (mirid expression tsv format)");
args = parser.parse_args();

seed_start = 1;
seed_stop = 7
maxbin=11


def get_num_of_matches(families, tseq):
	nm = 0;
	for fam in families:
		if(fam.match in tseq):
			nm += 1;
	return min(nm, maxbin);


mirdict = fasta2mirnas(args.mir, seed_start, seed_stop);
assign_expression(args.expr, mirdict, sep="\t");
families = mirnas2families(mirdict.values())
families = filter(lambda x: x.expression>10, families)
families.sort(key=lambda x: x.expression, reverse = True)

result = defaultdict(int);
result10 = defaultdict(int)
#print len(families)


total_clusters= 0;
for seqrecord in SeqIO.parse(args.path, "fasta"):
	tseq = str(seqrecord.seq.upper())
	if(len(tseq)<200):
		total_clusters += 1;
		result[get_num_of_matches(families, tseq)] +=1
		result10[get_num_of_matches(families[:10], tseq)] +=1
	
#print result;		
histogram(result10, title='number of seed matches to top 10 miRNA families', ylabel='number of clusters(total %d)' % total_clusters, xlabel='number of seed-mathces', xticks=range(5), xticklabels=None, xticksrotation = 0, output='num_of_seed_matches_10_families.pdf', color='skyblue', align=u'left', rwidth=0.5)
	
	
histogram(result, title='number of seed matches to all miRNA families', ylabel='number of clusters(total %d)' % total_clusters, xlabel='number of seed-mathces', xticks=range(maxbin+1), 
xticklabels=[str(x) for x in range(maxbin)] + ['>%d' % maxbin], xticksrotation = 0, output='num_of_seed_matches_all_families.pdf', color='skyblue', align=u'left', rwidth=0.5)	
	#if(not find_family(mirid, families)):
		#print mirid;
	
