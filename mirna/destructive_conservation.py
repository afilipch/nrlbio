#! /usr/lib/python
'''Checks how well is destructive pattern conserved for the given regions'''

import argparse
import os
import sys
from collections import defaultdict

from pybedtools import BedTool, Interval
from Bio import SeqIO

from nrlbio.mirna import destructive_score, mirfasta2conservation
from nrlbio.rnahybrid import get_rnahybrid
from nrlbio.generators import targets_generator

parser = argparse.ArgumentParser(description='Tool for destructive miRNA sites prediction');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genomic intervals to check binding conservation, gff/bed file");
#parser.add_argument("--system", nargs = '?', type=str, required=True, help="model system/reference species (hg19|dm6|...)")
#parser.add_argument("--maf", nargs = '?', type=str, required=True, help="path to the indexed MAF files")
#parser.add_argument("--genomes", nargs = '?', type=str, required=True, help="path to the genomes")
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "Path to the miRNAs conservation file, fasta format");
parser.add_argument('--mir2ucsc', nargs = '?', required=True, type = str, help = "Path to the table which connects miRNA specie names with those from MAF files, tsv file");
parser.add_argument('--mafasta', nargs = '?', required=True, type = str, help = "Path to the targets conservation file, fasta format");
parser.add_argument('--refspecie', nargs = '?', required=True, type = str, help = "Name of the reference specie genome [hg19, mm10,..]");
parser.add_argument('--minscore', nargs = '?', default=18.0, type = float, help = "Min destructive score required to count an interaction as conserved");
args = parser.parse_args();


def targets_generator(consfasta):
	sequences = {}
	for seqrecord in SeqIO.parse(consfasta, 'fasta'):
		a = seqrecord.id.split("|");
		if(len(a) == 6):
			if(sequences):
				yield name, mirid, sequences
			sequences = {}
			name, mirid = a[4], a[5]
		else:
			sequences[seqrecord.id] = str(seqrecord.seq.upper()).replace('U', 'T')
	else:
		yield name, mirid, sequences
		
		

#################################################################################################################
#get translational table from MirBase to UCSC genome names
translational_table = {}
with open(args.mir2ucsc) as f:
	for l in f:
		a = l.strip().split("\t");
		translational_table[a[1]] = a[0]
		
		
mir2cons, stub = mirfasta2conservation(args.mir, translational_table=translational_table)

#for v in mir2cons.values():
	#print v.keys();
	

#get conservation destrucive scores for aligned regions
coord2score = {}
for name, mirid, targets in targets_generator(args.mafasta):
	score = 0;
	scores = []
	targets.pop(args.refspecie, None)
	tdict = mir2cons[mirid]
	
	for specie, mirseq in tdict.items():
		tseq = targets.get(specie, None)
		if(tseq):
			energy, pattern, basepairing, pval, pos = get_rnahybrid(tseq.replace("-", ""), mirseq, extended=True);
			dscore = destructive_score(basepairing)
			if(dscore>=args.minscore):
				score+=1;
			scores.append(dscore)
		
	coord2score[name] = (score, ",".join([str(x) for x in scores]))
	
	
for interval in BedTool(args.path):
	score, scores = coord2score.get(interval.name, (0, ''))
	interval.attrs['cons_dscore'] = str(score)
	interval.attrs['cons_dscores'] = scores
	if(scores):
		sys.stdout.write(str(interval))
		
		
		
		
		
		
		
	
	
	
	
	