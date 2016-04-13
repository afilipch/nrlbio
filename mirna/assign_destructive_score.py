#! /usr/lib/python
'''Assignes to each miRNA:target interaction destructive score. Does not produce html report'''

import argparse
import os
import sys
import math
from collections import defaultdict
from multiprocessing import Pool

from Bio import SeqIO;
from pybedtools import BedTool

from nrlbio.generators import generator_doublebed
from nrlbio.mirna import slicing_score, destructive_score
from nrlbio.rnahybrid import get_rnahybrid


system_choices = ['hg19', 'hg38', 'mm9', 'mm10', 'ce6', 'circ']
sys2rhsys = {'hg19': '3utr_human', 'hg38': '3utr_human', 'mm9': '3utr_human', 'mm10': '3utr_human', 'ce6': '3utr_worm', 'ce10': '3utr_worm', 'circ': '3utr_human'}

parser = argparse.ArgumentParser(description='Assignes to each miRNA:target interaction destructive score. Does not produce html report');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the miRNA:target interactions, gff/bed file");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "Path to the miRNAs, fasta format");
parser.add_argument('--ftype', nargs = '?', default='single', choices=['single', 'double'], type = str, help = "Type of the input file, single or double gff/bed file. Choices [double|single] (default=single)");
parser.add_argument('--system', nargs = '?', required = True, choices = system_choices, type = str, help = "Genome system. Can be set to %s" % "|".join(system_choices))

parser.add_argument('--threads', nargs = '?', default = 2, type = int, help = "Number of threads to use")
parser.add_argument('--minscore', nargs = '?', default = 20.0, type = float, help = "Only the regions with destructive score greater or equal to [--minscore] will be selected")
args = parser.parse_args();

rhsys = sys2rhsys[args.system]

mirdict = {};
for seqrecord in SeqIO.parse(args.mir, "fasta"):
	mirdict[seqrecord.id] = str(seqrecord.seq.upper())
	
	

def generator_single(path):
	for interval in BedTool(args.path):
		yield interval, interval.attrs['seq'], interval.attrs['mseq'], float(interval.attrs['n_uniq'])
	
	
def generator_double(path):
	for count, (i1, i2) in enumerate(generator_doublebed(args.path)):
		i2.attrs['mseq'] = mirdict[i1.chrom]
		i2.attrs['mirid'] = i1.chrom
		i2.name = "ds%d" % (count+1)
		yield i2, i2.attrs['seq'], i2.attrs['mseq'], float(i2.attrs['n_uniq'])
		
	
if(args.ftype == 'single'):
	cgenerator = generator_single(args.path);
elif(args.ftype == 'double'):
	cgenerator = generator_double(args.path);
else:
	sys.exit('ftype is not set correctly\n')
	
	
def getscore(a):
	cid, tseq, mseq, n_uniq = a
	energy, pattern, basepairing, pval = get_rnahybrid(tseq, mseq, system = rhsys, extended=True);
	return cid, destructive_score(basepairing) + math.log(n_uniq, 2)


#for entry in cgenerator:
	#print getscore(entry)


pool = Pool(args.threads);
res = []
for ri, dscore in pool.imap(getscore, cgenerator):
	if(dscore>=args.minscore):
		ri.attrs['dscore'] = "%1.2f"  % dscore;
		sys.stdout.write(str(ri))
		

	
	
	
