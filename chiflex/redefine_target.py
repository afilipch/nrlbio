#! /usr/lib/python
'''Redefines genomic regions of miRNA targets, based on seed-anchored hybridization'''

import argparse
import sys
import copy
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import reverse_complement

from nrlbio.generators import generator_doublebed
from nrlbio.rnahybrid import get_rnahybrid
from nrlbio.pybedtools_extension import interval2seq
from nrlbio.mirna import MODES_ORDER, fasta2mirnas


parser = argparse.ArgumentParser(description='Redefines genomic regions of miRNA targets, based on seed-anchored hybridization');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to doublebed/gff file");
parser.add_argument('-g', '--genome', nargs = '?', required = True, type = str, help = "Path to the target reference (genome), fasta file");
parser.add_argument('-m', '--mir', nargs = '?', required = True, type = str, help = "Path to the miRNA reference (mirBase), fasta file");
parser.add_argument('--pval_cutoff', nargs = '?', default = 0.1, type = float, help = "Max p value required for additional(split) interactions to be output");
parser.add_argument('--upstream', nargs = '?', default = 10, type = int, help = "Maximum upstream extension of the target site allowed");
parser.add_argument('--downstream', nargs = '?', default = 10, type = int, help = "Maximum downstream extension of the target site allowed");
parser.add_argument('--addlength', nargs = '?', default = 25, type = int, help = "Length of the target site");
args = parser.parse_args();



mirnas = fasta2mirnas(args.mir, seed_start=1, seed_stop=7)
genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))


def seqgenerator(path):
	for i1, i2 in generator_doublebed(path):
		if(i2.strand == '+'):
			i2.start = i2.start - args.upstream
			i2.end = i2.end + args.downstream
		else:
			i2.start = i2.start - args.downstream
			i2.end = i2.end + args.upstream	
		
		tseq = interval2seq(i2, genome);	
		i2.attrs['seq'] = tseq

		yield i1, i2
		
#def get_mode(targetseq, mirna):
	#temp = mirna.find_fixed_types(targetseq)
	#for mode in MODES_ORDER:
		#if(temp[mode]):
			#return mode
	#else:
		#return 'none'


def seed_aware(targetseq, mirna, addlength):
	fp = mirna.find_fixed_positions(targetseq);
	for mode in MODES_ORDER:
		for p in fp.get(mode, []):
			return -30, 0, max(p - addlength + 10, 0), mode
	else:
		energy, stub, stub, pval, hybstart = get_rnahybrid(targetseq, mirna.seq, extended=True)
		return energy, pval, hybstart, 'none'


def split_target(targetseq, mirna, addlength, correction, pval_cutoff):
	#energy, stub, stub, pval, hybstart = get_rnahybrid(targetseq, mirseq, extended=True);
	energy, pval, hybstart, mode = seed_aware(targetseq, mirna, addlength);	
	
	hybend = hybstart+addlength
	aintervals  = [];
	if(pval <= pval_cutoff):
		if(hybstart>=addlength):
			aintervals.append((targetseq[:hybstart], correction));
		if(len(targetseq)-hybend>=addlength):
			aintervals.append((targetseq[hybend:], correction+hybend));
		return (targetseq[hybstart:hybend], correction+hybstart, energy, mode), aintervals
	else: 
		return None, []
	

def find_distinct_targets(i1, i2, addlength, pval_cutoff):
	targetseq = i2.attrs['seq']
	mirna = mirnas[i1.chrom]
	i1.attrs['seq'] = mirna.seq
	tsite, aintervals = split_target(targetseq, mirna, args.addlength, 0, 2);
	tsites = [tsite];
	
	while(aintervals):
		tais = [];
		for ai in aintervals:
			#print ai;
			tsite, tai = split_target(ai[0], mirna, addlength, ai[1], pval_cutoff);
			if(tsite):
				tsites.append(tsite);
				tais.extend(tai);
		aintervals = tais;
	
	return tsites


def adjust_coordinates(interval, correction):
	l = len(interval.attrs['seq'])
	if(interval.strand == '+'):
		interval.start = interval.start + correction;
		interval.end = interval.start+l
	else:
		interval.end = interval.end - correction
		interval.start = interval.end-l;


modes = defaultdict(int)

for i1, i2 in seqgenerator(args.path):
	mirna = mirnas[i1.chrom];
	
	tsites = find_distinct_targets(i1, i2, args.addlength, args.pval_cutoff)
	if(len(tsites)==1):
		tsite = tsites[0]
		modes[tsite[3]]+=1
		#modes[get_mode(tsite[0], mirna)] +=1
		i2.attrs['seq'] = tsite[0]		
		adjust_coordinates(i2, tsite[1]);
		sys.stdout.write("%s%s" % (i1, i2));
	else:
		nuniq = float(i1.attrs['n_uniq'])/len(tsites)
		for tsite in tsites:
			modes[tsite[3]]+=1
			#modes[get_mode(tsite[0], mirna)] +=1
			tempinterval = copy.copy(i2);
			tempinterval.attrs['n_uniq'] = str(nuniq);
			tempinterval.attrs['seq'] = tsite[0]
			adjust_coordinates(tempinterval, tsite[1]);
			sys.stdout.write("%s%s" % (i1, tempinterval));
			
			
total = float(sum(modes.values()))
sys.stderr.write("mode\tsignal\tsignal_fraction\n")
om = ['none'] + list(MODES_ORDER) 
for mode in om:
	sys.stderr.write("%s\t%d\t%1.2f\n" % (mode, modes[mode], modes[mode]/total))
		
		
	
	
	
	
	



	
	
	
	
	
	
	
	
	
	
	