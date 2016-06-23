#! /usr/bin/python
'''Checks how well is destructive pattern conserved for the given regions with the special emphasis on central bulge site'''
import argparse
import os
import sys
import copy
from collections import defaultdict

from pybedtools import BedTool, Interval
from Bio import SeqIO

from nrlbio.mirna import destructive_score, mirfasta2conservation
from nrlbio.rnahybrid import get_rnahybrid
from nrlbio.generators import targets_generator


parser = argparse.ArgumentParser(description='Tool for destructive miRNA sites prediction');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the genomic intervals to check binding conservation, gff/bed file");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "Path to the miRNAs conservation file, fasta format");
parser.add_argument('--mir2ucsc', nargs = '?', required=True, type = str, help = "Path to the table which connects miRNA specie names with those from MAF files, tsv file");
parser.add_argument('--mafasta', nargs = '?', required=True, type = str, help = "Path to the targets conservation file, fasta format");
parser.add_argument('--refspecie', nargs = '?', required=True, type = str, help = "Name of the reference specie genome [hg19, mm10,..]");
parser.add_argument('--minscore', nargs = '?', default=25.0, type = float, help = "Min destructive score required to count an interaction as conserved, and push it to the central bulge conservation analysis");
args = parser.parse_args();

#mm10tseq = 'CACCATACACT-TTGGGAACTTGCCAAGCTCTCTCCAGCCTCTGTGC---C'
#rn5tseq =  'CACCATACACT-TTGGGAACTTGCCGAGCTCTCTCCAGCCTCTGTGC---C'

#mm10mseq = 'TGGAGAGAAAGGCAGTTCCTGA'
#rn5mseq = 'TGGAGAGAAAGGCAGTTCCTGA'
#tseq1 = "GAACA----TTCTT" + "GT-CA--TGCGG----CTAGCA--ACCAA-AGAG" + "CAGGA--GATGC"

#mseq1 = "TCTTTGGTTATCTAGCTGTATGA"


def get_hybrid_on_mirna(basepairing, mseq, tseq, pos):
	curbp = copy.copy(basepairing)
	
	#print 
	#print "-"*120
	#print
	#print curbp, pos

	if(curbp[1][0] == ' ' and curbp[3][0] == ' '):
		curbp = [x[1:] for x in curbp]
	if(curbp[1][-1] == ' ' and curbp[3][-1] == ' '):
		curbp = [x[:-1] for x in curbp]
		
	if(curbp[1][0] == ' ' and curbp[2][0] == ' '):
		tpos = -2;
	else:
		tpos = -1;
		

	#print tseq
	#print curbp

	mpos = 0;
	mlen = len(mseq);
	matchdict = defaultdict(list);
	
	added = set()
	
	for tb, th, mh, mb in zip(*curbp):
		#print tb, th, mh, mb, tseq[tpos+pos]
		if(th != ' '):
			tpos += 1;
			mpos += 1;
			if(tpos+pos not in added):
				matchdict[mlen-mpos].append((True, tpos+pos, tseq[tpos+pos]))
				added.add(tpos+pos)
		else:
			if(tb != ' '):
				tpos += 1;
			if(mb != ' '):
				mpos += 1;
			if(tpos+pos not in added):
				matchdict[mlen-mpos].append((False, tpos+pos, tseq[tpos+pos]));
				added.add(tpos+pos)
		#print tb, th, mh, mb, tseq[tpos+pos]
		
	return matchdict


def lift_coordinates(tseq, matchdict):
	ngaps = 0;
	shift = [];
	for p, n in enumerate(tseq):
		if(n!='-'):
			shift.append(p)
	
	lifted = {}
	for k, v in matchdict.items():
		lifted[k] = [(x[0], shift[x[1]], x[2]) for x in v]
		
	return lifted


#get_hybrid_on_mirna(tseq1, mseq1)
#get_gapped_hybrid_on_mirna(tseq1, mseq1)
#refmatch = get_gapped_hybrid_on_mirna(mm10tseq, mm10mseq)

def get_bulge_score(cmatch, refmatch, bulge=(8,12)):
	#for k, v in refmatch.items():
		#print k, v
	bscore = 0.0;
	norm_score = 0;
	for pos in range(bulge[0], bulge[1]):
		#print pos, refmatch[pos], cmatch[pos]
		ref = dict([ (x[1], (x[0], x[2])) for x in refmatch.get(pos, [])])
		cur = dict([ (x[1], (x[0], x[2])) for x in cmatch.get(pos, [])])
		
		localscores = []
		for rpos, rhyb in ref.items():
			if(not rhyb[0]):
				chyb = cur.get(rpos);
				if(chyb and chyb[0]==False):
					if(chyb[1] != rhyb[1]):
						localscores.append(1);
					else:
						localscores.append(0);

		if(localscores):
			norm_score += 1;
			bscore += sum(localscores)/float(len(localscores))
			
	if(norm_score):
		return bscore/norm_score, norm_score
	else:
		return 0.0, 0;


def get_match_score(cmatch, refmatch, seed=(1,9), down=12):
	score = 0.0;
	norm_score = 0;
	for pos in list(range(seed[0], seed[1])) + list(range(down, min(len(refmatch), len(cmatch)))):
		#print pos, refmatch[pos], cmatch[pos]
		ref = refmatch.get(pos, [])
		cur = cmatch.get(pos, [])
		
		if(len(ref)==1 and ref[0][0]):
			norm_score+=1
			if(len(cur)==1 and cur[0][0] and ref[0][1] == cur[0][1]):
				score += 1;
	if(norm_score):
		return score/norm_score
	else:
		return 0.0;
	
	
#get_bulge_score(rn5tseq, rn5mseq, refmatch)	


#################################################################################################################
#get translational table from MirBase to UCSC genome names
translational_table = {}
with open(args.mir2ucsc) as f:
	for l in f:
		a = l.strip().split("\t");
		translational_table[a[1]] = a[0]
		
		
mir2cons, refmir = mirfasta2conservation(args.mir, translational_table=translational_table)

#for k, v in mir2cons.items():
	#print k;
	#print v
	#print
	#print "-"*140
	

#get conservation destrucive scores for aligned regions
coord2score = {}
for name, mirid, targets in targets_generator(args.mafasta):
	scores = [];
	bscores = [];
	mscores = [];
	bls = [];
	target = targets.pop(args.refspecie, None)
	mdict = mir2cons[mirid]

	mseq = refmir[mirid]
	wogaps = target.replace('-', '')
	energy, pattern, basepairing, pval, pos = get_rnahybrid(wogaps, mseq, extended=True);
	dscore = destructive_score(basepairing)
	
	if(dscore>=args.minscore):
		matchdict = get_hybrid_on_mirna(basepairing, mseq, wogaps, pos);
		refmatch = lift_coordinates(target, matchdict);
		#for k, v in refmatch.iteritems():
			#print k, v
	
		for specie, mseq in mdict.items():
			tseq = targets.get(specie, None)
			if(tseq):
				wogaps = tseq.replace('-', '')
				energy, pattern, basepairing, pval, pos = get_rnahybrid(wogaps, mseq, extended=True);
				dscore = destructive_score(basepairing)
				if(dscore>=args.minscore):
					matchdict = get_hybrid_on_mirna(basepairing, mseq, wogaps, pos);
					cmatch = lift_coordinates(tseq, matchdict);
					bscore, blength = get_bulge_score(cmatch, refmatch, bulge=(8,13))
					mscore = get_match_score(cmatch, refmatch, seed=(1,9), down=12)
					
					scores.append(dscore)
					bscores.append(bscore)
					mscores.append(mscore)
					bls.append(blength)
					
	coord2score[name] = (",".join([str(x) for x in scores]), ",".join([str(x) for x in bscores]), ",".join([str(x) for x in mscores]), ",".join([str(x) for x in bls]) )
					  
					  
					  
for interval in BedTool(args.path):
	scores, bscores, mscores, bls = coord2score.get(interval.name, ('', '', '', ''))
	interval.attrs['cons_dscores'] = scores
	interval.attrs['cons_bscores'] = bscores
	interval.attrs['cons_mscores'] = mscores
	interval.attrs['cons_bls'] = bls
	if(scores):
		sys.stdout.write(str(interval)); 
	
	
	
	

