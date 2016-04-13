#! /usr/bin/python
'''Finds miRNAs consereved in other species along with their respective counterparts'''
import argparse
import sys;
from collections import defaultdict

from Bio import SeqIO, pairwise2


parser = argparse.ArgumentParser(description='Finds miRNAs consereved in other species along with their respective counterparts');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the MirBase miRNA mature sequences, fasta format");
parser.add_argument('--specie', nargs = '?', required = True, type = str, help = "Name of the reference specie (MirBase coding: cel, hsa, mmu and so on)")
parser.add_argument('--minscore', nargs = '?', default = 45, type = float, help = "Min conservation score required to count miRNAs as related")
parser.add_argument('--table', nargs = '?', default = '', type = str, help = "Path to the MirBase:UCSC conservation translation table. If set, conservation will be computed only for those species in MirBase which have counterparts in UCSC multi-alignments")
args = parser.parse_args();

select_only = set()
if(args.table):
	with open(args.table) as f:
		for l in f:
			select_only.add(l.strip().split("\t")[1])



START = 1
END = 8



#Get a miRNA conservation score for the query and reference miRNAs
def get_score(rseq, qseq, rid, qid):
	scores = [0,];
	
	for alignment in pairwise2.align.localms(rseq, qseq, 2, -3, -5, -2):
		score = alignment[2]
		if( tuple(rid.split("-")[1:3]) ==  tuple(qid.split("-")[1:3]) ):
			score += 5
		
		for r, q in zip(rseq[START:END], qseq[START:END]):
			if(r==q):
				score += 2
			else:
				score -= 2
		scores.append(score)
		
	return max(scores)


#Get the most similar miRNA to the reference one from a particular specie
def get_best_alignment(rseq, rid, qmirnas):
	selected = [];
	for qid, qseq in qmirnas.items():
		score = get_score(rseq, qseq, rid, qid)
		selected.append((qid, qseq, score))
	return max(selected, key = lambda x: x[2])



#Complile mature miRNAs into dictionaries accroding to the specie name and miRNA name
specie2mir = defaultdict(dict);
refspecie = {};

for seqrecord in SeqIO.parse(args.path, 'fasta'):
	specie = seqrecord.name.split("-")[0]
	if(specie == args.specie):
		refspecie[seqrecord.id] = str(seqrecord.seq.upper()).replace('U', 'T');
	elif(not select_only or specie in select_only):
		specie2mir[specie][seqrecord.id] = str(seqrecord.seq.upper()).replace('U', 'T');
		
		
#Run the analysis for all miRNA
for count, (rmirid, rseq) in enumerate(refspecie.iteritems()):
	print ">%s\n%s" % (rmirid, rseq)
	for qmirnas in specie2mir.values():
		bmirid, bseq, bscore = get_best_alignment(rseq, rmirid, qmirnas)
		if(bscore>=args.minscore):
			print ">%s %1.1f\n%s" % (bmirid, bscore, bseq)
	if((count+1) % 10 == 0):
		sys.stderr.write("%d miRNAs processed\n" % (count+1));
		#sys.exit()
