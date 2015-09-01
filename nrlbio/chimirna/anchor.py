#! /usr/bin/python	
'''Script searches for reference(miRNA sequences) anchors in reads from fastq files and extends them using alignments. Outputs all matches in a format read_id, miRNA id, start of match in read(0based), end of match in reads(not included in match), length of read, start of match in reference(0based), end of match in reference(not included in match), length of reference, cut left, cut right, conversions in form (position relative to reference, nucl/gap in reference, nucl/gep in read) "|" separated if more than 1, weight of match(1 divided to the number of matches fo this read), alignment score. This script can work in control mode, shuffling sequences in situ. Also outputs candidates.fastq file to narrow downstream loading of sequences containing chimeras'''


import alignments;
import pieces;
import generators;
import parsing;
import sys;
import re;
import os;
import copy;
import argparse;
import multiprocessing;
from collections import Counter, defaultdict;
import itertools;
import logging;
from time import gmtime, strftime;

wf = open(os.path.join("log", "workflow.txt"), 'a');
wf.write("python " + " ".join(sys.argv) + "\n\n");
wf.close


logger = logging.getLogger(__name__);
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler(os.path.join("log", "anchors.txt"), 'a')
sh = logging.StreamHandler()
logger.addHandler(fh);
logger.addHandler(sh);



parser = argparse.ArgumentParser(description='Script searches for reference(miRNA sequences) anchors in reads from fastq files and extends them using alignments. Outputs all matches in a format read_id, miRNA id, start of match in read(0based), end of match in reads(not included in match), length of read, start of match in reference(0based), end of match in reference(not included in match), length of reference, cut left, cut right, conversions in form (position relative to reference, nucl/gap in reference, nucl/gep in read) "|" separated if more than 1, weight of match(1 divided to the number of matches fo this read), alignment score. This script can work in control mode, shuffling sequences in situ');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to fastq file");
parser.add_argument('-r', '--ref', nargs = '?',  type = str, required = True, help = "path to reference sequence - miRNAs in our case");
parser.add_argument('--control', nargs = '?',  default = False, const = True, type = bool, help = "swith to control mode - shuffle sequences from provided fastq file in situ");

## options controlling division into pieces
parser.add_argument('-l', '--length', nargs = '?', default = 12, type = int, help = "length of pieces");
parser.add_argument('-fr', '--first_ref', nargs = '?', default = 3, type = int, help = "use only first n pieces counting from the beginning of the reference sequences");
parser.add_argument('-fs', '--first_seq', nargs = '?', default = 3, type = int, help = "use only first n pieces counting from the beginning of the read sequences");
parser.add_argument('-m', '--mode', nargs = '?', default = "exp", choices = ["exp", "hc", "pc"], type = str, help = "this mode controls type of mutations allowed in pieces");

## options controlling alignment
parser.add_argument('--match', nargs = '?', default = 2, type = int,  help = "award for match");
parser.add_argument('--missmatch', nargs = '?', default = -5, type = int,  help = "penalty for missmatch");
parser.add_argument('--opening', nargs = '?', default = -6, type = int,  help = "penalty for indel opening");
parser.add_argument('--extension', nargs = '?', default = -4, type = int,  help = "penalty for indel extension");

## options controlling output
parser.add_argument('--fastq', nargs = '?', default = "left/candidates", type = str, help = "path to the output fastq file with miRNA part inside");
parser.add_argument('-o', '--output', nargs = '?', default = "left/anchors", type = str, help = "name of the output file without extension")

## other options 
parser.add_argument('-t', '--threads', nargs = '?', default = 8, type = int,  help = "number of threads");

args = parser.parse_args();

## get anchors(pieces)
ref_dict = parsing.fasta2dict(args.ref)
c_mode = pieces.get_mode(args.mode)
pieces = pieces.get_pieces(ref_dict, args.first_ref, args.length, mode = c_mode)

def get_hits(seq, pieces, first, length):
	'''outputs which reference items (ids) have piece inside the part of read. Which part is controlled by args.length and args.first_ref'''
	hits = set();
	part = seq[:first + length];
	for k, v in pieces.iteritems():
		if(k in part):
			hits.update(v)
	return hits;	
	
def align(hits, ref_dict, sid, seq, m, mm, g, e):
	'''outputs all best alignments to all anchoring reference sequences on base of hits found by get_hits function'''
	refseq_dict = {};
	for el in hits:
		refseq_dict[el] = ref_dict[el];
	ans = [];	
	for a in alignments.get_best_alignments(refseq_dict, sid, seq, m, mm, g, e):
		ans.append(alignments.analyse_alignment(a, True))
	return ans	
	
def execute(arr):
	'''finds anchors and further extend them to best alignments for a chunk (list of line_sets) of reads'''
	seq_list, pieces, first, length, ref_dict, m, mm, g, e = arr;
	ans = [];
	fastq = [];
	
	for line_set in seq_list:
		hits = get_hits(line_set[1], pieces, first, length)
		if(hits):
			fastq += line_set;			
			for string in align(hits, ref_dict, line_set[0], line_set[1], m, mm, g, e):
				ans.append(string)			
	return ans, fastq	
	
chunk = 10000	
pool = multiprocessing.Pool(processes = 8)
res = pool.imap(execute, ((seq_list, pieces, args.first_seq, args.length, ref_dict, args.match, args.missmatch, args.opening, args.extension) for seq_list in generators.grouper(generators.generator_fastq(args.path, ["id", "seq", "sign", "qual"], shuffle = args.control), chunk)) )   


##>>> outputing
if(args.control):
	ff = args.fastq + "_control.fastq"
	af = args.output + "_control.tsv"
else:
	ff = args.fastq + ".fastq"
	af = args.output + ".tsv"	

f_handler = open(ff, 'w')
a_handler = open(af, 'w')
c = 0;
a = 0;
found = 0;
for arr, fastq in res:
	c+= chunk;
	sys.stderr.write("%d\n" % c);
	if arr:
		a += len(arr);
		a_handler.write("\n".join(arr)+"\n");
	if(fastq):
		found += len(fastq)/4;
		f_handler.write("\n".join(fastq) + "\n");
logger.info(strftime("time of analysis %Y-%m-%d %H:%M:%S", gmtime()))		
logger.info("in %d reads anchors were found in %d total reads\n" % (found, c))
logger.info("%d anchors were found in %d total reads\n" % (a, c))
f_handler.close()	
a_handler.close()







