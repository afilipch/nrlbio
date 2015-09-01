#! /usr/bin/python	
'''Script tries to predict miRNA bound to provided targets'''
#import sam_lib;
import rnahybrid;
import argparse;
from Bio import SeqIO;



parser = argparse.ArgumentParser(description='Script tries to predict miRNA bound to provided targets');
# input files
parser.add_argument('--targets', nargs = '?', type = str, help = "potential targets as fasta file");
parser.add_argument('--mir', nargs = '?', type = str, help = "expressed (or all) miRNA for organism of interest as fasta file");
parser.add_argument('--seed', nargs = '?', default = False, const = True, help = "do require seed 3-7");
parser.add_argument('-n', '--best', nargs = '?', default = 10, type = int, help = "show n best miRNA");
args = parser.parse_args();

if(args.seed):
	seed = "-f 3,7 "
else:
	seed = ""

def find_best(target, mirdict, n, seed):
	result = [];
	for k, v in mirdict.items():
		result.append([k] + list(rnahybrid.represent(str(target.seq), str(v.seq), arguments = seed)));
	result.sort(key = lambda x: x[2])
	return ["miRNA: %s\n%s" % (x[0],x[1]) for x in result[:n]]
	
	


mirdict = SeqIO.to_dict(SeqIO.parse(args.mir, "fasta"))


for target in SeqIO.parse(args.targets, "fasta"):
	print target.description + "\n"
	for s in find_best(target, mirdict, args.best, seed):
		print "-"*100
		print s;

	print "\n"*5	

	
#for k, v in mirdict.items():
	#print k, dir(v.seq)