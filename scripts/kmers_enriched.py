#! /usr/bin/python
import sys;
import os;
import argparse;
from collections import defaultdict, Counter;

from multiprocessing import Pool

from nrlbio.generators import generator_seq;

parser = argparse.ArgumentParser(description='search for highly enriched mers in the reads, outputs results to  STDOUT');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to fasta/fastq file of reads");
parser.add_argument('--ftype', nargs = '?', default = 'fasta', choices = ['fasta', 'fastq'], type = str, help = "format of input file. Can be fasta or fastq");
parser.add_argument('-l', '--length', nargs = '?', default = 10, type = int, help = "length of mer");
parser.add_argument('-m', '--mode', nargs = '?', default = "all", choices = ['all', 'start', 'end'], type = str, help = "search for all kmers, from start, from end; can be all|start|end");
parser.add_argument('-n', '--steps', nargs = '?', default = 80, const = 1, type = int, help = "number of steps from start or end position");
parser.add_argument('--sparse_coefficient', nargs = '?', default = 1, type = int, help = "If sparse_coefficient set to 3, then only each 3rd record will be analyzed, saves time in a case of large files")
args = parser.parse_args();




def count_start(seq, length, steps):
	return set([seq[x:x+length] for x in range(0, min(steps, len(seq) - length + 1) )]);
	
def count_end(seq, length, steps):
	stop = len(seq) - length + 1
	start = max(0, stop-steps)
	return set([seq[x:x+length] for x in range(start, stop)]);
	
def count_all(seq, length, steps):
	return set([seq[x:x+length] for x in range(0, len(seq) - length + 1)]);
	
length = args.length
steps = args.steps	
sparse = args.sparse_coefficient;
	
if(args.mode == "start"):
	countfun = count_start
elif(args.mode == "end"):
	countfun = count_end
else:
	countfun = count_all
	
	
	


counts = defaultdict(int);
if(sparse==1):
	for seq in generator_seq(args.path, args.ftype):
		for el in countfun(seq, length, steps):
			counts[el] += 1
else:
	for i, seq in enumerate(generator_seq(args.path, args.ftype)):
		if(i % sparse == 0):
			for el in countfun(seq, length, steps):
				counts[el] += 1	
  

  

for mer in sorted(counts.keys(), key = lambda x: counts[x], reverse = True)[:30]:
	print "%s\t%d" % (mer, counts[mer])	



	

	
#pool = Pool(6)
#res = pool.imap(countfun, generator_seq(args.path, args.ftype), chunksize=140)

#for m in res:
	#pass;
	
	
#sys.exit();



 
