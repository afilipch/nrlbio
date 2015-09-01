import parsing;
import sys;
import os;
import argparse;
from collections import Counter, defaultdict;
import logging;

logger = logging.getLogger(__name__);
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler(os.path.join("log", "gmode_evalute.txt"))
sh = logging.StreamHandler()
logger.addHandler(fh);
logger.addHandler(sh);


parser = argparse.ArgumentParser(description='exploring productivity of gmode, assignes new left chipart.read_end values');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to gmode file which comes from left chipart filtering");
parser.add_argument('-b', '--bed', nargs = '?', type = str, required = True, help = "path file of filtered right chiparts(sam/mapreads.tsv)");
parser.add_argument('-l', '--left', nargs = '?', type = str, required = True, help = "path file of filtered left chiparts(left/filtered_chiparts.tsv)");
parser.add_argument('--mir', nargs = '?', required = True, type = str, help = "path to fasta file of micro RNAs (mir.fa)");
parser.add_argument('-o', '--output', nargs = '?', default = 'gmode/chiparts.tsv', type = str, help = "name of the output file")
args = parser.parse_args();

def evaluate(garr, barr):
	gadd = int(garr[3]) - int(garr[4])
	gmap = int(barr[6])
	if(not gmap):
		return "worthwhile", gadd
	elif(gmap == gadd):
		return "worthless", 0
	else:
		return "unclear", 0
	
stat = defaultdict(int)
gmode = parsing.tsv2dict(args.path[0], 0, sep = "\t");

negate = {};

f = open(args.bed)
for line in f:
	arr = line.strip().split("\t");
	garr = gmode.get(arr[3], None)
	if(garr):
		label, neg = evaluate(garr ,arr)
		stat[label] += 1;
		negate[arr[3]] = neg;
f.close()

for k in ["worthwhile", "unclear", "worthless"]:
	logger.info("%s\t%d" % (k, stat[k]))
	
	
mirdict = parsing.fasta2dict(args.mir)	
cut = defaultdict(int)	
cut_after = defaultdict(int)
cut_length = defaultdict(int);
mir_length = defaultdict(int);

f = open(args.left)
o = open(args.output, 'w');
for line in f:
	arr = line.strip().split("\t");
	neg = negate.get(arr[0], "none")
	if(neg != "none"):
		arr[6] = str(int(arr[6]) - negate[arr[0]])
		mirseq = mirdict.get(arr[1], "n");
		if(mirseq != 'n'):
			arr[9] = mirseq[int(arr[6])-1:int(arr[6])+1]
	cut[arr[9]] += 1;
	cut_after[arr[9][0]] += 1; 
	cut_length[int(arr[7]) - int(arr[6])] += 1;
	mir_length[arr[6]] += 1;
	o.write("\t".join(arr) + "\n");	
o.close()	
f.close()	
	
	
o = open('gmode/cut.tsv', 'w');
keys = filter(lambda x: x != "NN", cut)
t = sum([cut[x] for x in keys])	
for k in sorted(keys):
	o.write("%s\t%1.4f\n" % (k, float(cut[k])/t)); 
o.close()
	
o = open('gmode/cut_after.tsv', 'w');
keys = filter(lambda x: x != "N", cut_after)
t = sum([cut_after[x] for x in keys])	
for k in sorted(keys):
	o.write("%s\t%1.4f\n" % (k, float(cut_after[k])/t)); 
o.close()	
	
o = open('gmode/cut_length.tsv', 'w');
o.write(parsing.counter2string(cut_length, fraction = True)); 
o.close()	
	
o = open('gmode/mir_length.tsv', 'w');
o.write(parsing.counter2string(mir_length, fraction = True)); 
o.close()		
	
