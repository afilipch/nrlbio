#! /usr/bin/python	
'''Script filter out all sam records for the same read with alignment score less or equal to maximum. Output is bed-like file'''
import sam_lib;
#import miscellaneous
import parsing;
import sys;
import os;
import argparse;
import logging;
from collections import *;
import itertools;

wf = open(os.path.join("log", "workflow.txt"), 'a');
wf.write("python " + " ".join(sys.argv) + "\n\n");
wf.close


parser = argparse.ArgumentParser(description='Script filter out all sam records for the same read with alignment score less or equal to maximum. Output is bed-like file');
# input files
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to sam file(result of mapping true reads)");
parser.add_argument('-c', '--control', nargs = '?', type = str, help = "path to sam file(result of mapping control reads)");
parser.add_argument('-r', '--ref', nargs = '?', required = True, type = str, help = "path to the bed file which was the reference for mapping")
parser.add_argument('-o', '--output', nargs = '?', default = "sam/mapreads.tsv", type = str, help = "name of the output file without extension")
parser.add_argument('-s', '--system', nargs = '?', required = True, type = str, help = "genome system")
parser.add_argument('--fdr', nargs = '?', default = 0.05, type = float, help = "maximum false discovery rate allowed");
args = parser.parse_args();

#>>>> logging output is duplicated into stderr and file log/exp_chipart_length{ars.minlength}.txt 
logger = logging.getLogger(__name__);
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler(os.path.join("log", "mapreads.txt"))
sh = logging.StreamHandler()
logger.addHandler(fh);
logger.addHandler(sh);


def collapse_same_read(same_read):
	''' provides only those mappings for the same read which have the best alignment score. suboptimal in terms that can report mappings to the same locus, but in a way which will not violate uniqness of mapping'''
	if(not same_read):
		return None;		
	maximum = max([int(x[11][5:]) for x in same_read])
	top = [];
	start = -1000; 
	refname = "undef";	
	for arr in same_read:
		if(int(arr[11][5:]) == maximum and (arr[2] != refname or abs(int(arr[3]) - start) > 15)):
			top.append(arr);
			refname = arr[2]
			start = int(arr[3]);				
	return top;
	
		
def arr2sam_record(arr):
	'''parsing of sam record into an object of Sam_record class'''
	n = arr[0].split("_");
	read_id = "@" + "_".join(n[:-2]);
	collapse = int(n[-3][1:])
	conversion = int(n[-1])
	ref_id = arr[2];
	start = int(arr[3]) - 1 
	map_quality = int(arr[4])
	cigar = arr[5]
	sequence = arr[9]
	score = int(arr[11][5:])
	#print arr;
	#print read_id, ref_id, collapse, conversion, start, map_quality, cigar, sequence, score
	return sam_lib.Sam_record(read_id, ref_id, collapse, conversion, start, map_quality, cigar, sequence, score)
	
def get_sam_records(path, prefix):
	'''converts sam file into list of Sam_record objects, based on uniquely mapped reads'''
	trash = [];
	sam_records = []		
	sh = open(path);
	same_read = [];
	readid = "";
	t, u = -1,0
	
	##>> all reads processing
	for line in sh:
		arr = line.split("\t")
		tid = "_".join(arr[0].split("_")[:-3]);
		if(tid == readid):
			same_read.append(arr)
		else:
			t += 1;
			carr = collapse_same_read(same_read);
			same_read = [arr];
			readid = tid;
			if(carr and len(carr) == 1):
				u+=1;
				sam_records.append(arr2sam_record(carr[0]))
			elif(carr):
				for el in carr:
					trash.append(arr2sam_record(el))
	
	###>> last read processing
	t += 1;
	carr = collapse_same_read(same_read);	
	if(carr and len(carr) == 1):
		u+=1;
		sam_records.append(arr2sam_record(carr[0]))
	elif(carr):
		for el in carr:
			trash.append(arr2sam_record(el))		
		
	sh.close()
	logger.info("in %s %d total reads mapped %d mapped uniquely" % (prefix, t, u))
	
	#salvage = sam_lib.intersect(trash, sam_records);
	return sam_records;
	
def score_filtering(real, control, fdr):
	'''apply first step of filtering based on alignment score, Outputs cutoff, sam_records with score more or equal than that are considered to be reliable. Also searches for max reliable gap between 
	left and right chiparts'''
	logger.info("\nfirst step of filtering based on alignment score")
	logger.info("%s\t%s\t%s\t%s" % ("score", "real", "control", "fdr"))	
	filtered = [];
	cutoff = 1000;
	cont = True;
	real_scores = [x.score for x in real]
	control_scores = [x.score for x in control]
	
	for score in sorted(set(real_scores + control_scores)):
		rn = len(filter(lambda x: x >= score, real_scores))
		cn = len(filter(lambda x: x >= score, control_scores))
		myfdr = cn/float(rn)
		logger.info("%d\t%d\t%d\t%1.3f" % (score, rn, cn, myfdr))
		
		if(cont and fdr > myfdr):
			filtered = filter(lambda x: x.score >= score, real)
			cutoff = score;
			cont = False
		else:
			pass;
	
	#search for max_gap: gap between left and right part of chimera, here it equals to read.start
	logger.info("\nsearch for optimal max gap requirments for reads")
	logger.info("%s\t%s\t%s\t%s" % ("gap", "real", "total", "fraction"))			
	max_gap = 0 
	cont = True;
	total = float(len(filtered));
	for i in range(20):
		n = len(filter(lambda x: x.read_start <= i, filtered))
		fraction = n/total;
		logger.info("%d\t%d\t%d\t%1.3f" % (i, n, total, fraction))
		if(cont and fraction > 1 - 2*fdr):
			max_gap = i;
			cont = False;
		else:
			pass;
					
	return cutoff, max_gap;	
	
def gap_filtering(real, control, max_gap, fdr):
	'''try to salvage reads with low alignment score based on gap filtering'''
	logger.info("\nsecond step of filtering based on max gap %d" % max_gap)
	logger.info("%s\t%s\t%s\t%s" % ("score", "real", "control", "fdr"))	
	soft_cutoff = 1000;
	cont = True;
	
	for score in sorted(set([x.score for x in real + control])):
		rn = len(filter(lambda x: x.score >= score and x.read_start <= max_gap, real))
		cn = len(filter(lambda x: x.score >= score and x.read_start <= max_gap, control))
		myfdr = cn/float(rn)
		logger.info("%d\t%d\t%d\t%1.3f" % (score, rn, cn, myfdr))
		
		if(cont and fdr > myfdr):
			soft_cutoff = score;
			cont = False
		else:
			pass;	
	return soft_cutoff;
	
	
sr = get_sam_records(args.path[0], "signal")	
sr_control = get_sam_records(args.control, "control")	
cutoff, max_gap = score_filtering(sr, sr_control, args.fdr)

rel1 = filter(lambda x: x.score >= cutoff, sr)
unrel1 = filter(lambda x: x.score < cutoff, sr)
rel_control1 = filter(lambda x: x.score >= cutoff, sr_control)
unrel_control1 = filter(lambda x: x.score < cutoff, sr_control)


soft_cutoff = gap_filtering(unrel1, unrel_control1, max_gap, args.fdr*3);
rel2 = filter(lambda x:  soft_cutoff <= x.score < cutoff and x.read_start <= max_gap, sr)	
rel_control2 = filter(lambda x: soft_cutoff <= x.score < cutoff and x.read_start <= max_gap, sr_control)

s1 = rel1 + rel2;
sc1 = rel_control1 + rel_control2;
s2 = filter(lambda x:  soft_cutoff <= x.score and x.read_start <= max_gap, sr)	
sc2 = filter(lambda x: soft_cutoff <= x.score and x.read_start <= max_gap, sr_control)
w1 = len(s1);
w2 = len(s2);
if(w2*1.2 > w1 and w2/(len(sc2) +0.001) >= 20):
	reliable = s2;
	reliable_control = sc2;
	unreliable = filter(lambda x:  soft_cutoff > x.score or x.read_start > max_gap, sr)
else:
	reliable = s1;
	reliable_control = sc1;	
	unreliable = unrel1 + filter(lambda x:  not(soft_cutoff <= x.score < cutoff and x.read_start <= max_gap), sr)	

	
	
w1, w2 = len(reliable), len(reliable_control);
myfdr = w2/float(w1); 
logger.info("\n%d weight of signal %d weight of control, resulted fdr %1.3f" % (w1, w2, myfdr))

## salvage the reads which are unreliable but intersect reliable;

salvage = sam_lib.intersect(unreliable, reliable);
logger.info("\n%d of unreliable mappings were salvaged, because of intersection with reliable" % len(salvage))
reliable += salvage;
#print len(set([x.read_id for x in reliable]));
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> assign to genome
system = args.system
exec("from sequence_data.systems import %s as gsys" % system);
def assign_genome_coord(srlist, arr, system):
	for sr in srlist:
		if(arr[5] == "+"):
			sr.genome_start = int(arr[1])+ sr.ref_start;
			sr.genome_end = int(arr[1])+ sr.ref_end;
		else:
			sr.genome_start = int(arr[2]) - sr.ref_end;
			sr.genome_end =  int(arr[2]) - sr.ref_start; 		
		sr.chromosome = arr[0];
		sr.strand = arr[5];
		sr.refseq = gsys.genome.get_oriented(sr.chromosome, sr.genome_start, sr.genome_end, sr.strand).upper();
		#print sr.chromosome, sr.genome_start, sr.genome_end, sr.read_id, sr.score, sr.strand, sr.refseq, sr.match_part, sr.cigar
		
srdict = defaultdict(list)
for s in reliable:
	srdict[s.ref_id].append(s);
	
srdict_control = defaultdict(list)	
for s in reliable_control:
	srdict_control[s.ref_id].append(s);
	
r = open(args.ref);
for line in r:
	arr = line.strip().split("\t")
	srlist = srdict[arr[3]]
	if(srlist):
		assign_genome_coord(srlist, arr, system)	
	srlist_control = srdict_control[arr[3]]
	if(srlist_control):
		assign_genome_coord(srlist_control, arr, system)
r.close();			
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> get statistics: conversions, num of conversions 
def conversions(srlist):
	types = defaultdict(float);
	for n in itertools.permutations("ACTG-",2):
		types["".join(n)] = 0.1
	num = defaultdict(int);
	for sr in srlist:
		a = sr.get_conversions();
		num[len(a)] += 1;
		for e in a:
			types[e[1]+e[2]] += 1.0/len(a);
		if(not a):
			types['no'] += 1;
	return types, num;
	
types, num = conversions(reliable);
types_control, num_control = conversions(reliable_control);

sh = open(os.path.join("rstatistics", "conv_numbers.tsv"), 'w');
sh.write(parsing.counter2string(num))
sh.close();
sh = open(os.path.join("rstatistics", "conv_types.tsv"), 'w');
sh.write(parsing.counters2string(types, types_control))
sh.close();

#per read conversion
sh = open(os.path.join("rstatistics", "read2conv.tsv"), 'w');
for sr in reliable:
	a = sr.get_conversions();
	sc = "|".join([x[1]+x[2] for x in a])
	if(not sc):
		sc = 'no'
	sh.write("%s\t%s\n" % (sr.read_id, sc))
sh.close();

		
		
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> outputing
o = open(args.output, 'w')
for sr in reliable:
	sarr = sr.chromosome, sr.genome_start, sr.genome_end, sr.read_id, sr.score, sr.strand, sr.read_start, sr.read_end, sr.collapse, sr.conversion, sr.map_quality, sr.cigar, sr.match_part;
	o.write("\t".join([str(x) for x in sarr]) + "\n");
o.close();	

o = open(os.path.join('sam', 'right.tsv'), 'w')
for sr in reliable:
	sarr = sr.chromosome, sr.genome_start, sr.genome_end, sr.read_id, sr.score, sr.strand, sr.read_start, sr.read_end, sr.collapse, sr.conversion, sr.get_alignment()[0]
	o.write("\t".join([str(x) for x in sarr]) + "\n");
o.close();	




