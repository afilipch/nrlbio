#! /usr/bin/python	
'''Script produce statistics and filtering of chiparts. Chipart is a string in form: read_id, ref_id, read_start, read_end, read_length, ref_start, ref_end, ref_length, cut_left, cut_right, conversions, weight. Chipart represents found part of reference sequence(miRNA) in read sequence. Outputs right/right.fastq file which contains right parts of sequences and is used for downstream analysis. Outputs filtered chiparts which connects miRNA to the read id. If gmode is enabled tries to go three nucleotides upstream from the end of miRNA part to reach nearest G, not applicable for fulllength miRNA'''
from collections import *;
import parsing;
import chipart_lib;
import argparse;
import logging;
import os;
import sys;
import copy
import logging;

logger = logging.getLogger(__name__);
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler(os.path.join("log", "chipart.txt"), 'a')
sh = logging.StreamHandler()
logger.addHandler(fh);
logger.addHandler(sh);

wf = open(os.path.join("log", "workflow.txt"), 'a');
wf.write("python " + " ".join(sys.argv) + "\n\n");
wf.close




parser = argparse.ArgumentParser(description='Script produce statistics and filtering of chiparts');
# input files
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to chiparts file (left/long.tsv)");
parser.add_argument('-c', '--control', nargs = '+', type = str, help = "path to control chiparts file (left/long_control.tsv)");
parser.add_argument('--fastq', nargs = '?', default = "left/candidates.fastq", type = str, help = "path to the fastq file with miRNA part inside");

#control workflow
parser.add_argument('--fdr', nargs = '?', default = 0.05, type = float, help = "maximum false discovery rate allowed");

#control of the lenght of right chimeric part
parser.add_argument('-l', '--minlength', nargs = '?', default = 15, type = int, help = "minimum length of the right chipart to produce right part further")
parser.add_argument('-g', '--gmode', nargs = '?', default = False, const = True, help = "account for prefernce of T1 RNAse to cut after G in PAR-CLIP protocol");
parser.add_argument('-ms', '--mixture', nargs = '?', default = False, const = True, help = "account for a mixture of chimeras produced by exogeneous and endogeneous ligases");

args = parser.parse_args();




def gethallmark(chipart):
	'''produce hallmark of the chipart. tuple of read start and ref start supposed to be the best differntiator between signal and noise'''
	return chipart.read_start, chipart.ref_start 

	
def weighted_length(chiparts):
	'''because it is possible to get different alignments to the different references for one read, each chipart has an attribute weight which is reciprocal to the number of chiparts derived form the same read. So to get trustworthy evaluation of any features of chiparts, one should account for weights. These function outputs the sum of weights for given list of chiparts'''
	return sum([x.weight for x in chiparts])
	

	
def score_filtering(strata_real, strata_control, fdr):
	'''Apply first step of filtering based on alignment score. Outputs lists of reliable scores - that is chiparts with these scores are considered to be true positives'''
	logger.info("\nfirst step of filtering based on alignment score")
	logger.info("%s\t%s\t%s" % ("score", "real", "control"))		
	
	rel_scores = [];#chiparts which passed alignment score cutoff
	unrel_scores = [];#control chiparts which passed alignment score cutoff
	
	for score in sorted(strata_real.keys()):
		control = weighted_length(strata_control[score])
		real = weighted_length(strata_real[score])
		logger.info("%d\t%d\t%d" % (score, real, control))
		if(real*fdr > control):
			rel_scores.append(score);
		else:
			unrel_scores.append(score);	
	return rel_scores, unrel_scores;	
	
def hallmark_filtering(strata_real, strata_control, unrel_scores, hallmarks, fdr):
	'''Apply second step of filtering to reads with unreliable alignment scores based on hallmarks. Hallmarks are infered on enriched fiatures in reliable chiparts'''
	logger.info("\nsecond step of filtering based on hallmarks")
	logger.info("%s\t%s\t%s\t%s" % ("score", "hallmark", "real", "control"))	
	rescued = [];
	rescued_control = [];
	unrel = [];
	for score in unrel_scores:
		for h in hallmarks:
			real_chiparts = filter(lambda x: gethallmark(x) == h, strata_real[score])
			control_chiparts = filter(lambda x: gethallmark(x) == h, strata_control[score])
			real = weighted_length(real_chiparts);
			control = weighted_length(control_chiparts);
			if(real*fdr > control):
				rescued += real_chiparts;
				rescued_control += control_chiparts;
				logger.info("%1.1f\t%d in read %d in ref\t%d\t%d" % (score, h[0], h[1], real, control))	
			else:
				unrel += real_chiparts;
	return rescued, unrel, rescued_control
				
			

#>>>>>> getting strata dictionaries. Each stratum represent chiparts aligned with the same score. it is done for both real and control.		
chiparts = chipart_lib.read(args.path);
control = chipart_lib.read(args.control);

strata_real = defaultdict(list)
for el in chiparts:
	strata_real[el.score].append(el);
	
strata_control = defaultdict(list)
for el in control:
	strata_control[el.score].append(el)
	
rel_scores, unrel_scores = score_filtering(strata_real, strata_control, args.fdr)	

#>>>>>>>>>>>>>>>hallmarks finding and outputting as statistics
rel1 = []
for sc in rel_scores:
	rel1 += strata_real[sc];
unrel1 = []
for sc in unrel_scores:
	unrel1 += strata_real[sc];
control1 = []
for sc in rel_scores:
	control1 += strata_control[sc];
		
		
s1, s2 = 	chipart_lib.get_starts(rel1), chipart_lib.get_starts(unrel1)
chipart_lib.repr_starts(s1, s2, os.path.join("lstatistics", "start_positions.tsv"));
w = sum(s1.values())
th = [];
logger.info("\nreliable hallmarks(start positions)")
for k,v in s1.iteritems():
	fraction = v/w;
	if(fraction > 0.01):
		th.append((k,fraction));
		logger.info("%d in read %d in ref\t%1.3f fraction" % (k[0],k[1],fraction))
hallmarks = [x[0] for x in sorted(th, key = lambda x: x[1], reverse = True)]			
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
#>>>>>>>>>>>>>filtering according to the hallmarks 
rescued, unrel2, rescued_control = hallmark_filtering(strata_real, strata_control, unrel_scores, hallmarks, args.fdr)
rel2 = rescued + rel1;
control2 = control1 + rescued_control;
 
#>>>>>>>>>>>>>>other features processing and outputting as statistics
w2, w1 = weighted_length(control2), weighted_length(rel2);
myfdr = w2/w1; 
logger.info("\n%d weight of signal %d weight of control, resulted fdr %1.3f" % (w1, w2, myfdr))


stat1, stat2 = chipart_lib.get_stats(rel2), chipart_lib.get_stats(unrel2)
chipart_lib.repr_stats(stat1, stat2, os.path.join("lstatistics", "conv_numbers.tsv"), "num_conversion");
chipart_lib.repr_stats(stat1, stat2, os.path.join("lstatistics", "conv_positions.tsv"), "pos_conversion");
chipart_lib.repr_stats(stat1, stat2, os.path.join("lstatistics", "conv_types.tsv"), "type_conversion");
chipart_lib.repr_stats(stat1, stat2, os.path.join("lstatistics", "cut_between.tsv"), "cut_right");
chipart_lib.repr_stats(stat1, stat2, os.path.join("lstatistics", "cut_after.tsv"), "cut_right_after");
chipart_lib.repr_stats(stat1, stat2, os.path.join("lstatistics", "reference_end.tsv"), "ref_end");
chipart_lib.repr_stats(stat1, stat2, os.path.join("lstatistics", "right_length.tsv"), "right_length");
chipart_lib.repr_stats(stat1, stat2, os.path.join("lstatistics", "cut_length.tsv"), "cut_length");

#>>>>>>>>>>>>>>collapse chiparts derived from the same read
def collapse(chiparts, hallmarks):
	'''first takes matches with the most abundant hallmark'''
	ans = [];
	for h in hallmarks:
		for c in chiparts:
			if(h == gethallmark(c)):
				ans.append(c);
		if(ans):
			return ans;
	return copy.copy(chiparts)

	

	
def solve_ambiguous(chiparts, famid, ref_ids):	
	fam = "fam" + str(famid) + "\t" + ",".join(ref_ids);
	new_chi = chipart_lib.Chipart(chiparts[0].read_id, "fam" + str(famid), chiparts[0].read_start, chiparts[0].read_end, chiparts[0].read_length, chiparts[0].ref_start, chiparts[0].ref_end, chiparts[0].ref_length, chiparts[0].cut_left, chiparts[0].cut_right, chiparts[0].conversions, 1.0, chiparts[0].score)
	return new_chi, fam;
	
uniq = defaultdict(list);
for ch in rel2:
	uniq[ch.read_id].append(ch);

fam_handler = open(os.path.join("left", "fam_ids.tsv"), 'w')	
collapsed = [];# final list of filtered and collapsed chiparts;	
famids = {};
cf = 1;
for l in uniq.values():
	if len(l) > 1:
		a = collapse(l, hallmarks)
		if(len(set([x.ref_id for x in a]))>1):
			ref_ids = tuple(set([x.ref_id for x in a]))
			if(ref_ids in famids):
				famid = famids[ref_ids]
				ch, fam = solve_ambiguous(a, famid, ref_ids);
			else:
				famid = cf;
				famids[ref_ids] = famid;
				cf += 1;
				ch, fam = solve_ambiguous(a, famid, ref_ids);
				fam_handler.write(fam + "\n");
			collapsed.append(ch);
		else:
			collapsed.append(a[0]);
	else:
		collapsed.append(l[0]);
fam_handler.close();		
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>		
		
#>>>>>>>>output right part fastq files for filtered chiparts

def get_pos(chipart, seq, gmode):
	if(gmode):
		if(seq[chipart.read_end-1] == "G" or chipart.ref_end == chipart.ref_length):
			return chipart.read_end
		else:
			for i in range(1,4):
				if(seq[chipart.read_end-1-i] == "G"):
					return chipart.read_end - i
			return 	chipart.read_end	
	else:
		return chipart.read_end

def output(chiparts, in_fastq, right_fastq, filtered_fastq, out_chiparts, gmode):
	tot = 0
	chi_dict = {};
	for ch in chiparts:
		chi_dict[ch.read_id] = ch;
	chh = open(out_chiparts, 'w');
	rif = open(right_fastq, 'w');
	fif = open(filtered_fastq, 'w');
	inf = open(in_fastq);
	if(gmode):
		gmo = open(os.path.join("gmode", "gmode.tsv"), 'w');
	line_set = [];
	for line in inf:
		line_set.append(line.strip());
		if(len(line_set) == 4):
			c = chi_dict.get(line_set[0], None)
			if(c):
				pos = get_pos(c, line_set[1], gmode);
				if(c.read_length - pos >= args.minlength):
					tot += 1;
					rif.write("\n".join([line_set[0], line_set[1][pos:], line_set[2], line_set[3][pos:]]) + "\n");
					fif.write("\n".join(line_set) + "\n");
					chh.write(chipart_lib.repr_chipart(c)+"\n");
					if(gmode and pos != c.read_end):
						gmo.write("\t".join([c.read_id, c.ref_id, str(c.read_start), str(c.read_end), str(pos)]) + "\n")
				else:
					pass;
			else:
				pass
			line_set = [];
	chh.close();
	rif.close();
	fif.close();
	inf.close();
	if(gmode):
		gmo.close()
	return tot	
		
tot = output(collapsed, args.fastq, os.path.join("right", "right.fastq"), os.path.join("left", "filtered_candidates.fastq"), os.path.join("left", "filtered_chiparts.tsv"), args.gmode)		
logger.info("\n%d total filtered %d filtered with min length%d\n" % (w1, tot, args.minlength))

