#! /usr/bin/python	
'''Script produces interactions and output the as bed file like in previous implementations. Also it makes some statistics: cut, cut_after, and crosslink relative to the seed site'''
#import sam_lib;
import mir_lib
import chipart_lib;
import cluster_lib
import parsing;
import sys;
import os;
import argparse;
import logging;
from collections import *;

wf = open(os.path.join("log", "workflow.txt"), 'a');
wf.write("python " + " ".join(sys.argv) + "\n\n");
wf.close

#>>>> logging output is duplicated into stderr and file log/exp_chipart_length{ars.minlength}.txt 
logger = logging.getLogger(__name__);
logger.setLevel(logging.DEBUG);
fh = logging.FileHandler(os.path.join("log", "clustering.txt"));
sh = logging.StreamHandler();
logger.addHandler(fh);
logger.addHandler(sh);


parser = argparse.ArgumentParser(description='Script produces interactions from filtered uniquely mapped reads given as input');
# input files
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to filtered uniquely mapped reads(sam/mapreads.tsv)");
parser.add_argument('-c', '--chiparts', nargs = '?', required = True, type = str, help = "path to filtered left chiparts (filtered_chiparts.tsv)");
parser.add_argument('--mir', nargs = '?', required = True, type = str, help = "path to fasta file of micro RNAs (mir.fa)");
parser.add_argument('--fam', nargs = '?', required = True, type = str, help = "path to tsv file of families (fam_ids.tsv)")
parser.add_argument('-s', '--system', nargs = '?', required = True, type = str, help = "genome system")
parser.add_argument('-n', '--name', nargs = '?', default = "i", type = str, help = "name of dataset")
args = parser.parse_args();


exec("from sequence_data.systems import %s as gsys" % args.system);
extdownstream = 20;
extupstream = 20;


#>>>>>>> read mapped_reads;
mapreads = [];
f = open(args.path[0]);
for line in f:
	chromosome, genome_start, genome_end, read_id, score, strand, read_start, read_end, collapse, conversion, map_quality, cigar, match_part = line.strip().split("\t");
	if(strand == "+"):
		extleft, extright = extupstream, extdownstream
	else:
		extleft, extright = extupstream, extdownstream
		
	extrefseq = gsys.genome.get_oriented(chromosome, int(genome_start) - extleft, int(genome_end) + extright, strand).upper();
	mapreads.append(cluster_lib.Mapread(chromosome, int(genome_start), int(genome_end), read_id, int(score), strand, int(read_start), int(read_end), int(collapse), int(conversion), int(map_quality), cigar, match_part, extrefseq, extupstream, extdownstream));
f.close();

#>>>>>>> 

#for m in mapreads:
	#print m;

#>>>>>> cut statistics;
lcut = defaultdict(int);
lcut_after = defaultdict(int);
rcut = defaultdict(int);
rcut_after = defaultdict(int);
for m in mapreads:
	lcut[m.extrefseq[m.extupstream-1:m.extupstream+1]]+=1;
	lcut_after[m.extrefseq[m.extupstream-1]]+=1;
	rcut[m.extrefseq[-m.extdownstream-1: -m.extdownstream+1]]+=1;
	rcut_after[m.extrefseq[-m.extdownstream-1]]+=1;
		
	
for name, counter in [("cut_right_after.tsv", rcut_after) , ("cut_right.tsv", rcut),	("cut_left_after.tsv", lcut_after), ("cut_left.tsv", lcut)]:
	sh = open(os.path.join("rstatistics", name), 'w');
	sh.write(parsing.counter2string(counter, fraction = True))
	sh.close;
#>>>>>> 

#>>>>>>>>> generate clusters 
clusters = cluster_lib.get_clusters(mapreads);
logger.info("%d clusters generated from %d uniquely mapped reads" % (len(clusters), len(mapreads)))

#>>>>>>>> get chimeras:
chimeras = [];
mdict = dict(zip([x.read_id for x in mapreads], mapreads));
chh = open(args.chiparts);
for line in chh:
	read_id, ref_id, read_start, read_end, read_length, ref_start, ref_end, ref_length, cut_left, cut_right, conversions, weight, score = line.strip().split("\t");
	m = mdict.get(read_id, None);
	if m:
		chimeras.append(cluster_lib.Chimera(chipart_lib.Chipart(read_id, ref_id, int(read_start), int(read_end), int(read_length), int(ref_start), int(ref_end), int(ref_length), cut_left, cut_right, conversions, int(float(weight)), float(score)), m));

#>>>>>>>> get interactions
interactions = cluster_lib.get_interactions(chimeras)
logger.info("%d interactions generated from %d uniquely mapped reads" % (len(interactions), len(chimeras)))	
#>>>>>>>>>>>>>>>>>>>>>>>>>

#>>>>>>> get seed and mirna information
mirdict = parsing.fasta2dict(args.mir)
## some of mirna parts are not identifiable uniquely so for them fam ids assigned, mirnas in the same fam should have the same seed which we wnt to use:
famtemp = parsing.tsv2dict(args.fam, 0, sep = "\t");
famdict = {};
for k, v in famtemp.items():
	seeds = set()
	for mirid in v[1].split(","):
		seeds.add(mirdict[mirid][:8]);                            #!!!!!! seed 2-8 is used here;
	if(len(seeds)>1):
		logger.warning("different seeds in family %s" % k);
	famdict[k] = list(seeds)[0];
	

undef = 0;
for inter in interactions:
	undef += inter.set_mirsseq(mirdict, famdict)
logger.info("%d of %d interactions are with uncertain miRNA family member" % (undef, len(interactions)))
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#>>>>>>> evaluate extension
extensions = {}
print len(set([x.cluster.mid for x in interactions]));	
for inter in interactions:
	if(inter.cluster.strand == "+"):
		extleft, extright = extupstream, extdownstream
	else:
		extleft, extright = extupstream, extdownstream	
		
	if(inter.mirseq == "undef"):	
		extensions[inter.cluster.mid] = [extupstream, extdownstream, gsys.genome.get_oriented(inter.cluster.chromosome, int(inter.cluster.start) - extleft, int(inter.cluster.end) + extright, inter.cluster.strand).upper(), inter.longseed]
	else:	
		extensions[inter.cluster.mid] = [extupstream, extdownstream, gsys.genome.get_oriented(inter.cluster.chromosome, int(inter.cluster.start) - extleft, int(inter.cluster.end) + extright, inter.cluster.strand).upper(), inter.mirseq]


extresult = {}
#seedlike = 0;
for k, v in extensions.iteritems():	
	pos = mir_lib.mm28(v[2], v[3]);
		
	if(pos == -1):
		extu, extd = -1, -1;
	else:
		#seedlike += 1;
		extu = max(0, 14 - pos + v[0])
		extd = max(0, pos + v[1] - len(v[2]) + 7)
		
	extresult[k] = [extu, extd];
	
rchance = (1.0+3*7)/4**7 # chance to find seedlike match at random while extending to one nt more	
def eval_extension(counter, rchance, rcutoff, acutoff):
	t = sum(counter.values())
	def cutoff(sa):
		return sum(sa)/(len(sa)*rchance*t)
		
	a = [0]*21;
	s = 0;
	add = [];
	stop = 0;
	for i in range(1, 21):
		a[i] = counter[i];
		add.append(a[i]);
		rc = cutoff(a[i:i+1])
		ac = cutoff(add)
		#print i, counter[i], rc, ac
		if(rc > rcutoff and ac > acutoff):
			add = [];
			stop = i;
		else:
			pass;
		
	return stop	
	
	
uewo = eval_extension(Counter([x[0] for x in extresult.values()]), rchance, 10, 6)
uews = eval_extension(Counter([x[0] for x in extresult.values()]), rchance, 3, 2)
dewo = eval_extension(Counter([x[1] for x in extresult.values()]), rchance, 10, 6)
dews = eval_extension(Counter([x[1] for x in extresult.values()]), rchance, 3, 2)

logger.info("upstream extension was chosen as %d for targets without seed and %d as maximum for targets with seeds" % (uewo, uews));		
logger.info(parsing.counter2string(Counter([x[0] for x in extresult.values()])));	
logger.info("downstream extension was chosen as %d for targets without seed and %d as maximum for targets with seeds" % (dewo, dews));	
logger.info(parsing.counter2string(Counter([x[1] for x in extresult.values()])));	

#>>> extend and get target sequences:)

for inter in interactions:
	extu, extd = extresult[inter.cluster.mid];
	if(extu == -1 and extd == -1):
		inter.extend(uewo, dewo);
	elif(extu > uews):
		if(extd > dews):
			inter.extend(uewo, dewo);
		else:
			inter.extend(uewo, extd);
	else:
		if(extd > dews):
			inter.extend(extu, dewo);
		else:
			inter.extend(extu, extd);
	inter.set_tseq(gsys.genome.get_oriented(inter.chromosome, inter.start, inter.end, inter.strand).upper())	
	
#>>> output as bed file	
oh = open(os.path.join("output", "interactions.bed"), 'w');
#t = len(interactions);
for inter in interactions:
	if(inter.mirseq == "undef"):
		ms = inter.longseed;
		mirid = famtemp[inter.mirid][1]
	else:
		ms = inter.mirseq
		mirid = inter.mirid
	
	arr = inter.chromosome, inter.start, inter.end, '%s%05d' % (args.name, inter.cluster.mid), inter.cluster.score, inter.strand, mirid, ms, inter.tseq, inter.cluster.map_quality, inter.cluster.totalreads, inter.cluster.indreads, inter.cluster.indseqs 
	oh.write("\t".join([str(x) for x in arr]) + "\n")
oh.close()	

#>>>output as file to connect interactions with read ids
oh = open(os.path.join("output", "read2int.tsv"), 'w');
for inter in interactions:
	for m in inter.cluster.mapreads:
		if(m.conversion > -1):
			c = m.genome_start - inter.start + m.conversion
		else:
			c = -1;
		oh.write('%s\t%s%05d\t%d\n' % (m.read_id, args.name, inter.cluster.mid, c))
oh.close()	

# find position relative to crosslink:
crosslink_relative_seed = defaultdict(int);
tcontrol = defaultdict(int)
for inter in interactions:
	pos = inter.tseq.rfind(inter.seed);
	if(pos > -1):
		for m in inter.cluster.mapreads:
			if(m.conversion > -1):
				if(inter.strand == "+"):
					rel = m.genome_start - inter.start + m.conversion - pos
				if(inter.strand == "-"):
					rel = inter.end - m.genome_end + m.conversion - pos
				crosslink_relative_seed[rel] += 1;
		for i in range(len(inter.tseq)):
			if(inter.tseq[i] == "T" and i >= inter.extu and i < len(inter.tseq) - inter.extd):
				tcontrol[i - pos] += len(inter.cluster.mapreads);

if(crosslink_relative_seed):				
	n = []
	for i in range(-15, 11):
		if(tcontrol[i] + crosslink_relative_seed[i] > 0):
			n.append(crosslink_relative_seed[i]/(tcontrol[i]+0.1));
	norma = len(n)/(sum(n)+1)
	
	for i in range(-15, 11):
		if(tcontrol[i] + crosslink_relative_seed[i] == 0):
			tcontrol[i] = norma;
			crosslink_relative_seed[i] = 1		
	crosslink_normalized = defaultdict(float);
	for i in range(-15, 11):
		crosslink_normalized[i] = norma*crosslink_relative_seed[i]/(0.1+tcontrol[i])

	ch = open(os.path.join("rstatistics", "crosslink.tsv"), 'w');
	ch.write(parsing.counter2string(crosslink_normalized));
	ch.close()
				


		


