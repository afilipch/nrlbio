import parsing;
import sys;
import os;
import argparse;
import chipart_lib;
import mir_lib;
import interaction_lib;
import gconst;
from collections import *;

parser = argparse.ArgumentParser(description='reanalize partition of the reads into miRNA and target parts');
parser.add_argument('-l', '--left', nargs = '+', type = str, required = True, help = "path to left chipart file (adjusted_chiparts.tsv)");
parser.add_argument('-i', '--interactions', nargs = '?', type = str, help = "path to interactions (output/interactions.bed)");
parser.add_argument('-r', '--reads', nargs = '?', type = str, required = True, help = "path to reads2int (output/reads2int.tsv)");
parser.add_argument('-mr', '--mir', nargs = '?', type = str, required = True, help = "path mirna in fasta format");
parser.add_argument('-f', '--fams', nargs = '?', type = str, required = True, help = "path file of families(left/fam_ids.tsv)");
# i have to decide further if i want to hardly link system to miRNA library
parser.add_argument('-s', '--system', nargs = '?', required = False, type = str, help = "genome system")
args = parser.parse_args();

def reassign(chipart, mir_dict, fams2mir):
	'''assignes to each left read possible hits to miRNAs according to new juncrtion annotation'''
	hits = [];
	if(chipart.ref_id in fams2mir):
		s = mir_dict[list(fams2mir[chipart.ref_id])[0]][chipart.ref_start:chipart.ref_end]
	elif(chipart.ref_id in mir_dict): 
		s = mir_dict[chipart.ref_id][chipart.ref_start:chipart.ref_end]
	
	for k, v in mir_dict.iteritems():
		ss = v[:chipart.ref_end];
		pos = ss.find(s);
		if (pos > -1):
			hits.append((k, pos));
	#if('mmu-miR-466f-3p' in [x[0] for x in hits] and 'mmu-miR-466i-3p' in [x[0] for x in hits]):
		#print s;
		#print hits;
		#m = min([x[1] for x in hits])
		#print m
		#print [x[0] for x in hits if x[1] == m]	
		#print
	m = min([x[1] for x in hits])
			
	return [x[0] for x in hits if x[1] == m]		

#print args.left
mir_dict = parsing.fasta2dict(args.mir);
int2reads = interaction_lib.int2reads(args.reads);
chiparts = chipart_lib.read(args.left);
fams2mir = mir_lib.read_families(args.fams);
interactions = interaction_lib.get(args.interactions, undef = True);

new_chiparts = [];
c = 1;
#k = 1;

for chipart in chiparts:
	if(chipart.ref_end - chipart.ref_start <11):
		continue;
	
	hits = reassign(chipart, mir_dict, fams2mir);
	#if('mmu-miR-466f-3p' in hits and 'mmu-miR-466i-3p' in hits):
		#print "bu"
	
	if(len(hits)>1):
		#print len(fams2mir)
		m = set(hits);
		for fam, mirs in fams2mir.items():
			if (m == mirs):
				#k+=1;
				nc = chipart_lib.Chipart(chipart.read_id, fam, chipart.read_start, chipart.read_end, chipart.read_length, chipart.ref_start, chipart.ref_end, chipart.ref_length, chipart.cut_left, chipart.cut_right, chipart.conversions, chipart.weight, chipart.score)
				break;
		else:
			fam = "fam%dra" % c;
			#print m;
			c += 1;
			nc = chipart_lib.Chipart(chipart.read_id, fam, chipart.read_start, chipart.read_end, chipart.read_length, chipart.ref_start, chipart.ref_end, chipart.ref_length, chipart.cut_left, chipart.cut_right, chipart.conversions, chipart.weight, chipart.score);
			fams2mir[fam] = m;
	else:
			nc = chipart_lib.Chipart(chipart.read_id, hits[0], chipart.read_start, chipart.read_end, chipart.read_length, chipart.ref_start, chipart.ref_end, chipart.ref_length, chipart.cut_left, chipart.cut_right, chipart.conversions, chipart.weight, chipart.score)
		
	new_chiparts.append(nc)
	
print c;
print len(new_chiparts), len(chiparts);

read2chipart = dict([(c.read_id, c) for c in new_chiparts])


#now we want to change the interactions.bed: reassign families according new chiparts. If an interactions consists of more than 1 discorcodant chiparts we SHOULD split it
new_interactions = [];
new_int2reads = {};
for inter in interactions:
	reads = int2reads[inter.iid];
	#print reads;
	mirids = set();
	mirid2reads = defaultdict(list)
	for read in reads:
		if(read in read2chipart):
			mirids.add(read2chipart[read].ref_id);
			mirid2reads[read2chipart[read].ref_id].append(read)
			
	for c, mirid in enumerate(mirids):
		if(mirid in fams2mir):
			sid = ",".join(fams2mir[mirid])
			#mirseq = mir_dict[list(fams2mir[mirid])[0]][1:8];
			common, pos = mir_lib.get_common_part([mir_dict[x] for x in fams2mir[mirid]], fr = 0, tolerate_first_mismatch = 4)
			mirseq = common + "XXX"
			#if(not pos):
				#print common, pos
				#for el in [mir_dict[x] for x in fams2mir[mirid]]:
					#print el;
		else:
			sid = mirid
			mirseq = inter.mirseq
		if(len(mirids) > 1):
			iid = "%s-%d" % (inter.iid, c+1);
		else:
			iid = inter.iid
			
		new_interactions.append(interaction_lib.Interaction(inter.chromosome, inter.start, inter.end, iid, inter.score, inter.strand, sid, mirseq, inter.tseq, inter.map_quality, inter.totalreads, inter.indreads, inter.indseqs))	
		new_int2reads[iid] = mirid2reads[mirid]
		
print len(new_interactions)

directory = "final"
if not os.path.exists(directory):
    os.makedirs(directory);
    
## output into "final" directory

w = open(os.path.join(directory, "interactions.bed"), 'w');
for inter in new_interactions:
	w.write("\t".join([str(x) for x in inter]) + "\n")
w.close()	
    
    
w = open(os.path.join(directory, "read2int.tsv"), 'w');
for k, v in new_int2reads.items():
	for r in v:
		w.write("%s\t%s\n" % (r,k));   
w.close();		
    
w = open(os.path.join(directory, "left.tsv"), 'w');
for nl in new_chiparts:
	larr = nl.read_id, nl.ref_id, nl.read_start, nl.read_end, nl.read_length, nl.ref_start, nl.ref_end, nl.ref_length, nl.cut_left, nl.cut_right, "|".join(nl.conversions), nl.weight, nl.score
	w.write("\t".join([str(x) for x in larr]) + "\n");	
w.close()	

w = open(os.path.join(directory, "fam_ids.tsv"), 'w');
for f, m in fams2mir.items():
	w.write("%s\t%s\n" % (f, ",".join(m)));
w.close()	












	