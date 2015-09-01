#! /usr/bin/python	
'''Script produces detailed html report for interactions recovered'''
from collections import *;
import argparse;
import os;
import sys;
import copy;
import chipart_lib;
import html_lib;
import interaction_lib;
import gconst;
import jinja2





parser = argparse.ArgumentParser(description='Script produce detailed html overview of interactions');
# input files
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to read2int table(output/read2int.tsv)");
parser.add_argument('-i', '--interactions', nargs = '?', type = str, help = "path to interactions (output/interactions.bed)");
parser.add_argument('-l', '--left', nargs = '+', type = str, help = "path to filtered left parts(left/filtered_chiparts.tsv or gmode/chiparts.tsv)");
parser.add_argument('-r', '--right', nargs = '+', type = str, help = "path to filtered right parts(sam/right.tsv)");
parser.add_argument('-f', '--fastq', nargs = '+', type = str, help = "path to potential chimeras in fastq format(left/filtered_candidates.fastq)");
parser.add_argument('-t', '--template', nargs = '?', default = "detailed.html", type = str, help = "path to jinja2 template in html format");
parser.add_argument('--title', nargs = '?', default = "interactions", type = str, help = "html title");
#parser.add_argument('-t', '--template', nargs = '?', default = os.path.join(gconst.troot, "detailed.html"), type = str, help = "path to jinja2 template in html format");
parser.add_argument('-s', '--system', nargs = '?', required = True, type = str, help = "genome system")
args = parser.parse_args();

Right = namedtuple("Right", 'chromosome, genome_start, genome_end, read_id, score, strand, read_start, read_end, collapse, conversion, aligned')

system = args.system
exec("from sequence_data.systems import %s as gsys" % system);

rid2rparts = {};
for p in args.right:
	f = open(p);
	for l in f:
		a = l.strip().split("\t");
		rid2rparts[a[3]] = Right(a[0],int(a[1]), int(a[2]), a[3], int(a[4]), a[5], int(a[6]), int(a[7]), int(a[8]), int(a[9]), a[10])
	f.close();


chiparts = chipart_lib.read(args.left); 
rid2chipart = {};
for ch in chiparts:
	rid2chipart[ch.read_id] = ch;

	

int2reads = defaultdict(set);
f = open(args.path[0]);
for l in f:
	a = l.strip().split("\t");
	int2reads[a[1]].add(a[0]);
f.close();

rid2seq = {};
for p in args.fastq:
	f = open(p);
	lset = [];
	for l in f:
		lset.append(l.strip());
		if(len(lset) == 4):
			rid2seq[lset[0]] = lset[1];
			lset = [];
	f.close();


interactions = interaction_lib.get(args.interactions, undef = True);
interactions.sort(key = lambda x: x.mirid)
html_interactions = [];
depth = 11 


for interaction in interactions:
	#sys.stderr.write(str(interaction) + "\n")
	lparts = [];
	rparts = [];
	seqs = []; 
	rids = [];
	for rid in int2reads[interaction.iid]:
		lparts.append(rid2chipart[rid]);
		rparts.append(rid2rparts[rid]);
		seqs.append(rid2seq[rid]);
		rids.append(rid);
	refseq = gsys.genome.get_oriented(interaction.chromosome, interaction.start - depth, interaction.end + depth, interaction.strand).upper()	
	
	#if( interaction.iid == 'ago2_kishore-pc01495-2'):
		#for lp in lparts:
			#sys.stderr.write(str(lp) + "\n");		
		#t = html_lib.HTML_interaction(interaction, lparts, rparts, seqs, rids, args.system, refseq, depth)
		#for read in t.reads:
			#sys.stderr.write(read.left_mapped + "\t"+ str(read.unmapped) + "\t"+ str(read.right_mapped) + "\n");
			
	html_interactions.append(html_lib.HTML_interaction(interaction, lparts, rparts, seqs, rids, args.system, refseq, depth));
#sys.stderr.write("%d\t%d\n" % (len(html_interactions),len(interactions)));
	
env = jinja2.Environment(loader=jinja2.FileSystemLoader(gconst.troot))
t = env.get_template(args.template)		
print t.render({"title": args.title, "interactions": html_interactions})	





