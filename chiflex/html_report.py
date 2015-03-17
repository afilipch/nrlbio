#! /usr/bin/python	
'''Script produces detailed html report for RNA-RNA interactions'''
import argparse;
import os;
import sys;
from collections import defaultdict;

import pysam
import jinja2

from nrlbio.generators import generator_doublebed
from nrlbio.samlib import ArWrapper
from nrlbio.interaction import Interaction


system_choices = ['hg19', 'hg38', 'mm9', 'ce6']


parser = argparse.ArgumentParser(description='Script produces detailed html report for RNA-RNA interactions');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to interactions");
parser.add_argument('--table', nargs = '+', type = str, help = "path to the table, which connects interaction ids and mapping hits(reads) ids. Normally can be found as \'interactions/interactions2readid.bed\'");
parser.add_argument('--sam', nargs = '+', type = str, help = "path to the sam/bam file with hits composing the interactions. Normally can be found as \'sam/[name].unique_chimera.bam\'");
parser.add_argument('--system', nargs = '?', required = True, choices = system_choices, type = str, help = "genome system. Can be set to %s" % "|".join(system_choices))
#parser.add_argument('-o', '--output', nargs = '?', default = None, type = str, help = "Path to the output report. If not set, report will be output to STDOUT");
parser.add_argument('-r', '--reference', nargs = '+', default = None, type = str, help = "path to the references used to extract sequences. If not set, '--extensions' will be set to None");
parser.add_argument('-e', '--extensions', nargs = '+', default = [], type = str, help = "Controls how far interacting intervals' sequences will be extended. For example extensions= 1,4 0,7 will extend first interval 1nt upstream and 4nts downstream, while the second interaction will be extended 0nt upstream and 7nts downstream. If '--reference' is not set, willl be set to None");
parser.add_argument('--template', nargs = '?', default = "interaction_detailed.html", type = str, help = "path to jinja2/django template");
parser.add_argument('--reassign', nargs = '?', default = False, const = True, type = bool, help = "Has to be set, if sam/bam hits have to be reassigned to genomic coordinates. That is, if nongenomic reference was used for mapping");
args = parser.parse_args();

#set jinja2 environment
environment = jinja2.Environment(loader=jinja2.PackageLoader('nrlbio', 'templates'))
template = environment.get_template(args.template)




#parse extensions
extensions = []
for e in args.extensions:
	extensions.append([int(x) for x in e.split(",")])
##################################################################################


#get reference dictionary, if needed
if(args.reference):
	from nrlbio.generators import generator_seqrecord
	from Bio import SeqIO
	reference = SeqIO.to_dict(generator_seqrecord(args.reference, 'fasta'));
else:
	reference = None
	extensions = None
##################################################################################


#connect sam/bam hits to the interactions
hid2iid = {};
for table in args.table:
	with open(table) as f:
		for l in f:
			iid, hids = l.strip().split("\t");
			for hid in hids.split(","):
				hid2iid[hid] = iid;
				
iid2hits = defaultdict(list);
for sampath in args.sam:
	samfile = pysam.Samfile(sampath);
	for aligned_read in samfile.fetch(until_eof=True):
		rname = samfile.getrname(aligned_read.tid)
		arw = ArWrapper(aligned_read, rname, add_nr_tag=False);
		arw.reassign_coordinates(reassign=args.reassign);
		iid = hid2iid.get(arw.qname);
		if(iid):
			iid2hits[iid].append(arw);
		#print "\t".join([str(x) for x in (arw.rname, arw.chrom, arw.start, arw.stop, arw.is_reverse)]);
		#print "\t".join([str(x) for x in (arw.aligned_read.pos, arw.aligned_read.aend)]);
		#print "_"*120
	samfile.close();	
##################################################################################


#compile and process interactions
for ipath in args.path:	
	for intervals in generator_doublebed(ipath):
		name = intervals[0].name.split("|")[0]
		inter = Interaction(name, intervals, aligned_reads = iid2hits[name])
		inter.set_extended_intervals(reference=reference, extensions=extensions)
		inter.set_html_attributes(args.system)
		print template.render({"interaction": inter})
		#sys.exit()
