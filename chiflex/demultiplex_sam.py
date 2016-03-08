#! /usr/lib/python
'''Assignes to each read, if it comes from mapping to real or control reference. Convolutes backward converted reads. Filters out non-unique and not the best hits'''
import argparse
import os;
import sys

import pysam;

from nrlbio.statistics import sam as sam_statistics
from nrlbio.samlib import demultiplex_read_hits
from nrlbio.generators import generator_segments




parser = argparse.ArgumentParser(description='Assignes to each read if it comes from mapping to real or control reference. Convolutes backward converted reads. Filters out non-unique and not the best hits');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
parser.add_argument('-o', '--output', nargs = '?', default = "sam", type = str, help = "path to the output folder");
parser.add_argument('-r', '--report', nargs = '?', default = "reports", type = str, help = "path to the report folder");
parser.add_argument('-n', '--name', nargs = '?', required = True, type = str, help = "name for output files, should reflect nature of mapping reference");
parser.add_argument('--backward', nargs = '?', default = False, const=True, type = bool, help = "This flag has to be set, if reads were backward converted");
parser.add_argument('--collapsed', nargs = '?', default = False, const=True, type = bool, help = "This flag has to be set, if reads were collapsed");
parser.add_argument('-s', '--score', nargs = '?', default = 'as', choices = ['as', 'as_qstart', 'as_qstart_entropy', 'as_qstart_rstart', 'as_qstart_rstart_entropy'], type = str, help = "score function for hits");
parser.add_argument('--bestdistance', nargs = '?', default = 10, type = float, help = "minimal distance allowed between the best and the second best hit. If the actual distance is less, than hit will be assigned as nonunique");
args = parser.parse_args();

#print args.backward
#print args.collapsed

# parse [score] argument
score2function = {'as': 'as_score', 'as_qstart': 'as_qstart_score', 'as_qstart_entropy': 'as_qstart_entropy_score', 'as_qstart_rstart': 'as_qstart_rstart_score', 'as_qstart_rstart_entropy': 'as_qstart_rstart_entropy_score'}
exec("from nrlbio.samlib import %s as key_score" % score2function[args.score]);


# open output and input sam/bam files
samfile = pysam.Samfile(args.path)
sam_unique = pysam.Samfile(os.path.join(args.output, "%s.unique.bam" % args.name), "wb", template=samfile)
sam_control = pysam.Samfile(os.path.join(args.output, "%s.control.bam" % args.name), "wb", template=samfile)
sam_nonunique = pysam.Samfile(os.path.join(args.output, "%s.nonunique.bam" % args.name), "wb", template=samfile)


# create statistics for real, control and nonunique mappings
stat_unique = sam_statistics.Stat(name="Uniquely mapped")
stat_control = sam_statistics.Stat(name="Mapped to control reference")
stat_nonunique = sam_statistics.Stat(name="Nonuniquely mapped")


#mappings derived from the same read are pulled together. Collapsed into one (or more, in case of non-unique mappings) wrapping read object. 
#then they are written to a new destination, according to their source: real, or control
counts = [0]*4
for arwlist in generator_segments(args.path, key_score=key_score, add_nr_tag=args.collapsed, converted=args.backward):
	hits = demultiplex_read_hits(arwlist, args.bestdistance, backward=args.backward)
	if(len(hits)==1):
		hit = hits[0];
		if(hit.control):
			counts[1]+=1
			sam_control.write(hit.aligned_read)
			stat_control.increment_basic(hit.aligned_read)
			stat_control.increment_short(hit.aligned_read)
		else:
			counts[0]+=1
			sam_unique.write(hit.aligned_read);
			stat_unique.increment_basic(hit.aligned_read)
			stat_unique.increment_short(hit.aligned_read)
	else:
		snu = filter(lambda x: not x.control, hits)
		if(snu):
			counts[2]+=1
			for nu in snu:
				sam_nonunique.write(nu.aligned_read);
				stat_nonunique.increment_basic(nu.aligned_read)
				stat_nonunique.increment_short(nu.aligned_read)
		else:
			counts[3]+=1;


# create html reports from statistics
stat_unique.tohtml(os.path.join(args.report, "%s.unique.html" % args.name))
stat_control.tohtml(os.path.join(args.report, "%s.control.html" % args.name))
stat_nonunique.tohtml(os.path.join(args.report, "%s.nonunique.html" % args.name))
			
sys.stderr.write("number of unique hits: %d\nnumber of control hits: %d\nnumber of nonunique hits: %d\nnumber of control nonunique hits: %d\n\n" % tuple(counts));