#! /usr/lib/python
'''Assignes to each read, if it comes from mapping to real or control reference. Convolutes backward converted reads. Filters out non-unique and not the best hits'''
import argparse
import os;

import pysam;

from nrlbio.statistics import sam as sam_statistics
from nrlbio.samlib import ArWrapper, demultiplex_read_hits




parser = argparse.ArgumentParser(description='Assignes to each read if it comes from mapping to real or control reference. Convolutes backward converted reads. Filters out non-unique and not the best hits');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
parser.add_argument('-o', '--output', nargs = '?', default = "sam", type = str, help = "path to the output folder");
parser.add_argument('-r', '--report', nargs = '?', default = "reports", type = str, help = "path to the report folder");
parser.add_argument('-n', '--name', nargs = '?', required = True, type = str, help = "name for output files, should reflect nature of mapping reference");
parser.add_argument('-s', '--score', nargs = '?', default = 'as', choices = ['as', 'as_qstart', 'as_qstart_entropy', 'as_qstart_rstart', 'as_qstart_rstart_entropy'], type = str, help = "score function for hits");
parser.add_argument('--bestdistance', nargs = '?', default = 10, type = float, help = "minimal distance allowed between the best and the second best hit. If the actual distance is less, than hit will be assigned as nonunique");
args = parser.parse_args();



def _iteration(arwlist, bestdistance):
	if(arwlist):
		unique, nonunique, control = demultiplex_read_hits(arwlist, bestdistance)
		if(unique):
			sam_unique.write(unique.aligned_read);
			stat_unique.increment_basic(unique.aligned_read)
			stat_unique.increment_short(unique.aligned_read)
		elif(control):
			sam_control.write(control.aligned_read)
			stat_control.increment_basic(control.aligned_read)
			stat_control.increment_short(control.aligned_read)
		for nu in nonunique:
			sam_nonunique.write(nu.aligned_read);
			stat_nonunique.increment_basic(nu.aligned_read)
			stat_nonunique.increment_short(nu.aligned_read)



			
# parse input parameters
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
arwlist = [];
current_name = '';
for aligned_read in samfile.fetch(until_eof=True):
	if(not aligned_read.is_unmapped):
		
		rname = samfile.getrname(aligned_read.tid)
		arw = ArWrapper(aligned_read, rname, score_function=key_score, add_nr_tag=True)
		
		if(current_name != arw.qname):
			_iteration(arwlist, args.bestdistance)
			arwlist = [arw];
			current_name = arw.qname;
			
		else:
			arwlist.append(arw);
else:
	_iteration(arwlist, args.bestdistance);
			
			
# create html reports from statistics			
stat_unique.tohtml(os.path.join(args.report, "%s.unique.html" % args.name))
stat_control.tohtml(os.path.join(args.report, "%s.control.html" % args.name))
stat_nonunique.tohtml(os.path.join(args.report, "%s.nonunique.html" % args.name))
			
			