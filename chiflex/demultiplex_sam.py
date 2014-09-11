#! /usr/lib/python
'''Assignes to each read, if it comes from mapping to real or control reference. Convolutes backward converted reads. Filters out non-unique and not the best hits'''
import argparse
import os;

import pysam;

from nrlbio import sam_statistics
from nrlbio import samlib
#import sys;
#import os;
#import argparse;
#import copy;
#from Bio import SeqIO
#from collections import *;
#from Bio.Seq import reverse_complement;
 
#import chimera_lib;
#import functions;



parser = argparse.ArgumentParser(description='Assignes to each read if it comes from mapping to real or control reference. Convolutes backward converted reads. Filters out non-unique and not the best hits');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
parser.add_argument('-o', '--output', nargs = '?', default = "sam", type = str, help = "path to the output folder");
parser.add_argument('-r', '--report', nargs = '?', default = "reports", type = str, help = "path to the report folder");
parser.add_argument('-n', '--name', nargs = '?', required = True, type = str, help = "name for output files, should reflect nature of mapping reference");
args = parser.parse_args();

# open output and input sam/bam files
samfile = pysam.Samfile(args.path)
sam_unique = pysam.Samfile(os.path.join(args.output, "%s.unique.bam" % args.name), "wb", template=samfile)
sam_control = pysam.Samfile(os.path.join(args.output, "%s.control.bam" % args.name), "wb", template=samfile)
sam_nonunique = pysam.Samfile(os.path.join(args.output, "%s.nonunique.bam" % args.name), "wb", template=samfile)

# create statistics for real, control and nonunique mappings
stat_unique = sam_statistics.Stat(name="Uniquely mapped")
stat_control = sam_statistics.Stat(name="Mapped to contol reference")
stat_nonunique = sam_statistics.Stat(name="Nonuniquely mapped")


#mappings derived from the same read are pulled together. Collapsed into one (or more, in case of non-unique mappings) wrapping read object. 
#then they are written to a new destination, according to their source: real, or control
arwlist = [];
current_name = '';
for aligned_read in samfile.fetch(until_eof=True):
	if(not aligned_read.is_unmapped):		
		rname = samfile.getrname(aligned_read.tid)
		arw = samlib.ArWrapper(aligned_read, rname) 
		
		if(current_name != arw.qname):			
			if(arwlist):
				unique, nonunique, control = samlib.demultiplex_read_hits(arwlist, samlib.key_alignment_score)
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
			arwlist = [arw];
			current_name = arw.qname;
		else:
			arwlist.append(arw);
			
# create html reports from statistics			
stat_unique.tohtml(os.path.join(args.report, "%s.unique.html" % args.name))
stat_control.tohtml(os.path.join(args.report, "%s.control.html" % args.name))
stat_nonunique.tohtml(os.path.join(args.report, "%s.nonunique.html" % args.name))
			
			