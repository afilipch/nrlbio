#! /usr/lib/python
'''Assignes to each read, if it comes from mapping to real or control reference. Convolutes backward converted reads. Filters out non-unique and not the best hits'''
import argparse
import os;
import sys;

import pysam;

import nrlbio.statistics.sam as sam_statistics
from nrlbio.samlib import ArWrapper, demultiplex_read_hits
from nrlbio.chimera import as_gap_score, arwlist2chimera
from nrlbio.chimera import demultiplex as demultiplex_ch



parser = argparse.ArgumentParser(description='Assignes to each read if it comes from mapping to real or control reference. Convolutes backward converted reads. Filters out non-unique and not the best hits');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
parser.add_argument('-o', '--output', nargs = '?', default = "sam", type = str, help = "path to the output folder");
#parser.add_argument('-r', '--report', nargs = '?', default = "reports", type = str, help = "path to the report folder");
parser.add_argument('-n', '--name', nargs = '?', required = True, type = str, help = "name for output files, should reflect nature of mapping reference");
parser.add_argument('-s', '--score', nargs = '?', default = 'as', choices = ['as', 'as_qstart', 'as_qstart_entropy', 'as_qstart_pos', 'as_qstart_pos_entropy'], type = str, help = "score function for hits");
parser.add_argument('-sh', '--score_chimera', nargs = '?', default = 'as', choices = ['as', 'as_gap', 'as_gap_entropy'], type = str, help = "score function for chimeras");
parser.add_argument('-mg', '--maxgap', nargs = '?', default = 8, type = int, help = "maxgap is used to calculate chimera score, the more the maxgap, the less important is gap between hits for chimera evaluation");
args = parser.parse_args();

testid = '22462|chimera|chr1:+:108451004:108451037&chr1:-:227677228:227677295|intron:intergenic|0:0:0|0:0'


def compare_single_chimera(arw, chimera, maxgap):
	'''Decides if the source of the read is continuous or not(chimeric)
	
		arw samlib.ArWrapper: wrapper for continuous hit
		chimera samlib.Chimera: wrapper for chimeric hits
		maxgap int: controls the stringency of gap/overlap requirement for chimeric reads. The bigger maxgap, the weaker requirement.
		
	Returns tuple (None, chimera)|(arw|None): returns input arw or chimera depending on their comparison
	'''
	chs = as_gap_score(chimera, maxgap=maxgap)
	sis = arw.AS*(0.75 + 0.5**min(5, (arw.aligned_read.qstart+1)))
	if(sis>=chs):
		return arw, None
	else:	
		return None, chimera
		
		

counts = [0]*6

def _iteration(arwlist):
	#if(arwlist):
		#if(arwlist[0].qname == testid):
			#for arw in arwlist:
				#segment = arw.aligned_read;
				#print "%s\t%s\t%d\t%s\t%d\t%s\n" % (segment.query_name.split("|")[2], arw.rname, segment.reference_start, segment.cigarstring, segment.get_tag("AS"), segment.query_name)
			#print 
			
			#unique, nonunique, control = demultiplex_read_hits(arwlist, key_score);
			#segment = unique.aligned_read;
			#print "BEST UNIQUE", "%s\t%s\t%d\t%s\t%d\t%s\n" % (segment.query_name.split("|")[2], unique.rname, segment.reference_start, segment.cigarstring, segment.get_tag("AS"), segment.query_name)
			
			#chimeras = arwlist2chimera(arwlist, gap = 0, overlap = 6, score_function = key_score_chimera)
			#print "\npossible chimeras: %d\n" % len(chimeras);
			
			#unique_chimera, nonunique_chimera, control_chimera = demultiplex_ch(chimeras);
			#for arw in unique_chimera.ar_wrappers:
				#segment = arw.aligned_read;
				#print "BEST CHIMERA", "%s\t%s\t%d\t%s\t%d\t%s\n" % (segment.query_name.split("|")[2], arw.rname, segment.reference_start, segment.cigarstring, segment.get_tag("AS"), segment.query_name)
			#print	
				
				
				
				
		#demultiplex single hits
		unique, nonunique, control = demultiplex_read_hits(arwlist, key_score);
		
		#gather hits into chimeras
		chimeras = arwlist2chimera(arwlist, gap = 0, overlap = 6, score_function = key_score_chimera)
		unique_chimera, nonunique_chimera, control_chimera = demultiplex_ch(chimeras);
		
		#assign read to be chimeric or single separately for control and signal
		if(unique and unique_chimera):
			unique, unique_chimera = compare_single_chimera(unique, unique_chimera, args.maxgap)
		elif(unique and control_chimera):
			unique, control_chimera = compare_single_chimera(unique, control_chimera, args.maxgap)
		elif(control and unique_chimera):
			control, unique_chimera = compare_single_chimera(control, unique_chimera, args.maxgap)
		elif(control and control_chimera):
			control, control_chimera = compare_single_chimera(control, control_chimera, args.maxgap)
			
		#output single reads	
		if(unique):
			counts[0] += 1
			sam_unique.write(unique.aligned_read);
		elif(control):
			counts[1] += 1
			sam_control.write(control.aligned_read)
		if(nonunique):
			counts[2] += 1
			for nu in nonunique:
				sam_nonunique.write(nu.aligned_read);
			
		#output chimeras
		if(unique_chimera):
			counts[3] += 1
			for arw in unique_chimera.ar_wrappers:
				sam_unique_chimera.write(arw.aligned_read);
		elif(control_chimera):
			counts[4] += 1
			for arw in control_chimera.ar_wrappers:
				sam_control_chimera.write(arw.aligned_read);
		if(nonunique):
			counts[5] += 1
			for nu in nonunique_chimera:
				for arw in nu.ar_wrappers:
					sam_nonunique_chimera.write(arw.aligned_read);



# parse input parameters
score2function = {'as': 'as_score', 'as_qstart': 'as_qstart_score', 'as_qstart_entropy': 'as_qstart_entropy_score', 'as_qstart_pos': 'as_qstart_pos_score', 'as_qstart_pos_entropy': 'as_qstart_pos_entropy_score'}
exec("from nrlbio.samlib import %s as key_score" % score2function[args.score]);

score_2function_chimera = {'as': 'as_score', 'as_gap': 'as_gap_score', 'as_gap_entropy': 'as_gap_entropy_score'}
exec("from nrlbio.chimera import %s as key_score_chimera" % score_2function_chimera[args.score_chimera]);



# open input sam/bam file 
samfile = pysam.Samfile(args.path)

#open output bam files for single reads
sam_unique = pysam.Samfile(os.path.join(args.output, "%s.unique.bam" % args.name), "wb", template=samfile)
sam_control = pysam.Samfile(os.path.join(args.output, "%s.control.bam" % args.name), "wb", template=samfile)
sam_nonunique = pysam.Samfile(os.path.join(args.output, "%s.nonunique.bam" % args.name), "wb", template=samfile)

#open output bam files for chimeric reads
sam_unique_chimera = pysam.Samfile(os.path.join(args.output, "%s.unique_chimera.bam" % args.name), "wb", template=samfile)
sam_control_chimera = pysam.Samfile(os.path.join(args.output, "%s.control_chimera.bam" % args.name), "wb", template=samfile)
sam_nonunique_chimera = pysam.Samfile(os.path.join(args.output, "%s.nonunique_chimera.bam" % args.name), "wb", template=samfile)



#mappings derived from the same read are pulled together. Collapsed into one (or more, in case of non-unique mappings) wrapping read object. 
#then they are written to a new destination, according to their source: real, or control
arwlist = [];
current_name = '';
for aligned_read in samfile.fetch(until_eof=True):
	if(not aligned_read.is_unmapped):
		
		rname = samfile.getrname(aligned_read.tid)
		arw = ArWrapper(aligned_read, rname, add_nr_tag=False)
		
		if(current_name != arw.qname):
			_iteration(arwlist)
			arwlist = [arw];
			current_name = arw.qname;
			
		else:
			arwlist.append(arw);
else:
	_iteration(arwlist);
	
	
sys.stderr.write("number of unique hits: %d\nnumber of control hits: %d\nnumber of nonunique hits: %d\n\nnumber of unique chimeras: %d\nnumber of control chimeras: %d\nnumber of nonunique chimeras: %d\n\n" % tuple(counts));
			

			
			