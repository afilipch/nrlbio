#! /usr/lib/python
'''collapse nonunique sam records. The group of hits with identical alignment will be collapsed into one of them. All locations on reference corresponding to the identical alignments will be stored in separate bed file''' 
import argparse
import sys;

import pysam;

from nrlbio.samlib import ArWrapper, remove_duplicates, as_score



parser = argparse.ArgumentParser(description='collapse nonunique sam records. The group of hits with identical alignment will be collapsed into one of them. All locations on reference corresponding to the identical alignments will be stored in separate bed file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
parser.add_argument('-s', '--output_sam', nargs = '?', required = True, type = str, help = "path to output collapsed sam file");
parser.add_argument('-b', '--output_bed', nargs = '?', required = True, type = str, help = "path to output bed file");
parser.add_argument('-m', '--minscore', nargs = '?', default = 0, type = float, help = "min alignment score to collapse hits");
args = parser.parse_args();





def _iteration(arwlist):
	result = remove_duplicates(arwlist, as_score, minscore=args.minscore)
	if(result):
		arw, bed = result
		arw.aligned_read.tags = arw.aligned_read.tags + [("CI", collapsed+1)];	
		output_sam.write(arw.aligned_read);
		for i in bed:
			if(i[3]):
				strand = '-'
			else:
				strand = '+'
			output_bed.write("%s\t%d\t%d\t%d\t0\t%s\n" % (i[0], i[1], i[2], collapsed+1, strand));
		return 0, 1
	else:
		for arw in arwlist:
			arw.aligned_read.tags = arw.aligned_read.tags + [("CI", 0)];	
			output_sam.write(arw.aligned_read);
		return 1, 0	
		

samfile = pysam.Samfile(args.path)
output_sam = pysam.Samfile(args.output_sam, "wb", template=samfile)
output_bed = open(args.output_bed, 'w')

noncollapsed, collapsed = 0,0		
		
		
arwlist = [];
current_name = '';
for aligned_read in samfile.fetch(until_eof=True):
	if(not aligned_read.is_unmapped):
		rname = samfile.getrname(aligned_read.tid)
		arw = ArWrapper(aligned_read, rname, add_nr_tag=True)
		
		if(current_name != arw.qname):
			r = _iteration(arwlist);
			noncollapsed+=r[0] 
			collapsed+=r[1] 
			arwlist = [arw];
			current_name = arw.qname;
			
		else:
			arwlist.append(arw);
else:
	r = _iteration(arwlist);
	noncollapsed+=r[0] 
	collapsed+=r[1] 
	
	
output_bed.close();	
sys.stderr.write("total hits: %d\ncollapsed nonunique hits: %d\n" % (noncollapsed + collapsed, collapsed))


