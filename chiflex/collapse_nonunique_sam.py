#! /usr/lib/python
'''collapse nonunique sam records. The group of hits with identical alignment will be collapsed into one of them. All locations on reference corresponding to the identical alignments will be stored in separate bed file''' 
import argparse
import sys;

import pysam;

from nrlbio.samlib import remove_duplicates, as_score
from nrlbio.generators import generator_segments


parser = argparse.ArgumentParser(description='collapse nonunique sam records. The group of hits with identical alignment will be collapsed into one of them. All locations on reference corresponding to the identical alignments will be stored in separate bed file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to sam/bam file");
parser.add_argument('-s', '--output_sam', nargs = '?', required = True, type = str, help = "path to output collapsed sam file");
parser.add_argument('-b', '--output_bed', nargs = '?', required = True, type = str, help = "path to output bed file");
parser.add_argument('-m', '--minscore', nargs = '?', default = 0, type = float, help = "min alignment score to collapse hits");
args = parser.parse_args();





def _iteration(arwlist, cid):
	result = remove_duplicates(arwlist, as_score, minscore=args.minscore)
	if(result):
		arw, bed = result
		arw.aligned_read.tags = arw.aligned_read.tags + [("CI", collapsed)];	
		output_sam.write(arw.aligned_read);
		for i in bed:
			output_bed.write("%s\t%d\t%d\t%d\t0\t%s\n" % (i[0], i[1], i[2], collapsed, i[3]));
		return 1
	else:
		for arw in arwlist:
			arw.aligned_read.tags = arw.aligned_read.tags + [("CI", 0)];
			output_sam.write(arw.aligned_read);
		return 0
		

samfile = pysam.Samfile(args.path)
output_sam = pysam.Samfile(args.output_sam, "wb", template=samfile)


total, collapsed = 0,0
		
	

	

with open(args.output_bed, 'w') as output_bed:
	for arwlist in generator_segments(args.path):
		collapsed+= _iteration(arwlist, collapsed);
		total +=1 
	
	
sys.stderr.write("total hits: %d\ncollapsed nonunique hits: %d\n" % (total, collapsed))


