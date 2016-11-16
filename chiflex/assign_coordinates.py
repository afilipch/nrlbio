#! /usr/lib/python
'''assigns genomic to the chimeric reads which were mapped to non-genomic reference(transcriptome, 3'utr, etc.)'''
import argparse;
import sys;

from pybedtools import BedTool

parser = argparse.ArgumentParser(description='assigns genomic coordinates to the regions which are on non-genomic reference(transcriptome, 3\'utr, etc.)');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to bed/gff file");
parser.add_argument('--quite', nargs = '?', const=True, default = False, type = bool, help = "If set, warning messages will be supressed. Useful for mirna:target identification");
args = parser.parse_args();

def reassign(interval):
	try:
		chrom, strand, start, stop = interval.chrom.split("|")[:4]
		#print "_"*120
		#sys.stdout.write(str(interval));
	except:
		#print "_"*120
		if(not args.quite):
			sys.stderr.write('Interval chrom has to be in [chrom]|[strand]|[start]|[stop] format. The actual name is %s\n' % interval.chrom)
		return interval
		
	interval.chrom 	= chrom
	interval.strand = strand
	
	if(strand == '+'):
		interval.stop = int(start) + interval.stop
		interval.start = int(start) + interval.start
	else:
		tstop = interval.stop
		interval.stop = int(stop) - interval.start
		interval.start = int(stop) - tstop

	return interval


for interval in BedTool(args.path):
	sys.stdout.write(str(reassign(interval)));

