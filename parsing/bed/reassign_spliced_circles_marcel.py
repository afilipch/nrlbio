#! /usr/lib/python
'''assigns genomic coordinates to the regions which are on non-genomic reference(transcriptome, 3\'utr, etc.) and splice(frst is was used for Marcel\'s data)'''
import argparse;
import sys;

from pybedtools import BedTool

parser = argparse.ArgumentParser(description='assigns genomic coordinates to the regions which are on non-genomic reference(transcriptome, 3\'utr, etc.) and splice(frst is was used for Marcel\'s data)');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to bed/gff file");
parser.add_argument('--quite', nargs = '?', const=True, default = False, type = bool, help = "If set, warning messges will be supressed. Useful for mirna:target identification");
args = parser.parse_args();

def reassign(interval):
	try:
		name, chrom, starts, stops, strand, score. length = interval.chrom.split("|")[:7]
	except:
		#print "_"*120
		if(not args.quite):
			sys.stderr.write('Interval chrom has to be in [name]|[chrom]|[starta]|[stops]|[strand]|[score]|[length] format. The actual name is %s\n' % interval.chrom)
		return interval
	
	srev = bool(strand=='-')
	starts = list(sorted([int(x) for x in starts], reverse = srev))
	stops = list(sorted([int(x) for x in stops], reverse = srev))
	
	curlen = 0;
	for start, end in zip(starts, stops):
		curlen += end - start;
		if(interval.start < curlen):
			return reassign_to_exon(interval, chrom, strand, start, stop)
	
	
	
def reassign_to_exon(interval, chrom, strand, start, stop)
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