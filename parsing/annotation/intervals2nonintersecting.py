#! /usr/lib/python
'''script extracts certain features from annotation file and creates non-intersecting bed file with meaningful ids''' 
import argparse
import sys;
from collections import defaultdict

from pybedtools import BedTool, Interval;






parser = argparse.ArgumentParser(description='script extracts certain features from annotation file and creates non-intersectimg bed file with meaningful ids');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to gff annotation");
parser.add_argument('-t', '--transcription', nargs = '+', default = None, choices = ['exon', 'intron', 'None'], type = str, help = "type of transcription of intervals to choose");
parser.add_argument('-mt', '--regulation', nargs = '+', default = None, type = str, help = "type of regulation of intervals to choose");
args = parser.parse_args();

def type_filter(feature, transcription, regulation):
	if(transcription and feature.attrs['transcription'] not in transcription):
		return False;
	elif(regulation and feature.attrs['regulation'] not in regulation):
		return False;
	else:
		return True;


genes_features = defaultdict(list)

for i in BedTool(args.path).filter(type_filter, transcription=args.transcription, regulation=args.regulation):
	genes_features[i.attrs['gene_id']].append(i);
	
	


for c, (k, v) in enumerate(genes_features.iteritems()):
	b = BedTool((Interval(i.chrom, i.start, i.stop, strand = i.strand, name = k, score = "0")  for i in v))
	a = b.merge(s=True, d=0)
	sys.stderr.write("%d\n%d\n\n" % (c, len(a)))
	for i in a:
		sys.stdout.write(str(i))
	if(c>10000):
		sys.exit()	

