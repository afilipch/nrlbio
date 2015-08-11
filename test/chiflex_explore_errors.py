#! /usr/lib/python
'''Provides detailed overview on reads which were mapped/demultiplexed incorrectly'''
import argparse
import os;
import sys;
from collections import defaultdict

from pysam import Samfile
from Bio.Seq import reverse_complement;

from nrlbio.samlib import ArWrapper

#from nrlbio.pyplot_extension import pie, histogram



parser = argparse.ArgumentParser(description='Provides detailed overview on reads which were mapped/demultiplexed incorrectly');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to unprocessed  hits(e.g. unprocessed Bowtie2 ouput). Sam/Bam format");
parser.add_argument('-e', '--errors', nargs = '?', required = True, type = str, help = "Errors plain file in tsv format")
parser.add_argument('-a', '--additional', nargs = '+', default = [], type = str, help = "Paths to processed hits of any kind(chimeras, single, nonunique and so on) required for testing. Sam/Bam format")
parser.add_argument('-t', '--etype', nargs = '+', default = [], type = str, help = "type of error to analyze. Has to be in a forma [read_type],[mapped_type], e.g. \'chimera,single\'")
parser.add_argument('-o', '--outdir', nargs = '?', default = 'evaluation_errors', type = str, help = "Name of the output directory. If directory exists, files will be (re)written there")
args = parser.parse_args();

if(not os.path.exists(args.outdir)):
    os.makedirs(args.outdir)
#________________________________________________________________________________________________________________
#get error dict;
#________________________________________________________________________________________________________________
errors = defaultdict(set)
with open(args.errors, 'r') as f:
	for l in f:
		a = l.strip().split("\t");
		if(not args.etype or ",".join(a[:2]) in args.etype):
			errors[tuple(a[:2])].add(a[2]);
			
errors2segments = defaultdict(lambda: defaultdict(list));
samfile = Samfile(args.path)
for segment in samfile.fetch(until_eof=True):
	num = segment.query_name.split("|")[0]
	for etype, eset in errors.iteritems():
		if(num in eset):
			errors2segments[etype][num].append(segment);
			break;
		
		
additional = defaultdict(list);
for fname in args.additional:
	tsamfile = Samfile(fname);
	for segment in tsamfile.fetch(until_eof=True):
		num = segment.query_name.split("|")[0]
		additional[num].append(ArWrapper(segment, tsamfile.getrname(segment.tid)))
	tsamfile.close();
		
		
		
		
for etype, d in errors2segments.iteritems():
	with open(os.path.join(args.outdir, "%s_%s_error.txt" % etype), 'w') as f:
		for num, segments in d.iteritems():
			if(segments[0].is_reverse):
				seq = reverse_complement(segments[0].seq);
			else:	
				seq = segments[0].seq
			
			f.write("%s\nnumber of read:\t%s\n\nSequence:\t%s\n\nSegments:\n\n" % ("_"*140, num, seq))
			for segment in segments:
				f.write("%s\t%s\t%d\t%s\t%d\t%s\n\n" % (segment.query_name.split("|")[2], samfile.getrname(segment.tid), segment.reference_start, segment.cigarstring, segment.get_tag("AS"), segment.query_name));
				
			annotated = additional[num];
			if(annotated):
				f.write("Annotated as:\n\n")
				for arw in annotated:
					segment = arw.aligned_read;
					f.write("%s\t%s\t%d\t%s\t%d\t%s\n\n" % (segment.query_name.split("|")[2], arw.rname, segment.reference_start, segment.cigarstring, segment.get_tag("AS"), segment.query_name));
			
			
#for k, v in errors.iteritems():
	#print k, len(v)