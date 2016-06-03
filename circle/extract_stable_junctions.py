#! /usr/bin/python
'''Exctracts those exon-exon linear junctions, whose parental transcripts are not changed in expression between different conditions(replicates)'''


import sys;
import argparse
from collections import defaultdict, namedtuple
from itertools import combinations

from pybedtools import BedTool
import numpy



parser = argparse.ArgumentParser(description='Exctracts those exon-exon linear junctions, whose parental transcripts are not changed in expression between different conditions(replicates)');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the ensembl annotation, gff3 format");
parser.add_argument('--replicates', nargs = '+', required = True, type = str, help = "Path to the expression of transcripts, kallisto formated tsv files of replicates");
parser.add_argument('--maxvarcoeff', nargs = '?', default = 0.1, type = float, help = "Maximum variation coefficient allowed");
parser.add_argument('--maxdev', nargs = '?', default = 0.2, type = float, help = "Maximum deviation from a mean allowed");
parser.add_argument('--minexpr', nargs = '?', default = 50, type = float, help = "Min expression allowed");
args = parser.parse_args();

rsize = len(args.replicates);

#Transcript = namedtuple('Transcript', ['name', 'meanexpr', 'varcoeff', 'maxdev'])

texpr = defaultdict(list);
for rpath in args.replicates:
	with open(rpath) as f:
		f.next()
		for l in f:
			a = l.strip().split("\t")
			texpr[a[0]].append(float(a[4]))
			
			
			
texons = defaultdict(list)
for interval in BedTool(args.path):
	if(interval[2] == 'exon'):
		tname = interval.attrs['Parent'].split(':')[1]
		texons[tname].append(interval);
		
		
junctions_expr = defaultdict(lambda: numpy.zeros(rsize));
for tname, exons in texons.iteritems():
	expr = texpr.get(tname)
	if(expr):
		expr = numpy.array(expr);
		chrom, strand = exons[0].chrom, exons[0].strand
		for e1, e2 in zip(exons, exons[1:]):
			junctions_expr[(chrom, e1.end, e2.start, strand)] += expr
				

c = 0;
for junction, expr in junctions_expr.iteritems():
	mean = numpy.mean(expr)
	varcoeff = numpy.std(expr)/mean
	maxdev = max(mean-min(expr), max(expr)-mean)/mean
	if(varcoeff<args.maxvarcoeff and maxdev<args.maxdev and mean>args.minexpr):
		c+=1;
		print "chr%s\t%d\t%d\tj%d\t%1.2f\t%s" % (junction[0], junction[1], junction[2], c, mean, junction[3])
			
			
			
			
#stable_transcripts = {}; 

#for name, expr in texpr.iteritems():
	#if(all(expr)):
		#vals = numpy.array(expr)
		#mean = numpy.mean(vals)
		#varcoeff = numpy.std(vals)/mean
		#maxdev = max(mean-min(vals), max(vals)-mean)/mean
		
		#if(varcoeff<args.maxvarcoeff and maxdev<args.maxdev and mean>args.minexpr):
			#t = Transcript(name, mean, varcoeff, maxdev)
			#stable_transcripts[name] = t;
			
			
#texons = defaultdict(list);
#for interval in BedTool(args.path):
	#if(interval[2] == 'exon'):
		#tname = interval.attrs['Parent'].split(':')[1]
		#if(tname in stable_transcripts):
			#texons[tname].append(interval);
			
			
#for tname, exons in texons.iteritems():
	#exons.sort(key = lambda x: x.start)
	#chrom, strand = exons[0].chrom, exons[0].strand
	#score = stable_transcripts[tname].meanexpr
	##for exon in exons:
		##print exon;
	#for p, (e1, e2) in enumerate(zip(exons, exons[1:])):
		#print "chr%s\t%d\t%d\t%s|%d\t%1.2f\t%s" % (chrom, e1.end, e2.start, tname, p+1, score, strand)
		

			
			
			
			
			
			
			
			
			