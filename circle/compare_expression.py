#! /usr/bin/python
'''Compares expression of circular RNAs between two conditions (replicates are allowed)'''
import argparse
import sys;
from collections import defaultdict
from math import log

import numpy
import pysam;
from pybedtools import BedTool
from nrlbio.pybedtools_extension import construct_gff_interval


parser = argparse.ArgumentParser(description='Compares expression of circular RNAs between two conditions (replicates are allowed)');
parser.add_argument('-c1', '--condition1', nargs = '+', required = True, type = str, help = "Path to the circular RNA expression files related to the replicates of the first experiment (condition), gff file");
parser.add_argument('-c2', '--condition2', nargs = '+', required = True, type = str, help = "Path to the circular RNA expression files related to the replicates of the second experiment (condition), gff file");
parser.add_argument('--names', nargs = 2, default = ('condition_1', 'condition_2'), type = str, help = "Descriptive names for the first and second conditions");
parser.add_argument('--minexpr', nargs = '?', default = 1, type = float, help = "Min number of reads supporting splice junctions threshold for circRNA to be analysed");
args = parser.parse_args();


size1 = len(args.condition1);
size2 = len(args.condition2);


circles = defaultdict(list)

cexpr1 = defaultdict(dict);
for r, path in enumerate(args.condition1):
	for interval in BedTool(path):
		cexpr1[(interval.chrom, interval.start, interval.end, interval.strand)][r] = float(interval.attrs['n_uniq']), float(interval.attrs['norm_expr'])
		circles[(interval.chrom, interval.start, interval.end, interval.strand)].append(float(interval.attrs['n_uniq']))
		

cexpr2 = defaultdict(dict);
for r, path in enumerate(args.condition2):
	for interval in BedTool(path):
		cexpr2[(interval.chrom, interval.start, interval.end, interval.strand)][r] = float(interval.attrs['n_uniq']), float(interval.attrs['norm_expr'])
		circles[(interval.chrom, interval.start, interval.end, interval.strand)].append(float(interval.attrs['n_uniq']))	
	

count=0
diff = []
for circle, expr in circles.iteritems():
	if(max(expr) >= args.minexpr):
		count+=1;
		support1 = []
		e1 = []
		d = cexpr1[circle]
		for r in range(size1):
			s, e = d.get(r, (0, 0.0))
			support1.append(s);
			e1.append(e);
			
		support2 = []
		e2 = []
		d = cexpr2[circle]
		for r in range(size2):
			s, e = d.get(r, (0, 0.0))
			support2.append(s);
			e2.append(e);
			
		mean1 = sum(e1)/size1 + 1
		mean2 = sum(e2)/size2 + 1
		lfc = log(mean1/mean2, 2)
		
		sys.stdout.write(str(construct_gff_interval(circle[0], circle[1], circle[2], 'circ', score='0', strand=circle[3], source='un', frame='.', attrs=[("ID", "c%d" % count), ('lfc2', "%1.4f" % lfc), ('csj_reads_%s_' % args.names[0], ",".join([str(x) for x in support1])), ('csj_reads_%s_' % args.names[1], ",".join([str(x) for x in support2]))] )))
		
		
		
		#diff.append((circle[0], circle[1], circle[2], count, circle[3], lfc, "\t".join([str(x) for x in support1]), "\t".join([str(x) for x in support2])))
  

#print "#chrom\tstart\tend\tname\tscore\tstrand\tlfc(%s/%s)2\t%s\t%s" % (args.names[0], args.names[1], "\t".join(['csj_reads_%s_%d' % (args.names[0], x+1) for x in range(size1)]), "\t".join(['csj_reads_%s_%d' % (args.names[1], x+1) for x in range(size2)]))

#diff.sort(key = lambda x: x[5], reverse=True) 
#for dcirc in diff:  
	#print "%s\t%d\t%d\tc%d\t0\t%s\t%1.2f\t%s\t%s" % dcirc
		
		
		
		
		
			
			
			
			
			
			
			
			
	