#! /usr/bin/python
'''Exctracts an expression of small RNAs from sam/bam file'''
import argparse
import sys;
from collections import Counter

import pysam;
from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Exctracts an expression of small RNAs from sam/bam file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to sam/bam file");
parser.add_argument('-mp', '--mirbase_precursors', nargs = '?', required = True, type = str, help = "Path to the mirbase mirna precursors gff file");
args = parser.parse_args();

#Read mirbase mirna precursors file
mirnas = {};
for interval in BedTool(args.mirbase_precursors):
	if(interval[2] == 'miRNA'):
		mirnas[interval.attrs['Name']] = interval


def get_stat(reflist, rname):
	interval = mirnas.get(rname)
	
	if(interval):
		interval.attrs['raw_expr'] = str(len(reflist))
		
		cutleft = Counter([x.reference_start for x in reflist]).most_common(5)
		interval.attrs['cutleft'] = ",".join(["%d:%d" % x for x in cutleft])
		
		righttail = Counter([x.query_sequence[x.qend:] for x in reflist if x.reference_end == len(interval)]).most_common(5)
		righttail = [list(x) for x in righttail]
		for el in righttail:
			if(el[0] == ''):
				el[0] = 'noadd'
		interval.attrs['righttail'] = ",".join(["%s:%s" % tuple(x) for x in righttail])
	
	else:
		sys.stderr.write("WARNING: miRNA %s was not found in mirbase precursors provided\n" % rname);
	
	return interval
	

	
samfile = pysam.Samfile(args.path)
reflist = [];
current_name = '';
totalreads = 0


expression = [];

for segment in samfile.fetch(until_eof=True):
	if(not segment.is_unmapped):
		totalreads += 1;
		rname = samfile.getrname(segment.tid)
		if(current_name and current_name != rname):
			expression.append(get_stat(reflist, current_name))
			reflist = [segment];
		else:
			reflist.append(segment);
		current_name = rname
	else:
		if(current_name):
			expression.append(get_stat(reflist, current_name));
			reflist = [];
			
		current_name = '';
			
else:
	expression.append(get_stat(reflist, current_name))
	
	
samfile.close();

expression = filter(bool, expression);
expression.sort(key = lambda x: float(x.attrs['raw_expr']), reverse = True)
totalreads = totalreads/1000000.0

for interval in expression:
	interval.attrs['norm_expr'] = '%1.2f' % (float(interval.attrs['raw_expr'])/totalreads)
	sys.stdout.write(str(interval))


	
	