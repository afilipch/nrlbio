#! /usr/bin/python
'''Exctracts an expression of small RNAs from sam/bam file'''
import argparse
import sys;
from collections import Counter, defaultdict

import pysam;
from pybedtools import BedTool

parser = argparse.ArgumentParser(description='Exctracts an expression of small RNAs from sam/bam file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to sam/bam file");
parser.add_argument('-mp', '--mirbase_precursors', nargs = '?', required = True, type = str, help = "Path to the mirbase mirna precursors gff file");
parser.add_argument('--collapsed', nargs = '?', default = False, const = True, type = str, help = "If set, reads are considered to be sequence-collapsed");
args = parser.parse_args();

def collapsed_number(segment):
	return int(segment.qname.split("_")[-1][1:])

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
		
		#righttail = Counter([x.query_sequence[x.qend:] for x in reflist if x.reference_end == len(interval)]).most_common(5)
		#righttail = [list(x) for x in righttail]
		#for el in righttail:
			#if(el[0] == ''):
				#el[0] = 'noadd'
		#interval.attrs['righttail'] = ",".join(["%s:%s" % tuple(x) for x in righttail])
	
	else:
		sys.stderr.write("WARNING: miRNA %s was not found in mirbase precursors provided\n" % rname);
	
	return interval
	
	
def get_stat_collapsed(reflist, rname):
	interval = mirnas.get(rname)
	
	if(interval):
		interval.attrs['raw_expr'] = str(sum([collapsed_number(x) for x in reflist]))
		
		cutleft = defaultdict(int)
		for sg in reflist:
			cutleft[sg.reference_start] += collapsed_number(sg)
			
		cutleft = list(sorted(cutleft.items(), key = lambda x: x[1], reverse=True))[:5]
		interval.attrs['cutleft'] = ",".join(["%d:%d" % x for x in cutleft])
		
		#righttail = defaultdict(int)
		#for sg in reflist:
			#righttail[sg.reference_start] = collapsed_number(sg)
		#righttail = list(sorted(righttail.items(), key = x[1], reverse=True))[:5]
		#interval.attrs['righttail'] = ",".join(["%s:%s" % tuple(x) for x in righttail])
	
	else:
		sys.stderr.write("WARNING: miRNA %s was not found in mirbase precursors provided\n" % rname);
	
	return interval	
	
	

	
samfile = pysam.Samfile(args.path)
reflist = [];
current_name = '';
totalreads = 0


expression = [];

if(args.collapsed):
	for segment in samfile.fetch(until_eof=True):
		if(not segment.is_unmapped):
			totalreads += collapsed_number(segment);
			rname = samfile.getrname(segment.tid)
			if(current_name and current_name != rname):
				expression.append(get_stat_collapsed(reflist, current_name))
				reflist = [segment];
			else:
				reflist.append(segment);
			current_name = rname
		else:
			sys.stderr.write('NNNN\n')
			if(current_name):
				expression.append(get_stat_collapsed(reflist, current_name));
				reflist = [];
				
			current_name = '';
				
	else:
		expression.append(get_stat_collapsed(reflist, current_name))
else:				   
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


	
	