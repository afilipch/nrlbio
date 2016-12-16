#! /usr/bin/python
'''Provides the longest exonic structure of given circles on basis of transcript annotation(gff3 file)'''


import sys;
import argparse
from collections import defaultdict

from pybedtools import BedTool, Interval

from nrlbio.pybedtools_extension import gff2bed


parser = argparse.ArgumentParser(description='Provides the longest exonic structure of given circles on basis of transcript annotation(gff3 file)');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the circles, bed/gff format");
parser.add_argument('--gff3', nargs = '?', required = True, type = str, help = "Path to the genome system annotation, gff3 format");
args = parser.parse_args();


def check_intersection(e1, e2):
	return (e1.end < e2.start) or (e2.end < e1.start)



def get_longest(exons):
	maxes = defaultdict(float);
	pathes = defaultdict(list);
	
	for i in range(len(exons)):
		localscores = [];
		for j in range(i):
			if(check_intersection(exons[i], exons[j])):
				localscores.append(len(exons[i]) + maxes[j])
			else:	
				localscores.append(-1);
		localscores.append(len(exons[i]))
		
		lmax = max(localscores);
		route = pathes[localscores.index(lmax)];
		
		maxes[i] = lmax;
		pathes[i] = [i]
		pathes[i].extend(route);
		
	beststop, stub = max(maxes.items(), key = lambda x: x[1])
	bestroute = pathes[beststop]
	
	return bestroute
	

def get_exons(interval, exons):
	'''Outputs all exons which may form a structure of the given interval(circle)
	
	interval pybedtools.Interval: genomic interval representing a circle (boundaries of the transcript)
	exons list: each element is a tuple representing genomic interval of an exon
	'''
	
	exons = [Interval(x[0], int(x[1]), int(x[2])) for x in set(exons)]
	exons.sort(key = lambda x: x.start)
	exons = filter(lambda x: int(x[1]) >= interval.start and int(x[2]) <= interval.end, exons)
	
	if(exons and interval.start == int(exons[0][1]) and  interval.end == int(exons[-1][2])):
		route = get_longest(exons)
		texons = [exons[x] for x in route[::-1]]
		size = len(texons);
		ls = ",".join([str(len(x)) for x in texons])
		starts = ",".join([str(x.start) for x in texons])
		print '%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%s,\t%s,' % (interval.chrom, interval.start, interval.end, interval.name, '0', interval.strand, interval.start, interval.end, '0', size, ls, starts) 
	else:
		print '%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%s,\t%s,' %(interval.chrom, interval.start, interval.end, interval.name, '0', interval.strand, interval.start, interval.end, '0', 1, len(interval), interval.start)
		
	return None


#def get_longest(exons):
	#maxes = defaultdict(float);
	#pathes = defaultdict(list);
	#for i in range(len(exons)):
		#localscores = [];
		#for j in range(i):
			#if(check_intersection(exons[i], exons[j])):
				#localscores.append[0];
			#else:
				#localscores.append(len(exons[i]) + maxes[j])



#get exons from annotation gff3 file
exons = [];
for interval in BedTool(args.gff3):
	if('ID' in interval.attrs and interval.attrs['ID'].split(':')[0] == 'gene'):
		curname = interval.attrs['gene_id']
		enames = set()
	if(interval[2] == 'exon'):
		if(interval.name not in enames):
			enames.add(interval.name)
			interval.name = curname
			exons.append(gff2bed(interval))
		

	
#Get an intersection between circles and expms
bed = BedTool(args.path);
offset  = bed.field_count();
intersection = bed.intersect(b=exons, s=True, wao=True);

curname = ''
cexons = []
for interval in intersection:
	if(curname == interval.name):
		cexons.append(tuple(interval[offset:offset+6]))
	else:
		if(curname):
			get_exons(cinterval, cexons)
		cinterval = interval
		cexons = [tuple(interval[offset:offset+6])]
		curname = interval.name
else:
	get_exons(cinterval, cexons)
		

	


