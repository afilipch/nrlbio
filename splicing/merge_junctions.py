#! /usr/bin/python
'''Collapses circular or linear splice junctions with exactly the same coordinates into a single entry''' 
import sys;
import argparse
from collections import defaultdict, Counter

from pybedtools import BedTool

from nrlbio.pybedtools_extension import construct_gff_interval


parser = argparse.ArgumentParser(description='Collapses circular or linear splice junctions with exactly the same coordinates into a single entry');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the file with circular or linear splice junctions, bed/gff file");
parser.add_argument('-od', '--dictionary', nargs = '?', required = True, type = str, help = "path to output \"interaction to read id\" file")
parser.add_argument('--jtype', nargs = '?', required = True, choices=['csj', 'lsj'], type = str, help = "type of junction [csj|lsj]")
parser.add_argument('--support', nargs = '?', default = 1, type = int, help = "min read support for a splice junction to pass")
args = parser.parse_args()


def merge(reads, jnum, jtype):
	score = str(max([int(x.score) for x in reads]));
	gaps = [x.attrs['gap'] for x in reads]
	gap = max(set(gaps), key=gaps.count)
	#sys.stderr.write("%s\t%s\n" % (gaps, gap));
	attrs = [('n_uniq', len(reads)), ('gap', gap), ("ID", "c%d" % jnum)];
	if('ktype' in reads[0].attrs):
		attrs.append(('ktype', reads[0].attrs['ktype']))
	#for attr_name in reads[0].attrs.keys():
		#values = set(x.attrs.get(attr_name) for x in reads);
		#if(len(values)==1 and attr_name!='gap'):
			#attrs.append((attr_name, list(values)[0]));
	return construct_gff_interval(reads[0].chrom, reads[0].start, reads[0].stop, feature=jtype, score=score, strand=reads[0].strand, source='.', frame='.', attrs=attrs)


jnum = 0;
passed = 0;
total = 0;
junctions = {};
current = [];
with open(args.dictionary, 'w') as fdict:
	for interval in BedTool(args.path):
		total+=1
		location = (interval.chrom, interval.start, interval.stop);
		if(location in junctions):
			junctions[location] += 1;
			current.append(interval);
		elif(len(current)>=args.support):
			passed+=len(current);
			jnum+=1
			sys.stdout.write(str(merge(current, jnum, args.jtype)));
			for read in current:
				fdict.write("c%d\t%s\n" % (jnum, read.name));			
			junctions[location] = 1;
			current = [interval]
		else:
			junctions[location] = 1;
			current = [interval]
	else:
		if(len(current)>args.support):
			passed+=len(current);
			jnum+=1
			sys.stdout.write(str(merge(current, jnum, args.jtype)));
			for read in current:
				fdict.write("c%d\t%s\n" % (jnum, read.name));

				
				
				
sys.stderr.write("total splicing events: %d\npassed splicing events: %d\nfraction passed %1.5f\n\n" % (total, passed, float(passed)/total));
sys.stderr.write("%d splice sites generated from %d splicing events\n\n" % (jnum, passed));

sys.stderr.write("\nsupport\tnum of splice junctions\n")
support_count = Counter(junctions.values())
top = 10;
for k in range(1, top+1):
	sys.stderr.write("%d\t%d\n" % (k, support_count[k]));
sys.stderr.write( ">%d\t%d\n" % (top, sum([x[1] for x in support_count.items() if x[0]>top])) );
