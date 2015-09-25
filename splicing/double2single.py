#! /usr/bin/python
'''Converts doublebed formatted linear/circular splice junctions into single bed format''' 
import sys;
import argparse

from nrlbio.generators import generator_doublebed;
from nrlbio.pybedtools_extension import construct_gff_interval


parser = argparse.ArgumentParser(description='Converts doublebed formatted linear/circular splice junctions into single bed format');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the chimeras, double bed/gff file");
parser.add_argument('--jtype', nargs = '?', required = True, choices=['csj', 'lsj'], type = str, help = "type of junction [csj|lsj]")
args = parser.parse_args()



def linear(i1, i2):
	global bu
	start = min(i1.stop, i2.stop)
	stop = max(i1.start, i2.start)
	if(stop<=start):
		#sys.stderr.write("%s%s%s\n" % (i1, i2, "_"*140));
		#bu+=1;
		return None
	score = str(int(i1.score) + int(i2.score));
	attrs = [('score1', i1.score), ('score2', i2.score), ('ID', i1.name.split("|")[0])]
	for k1, v1 in i1.attrs.items():
		if(v1==i2.attrs.get(k1)):
			attrs.append((k1,v1));
		
	
	return construct_gff_interval(chrom=i1.chrom, start=start, stop=stop, feature='lsj', score=score, strand=i1.strand, source='.', frame='.', attrs=attrs)


def circular(i1, i2):
	start = min(i1.start, i2.start)
	stop = max(i1.stop, i2.stop)
	score = str(int(i1.score) + int(i2.score));
	attrs = [('score1', i1.score), ('score2', i2.score), ('ID', i1.name.split("|")[0])]
	for k1, v1 in i1.attrs.items():
		if(v1==i2.attrs.get(k1)):
			attrs.append((k1,v1));
	
	return construct_gff_interval(chrom=i1.chrom, start=start, stop=stop, feature='csj', score=score, strand=i1.strand, source='.', frame='.', attrs=attrs)



if(args.jtype=='lsj'):
	func = linear;
else:
	func = circular;

overlaping = 0;
for i1, i2 in generator_doublebed(args.path):
	ni = func(i1, i2)
	if(ni):
		sys.stdout.write(str(ni));
	else:
		overlaping+=1
	
sys.stderr.write("\nNumber of removed overlaping intervals %d\n" % overlaping);	