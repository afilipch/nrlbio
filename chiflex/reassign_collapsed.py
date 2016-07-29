'''Reassign nonuniquely mapped part of chimeras back in a way that closest(on a reference) hits are selected'''
import argparse
import sys;
from collections import defaultdict

from pybedtools import BedTool

from nrlbio.generators import generator_doublebed


parser = argparse.ArgumentParser(description='Reassign nonuniquely mapped part of chimeras back in a way that closest(on a reference) bed are selected');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the chimeras, bed format");
parser.add_argument('--collapsed', nargs = '?', required = True, type = str, help = "Path to collapsed nonunique mappings, bed file");
args = parser.parse_args();



pairs = {}
cdict = {}

for i1, i2 in generator_doublebed(args.path):
	pairs[(i1.chrom, i1.start, i1.end, i1.strand)] = (i2.chrom, i2.start, i2.end, i2.strand)
	cdict[(i1.chrom, i1.start, i1.end, i1.strand)] = [i1];
	cdict[(i2.chrom, i2.start, i2.end, i2.strand)] = [i2];




curname = '';
nonunique = [];
c = 0;
for interval in BedTool(args.collapsed):
	
	if(interval.name == curname):
		nonunique.append(interval);
	else:
		c += 1
		if(nonunique):
			fint = nonunique[0]
			l = cdict.get( (fint.chrom, fint.start, fint.end, fint.strand), None);
			if(l):
				l.extend(nonunique)
		curname = interval.name;
		nonunique = [interval];
		
	if(c % 1000 == 0):
		sys.stderr.write("%d processed\n" % c);
		
else:
	if(nonunique):
		fint = nonunique[0]
		l = cdict.get( (fint.chrom, fint.start, fint.end, fint.strand), None);
		if(l):
			l.extend(nonunique)
			
			
			
for k, v in cdict.items():
	if(len(v) > 1):
		print k
		for i in v:
			sys.stdout.write(str(i))
		print
		print
	
	
#for l in clist:
	#print "_"*130
	#for i in l:
		#sys.stdout.write(str(i))