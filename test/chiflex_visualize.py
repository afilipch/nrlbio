#! /usr/lib/python
'''visualizes chiflex performance on artificially generated reads'''
import argparse
import os;
import sys;
from collections import defaultdict

import yaml

from nrlbio.pyplot_extension import pie



parser = argparse.ArgumentParser(description='visualizes chiflex performance on artificially generated reads');
#parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to initial reads. Fasta file");
parser.add_argument('-o', '--outdir', nargs = '?', default = 'chiflex_evaluation', type = str, help = "Name of the output directory. If directory exists, files will be (re)written there")
parser.add_argument('-ms', '--mapping_stat', nargs = '?', required = True, type = str, help = "Path to the mapping statistics in yaml format");
args = parser.parse_args();

mapped2real_dir =  os.path.join(args.outdir, "mapped2real")
real2mapped_dir =  os.path.join(args.outdir, "real2mapped")

if(not os.path.exists(args.outdir)):
    os.makedirs(args.outdir, mapped2real_dir, real2mapped_dir)
else:	
	sys.stderr.write("%s folder currently exists, files will be (re)written there\n" % os.path.abspath(args.outdir));
if(not os.path.exists(mapped2real_dir)):
	os.makedirs(mapped2real_dir)
if(not os.path.exists(real2mapped_dir)):
	os.makedirs(real2mapped_dir)
	
	

#get mapping statistics
with open(args.mapping_stat, 'r') as f:
	mapping_stat = yaml.load(f);
	
	
	
real2mapped = defaultdict(lambda: defaultdict(int))
mapped2real = defaultdict(lambda: defaultdict(int))
for (read_type, mapped_type, is_exact, ttype), counts in mapping_stat.iteritems():
	if(read_type==mapped_type):
		if(is_exact):
			mp = "%s(correct)" % mapped_type
		else:
			mp = "%s(incorrect)" % mapped_type
	else:
		mp = mapped_type;
		
	real2mapped[read_type][mp] += counts;
	real2mapped["%s(%s)" % (read_type, ttype)][mp] += counts;
	
	mapped2real[mapped_type][read_type] += counts
	
	
	
for title, data in real2mapped.items():
	output = os.path.join(real2mapped_dir, title)
	pie(data, top=10, min_fraction=0.05, title=title, output=output, colors=('yellowgreen', 'gold', 'lightskyblue', 'lightcoral', 'azure', 'seashell', 'darkorchid', 'chartreuse'), pctdistance=0.8, labeldistance=1.15)	
	
for title, data in mapped2real.items():
	output = os.path.join(mapped2real_dir, title)
	pie(data, top=10, min_fraction=0.05, title=title, output=output, colors=('yellowgreen', 'gold', 'lightskyblue', 'lightcoral', 'azure', 'seashell', 'darkorchid', 'chartreuse'), pctdistance=0.8, labeldistance=1.15)	
	
	
	
#print "\n".join(sorted(real2mapped.keys()));

#print
#for k, v in real2mapped['chimera(intergenic|intergenic)'].items():
	#print "%s\t%d" % (k, v)
	
			