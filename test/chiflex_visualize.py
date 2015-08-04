#! /usr/lib/python
'''visualizes chiflex performance on artificially generated reads'''
import argparse
import os;
import sys;
from collections import defaultdict

import yaml

from nrlbio.pyplot_extension import pie, histogram



parser = argparse.ArgumentParser(description='visualizes chiflex performance on artificially generated reads');
#parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to initial reads. Fasta file");
parser.add_argument('-o', '--outdir', nargs = '?', default = 'evaluation_plots', type = str, help = "Name of the output directory. If directory exists, files will be (re)written there")
parser.add_argument('-ms', '--mapping_stat', nargs = '?', required = True, type = str, help = "Path to the mapping statistics in yaml format");
parser.add_argument('--chimera_stat', nargs = '?', required = True, type = str, help = "Path to the chimera mapping statistics in yaml format");
parser.add_argument('--single_stat', nargs = '?', required = True, type = str, help = "Path to the single mapping statistics in yaml format");
parser.add_argument('--control_stat', nargs = '?', required = True, type = str, help = "Path to the control mapping statistics in yaml format");
args = parser.parse_args();

mapped2real_dir =  os.path.join(args.outdir, "mapped2real")
real2mapped_dir =  os.path.join(args.outdir, "real2mapped")

if(not os.path.exists(args.outdir)):
    os.makedirs(args.outdir)
else:	
	sys.stderr.write("%s folder currently exists, files will be (re)written there\n" % os.path.abspath(args.outdir));
if(not os.path.exists(mapped2real_dir)):
	os.makedirs(mapped2real_dir)
if(not os.path.exists(real2mapped_dir)):
	os.makedirs(real2mapped_dir)
		
	
	
#________________________________________________________________________________________________________________
#get basic mapping statistics
#________________________________________________________________________________________________________________
#with open(args.mapping_stat, 'r') as f:
	#mapping_stat = yaml.load(f);
	
	
#real2mapped = defaultdict(lambda: defaultdict(int))
#mapped2real = defaultdict(lambda: defaultdict(int))
#for (read_type, mapped_type, is_exact, ttype), counts in mapping_stat.iteritems():
	#if(read_type==mapped_type):
		#if(is_exact):
			#mp = "%s(correct)" % mapped_type
		#else:
			#mp = "%s(incorrect)" % mapped_type
	#else:
		#mp = mapped_type;
		
	#real2mapped[read_type][mp] += counts;
	#real2mapped["%s(%s)" % (read_type, ttype)][mp] += counts;
	
	#mapped2real[mapped_type][read_type] += counts
	
	
#for title, data in real2mapped.items():
	#output = os.path.join(real2mapped_dir, title)
	#pie(data, top=10, min_fraction=0.05, title=title, output=output, colors=('yellowgreen', 'gold', 'lightskyblue', 'lightcoral', 'azure', 'seashell', 'darkorchid', 'chartreuse'), pctdistance=0.8, labeldistance=1.15)	
	
#for title, data in mapped2real.items():
	#output = os.path.join(mapped2real_dir, title)
	#pie(data, top=10, min_fraction=0.05, title=title, output=output, colors=('yellowgreen', 'gold', 'lightskyblue', 'lightcoral', 'azure', 'seashell', 'darkorchid', 'chartreuse'), pctdistance=0.8, labeldistance=1.15)	
	
	
	
	
	
#________________________________________________________________________________________________________________
#get detailed(length, coversions, gaps) mapping statistics on chimeras
#________________________________________________________________________________________________________________	
min_support = 1000	
minlength_range = (20, 52)
maxlength_range = (50, 80)
color = 'skyblue'
	
	
#________________________________________________________________________________________________________________
#get detailed(length, coversions, gaps) mapping statistics on chimeras
#________________________________________________________________________________________________________________
read_type = 'chimera'

with open(args.chimera_stat, 'r') as f:
	chimera_stat = yaml.load(f);
	

real2mapped_minlength = defaultdict(lambda: defaultdict(int))
mapped2real_minlength = defaultdict(lambda: defaultdict(int))
real2mapped_maxlength = defaultdict(lambda: defaultdict(int))
mapped2real_maxlength = defaultdict(lambda: defaultdict(int))

real2mapped_minconv = defaultdict(lambda: defaultdict(int))
mapped2real_minconv = defaultdict(lambda: defaultdict(int))
real2mapped_maxconv = defaultdict(lambda: defaultdict(int))
mapped2real_maxconv = defaultdict(lambda: defaultdict(int))

for (mapped_type, is_exact, length1, length2, conv1, conv2), counts in chimera_stat.iteritems():
	minlength = min(length1, length2)
	maxlength = max(length1, length2)
	minconv = min(conv1, conv2)
	maxconv = max(conv1, conv2)
	
	if(read_type==mapped_type):
		if(is_exact):
			mp = "%s(correct)" % mapped_type
		else:
			mp = "%s(incorrect)" % mapped_type
	else:
		mp = mapped_type;
		
	real2mapped_minlength[(read_type, mp)][minlength] += counts;
	real2mapped_maxlength[(read_type, mp)][maxlength] += counts;
	
	real2mapped_minconv[(read_type, mp)][minconv] += counts;
	real2mapped_maxconv[(read_type, mp)][maxconv] += counts;

	mapped2real_minlength[(mapped_type, read_type)][minlength] +=counts
	mapped2real_maxlength[(mapped_type, read_type)][maxlength] +=counts
		
	mapped2real_minconv[(mapped_type, read_type)][minconv] +=counts
	mapped2real_maxconv[(mapped_type, read_type)][maxconv] +=counts
	


real2mapped_detailed_dir =  os.path.join(real2mapped_dir, read_type)
if(not os.path.exists(real2mapped_detailed_dir)):
	os.makedirs(real2mapped_detailed_dir)
	
mapped2real_detailed_dir =  os.path.join(mapped2real_dir, read_type)
if(not os.path.exists(mapped2real_detailed_dir)):
	os.makedirs(mapped2real_detailed_dir);
	
#plotting		
for title, data in real2mapped_minlength.items():
	if(sum(data.values())>=min_support):
		title = "%s(minlength)" % "->".join(title)
		output = os.path.join(real2mapped_detailed_dir, title)
		histogram(data, title=title, xlabel="length, nt", ylabel='num of reads', output=output, step = 1, range=minlength_range, align='mid', color=color)

for title, data in real2mapped_maxlength.items():
	if(sum(data.values())>=min_support):
		title = "%s(maxlength)" % "->".join(title)
		output = os.path.join(real2mapped_detailed_dir, title)
		histogram(data, title=title, xlabel="length, nt", ylabel='num of reads', output=output, step = 1, range=maxlength_range, align='mid', color=color)
	
	
for title, data in mapped2real_minlength.items():
	if(sum(data.values())>=min_support):
		title = "%s(minlength)" % "->".join(title)
		output = os.path.join(mapped2real_detailed_dir, title)
		histogram(data, title=title, xlabel="length, nt", ylabel='num of reads', output=output, step = 1, range=minlength_range, align='mid', color=color)

for title, data in mapped2real_maxlength.items():
	if(sum(data.values())>=min_support):
		title = "%s(maxlength)" % "->".join(title);
		output = os.path.join(mapped2real_detailed_dir, title)
		histogram(data, title=title, xlabel="length, nt", ylabel='num of reads', output=output, step = 1, range=maxlength_range, align='mid', color=color)
		
		
for title, data in real2mapped_minconv.items():
	if(sum(data.values())>=min_support):
		title = "%s(minconv)" % "->".join(title)
		output = os.path.join(real2mapped_detailed_dir, title)
		histogram(data, title=title, xlabel="num of conversions", ylabel='num of reads', output=output, step = 1, align='mid', color=color)

for title, data in real2mapped_maxconv.items():
	if(sum(data.values())>=min_support):
		title = "%s(maxconv)" % "->".join(title);
		output = os.path.join(real2mapped_detailed_dir, title)
		histogram(data, title=title, xlabel="num of conversions", ylabel='num of reads', output=output, step = 1, align='mid', color=color)	
		
		
for title, data in mapped2real_minconv.items():
	if(sum(data.values())>=min_support):
		title = "%s(minconv)" % "->".join(title)
		output = os.path.join(mapped2real_detailed_dir, title)
		histogram(data, title=title, xlabel="num of conversions", ylabel='num of reads', output=output, step = 1, align='mid', color=color)

for title, data in mapped2real_maxconv.items():
	if(sum(data.values())>=min_support):
		title = "%s(maxconv)" % "->".join(title);
		output = os.path.join(mapped2real_detailed_dir, title)
		histogram(data, title=title, xlabel="num of conversions", ylabel='num of reads', output=output, step = 1, align='mid', color=color)

		
		
		
#________________________________________________________________________________________________________________
#get detailed(length, coversions, gaps) mapping statistics on chimeras
#________________________________________________________________________________________________________________		
read_type = 'chimera'

with open(args.chimera_stat, 'r') as f:
	chimera_stat = yaml.load(f);
	

real2mapped_length = defaultdict(lambda: defaultdict(int))
mapped2real_length = defaultdict(lambda: defaultdict(int))

real2mapped_conv = defaultdict(lambda: defaultdict(int))
mapped2real_conv = defaultdict(lambda: defaultdict(int))


for (mapped_type, is_exact, length1, length2, conv1, conv2), counts in chimera_stat.iteritems():
	minlength = min(length1, length2)
	maxlength = max(length1, length2)
	minconv = min(conv1, conv2)
	maxconv = max(conv1, conv2)
	
	if(read_type==mapped_type):
		if(is_exact):
			mp = "%s(correct)" % mapped_type
		else:
			mp = "%s(incorrect)" % mapped_type
	else:
		mp = mapped_type;
		
	real2mapped_minlength[(read_type, mp)][minlength] += counts;
	real2mapped_maxlength[(read_type, mp)][maxlength] += counts;
	
	real2mapped_minconv[(read_type, mp)][minconv] += counts;
	real2mapped_maxconv[(read_type, mp)][maxconv] += counts;

	mapped2real_minlength[(mapped_type, read_type)][minlength] +=counts
	mapped2real_maxlength[(mapped_type, read_type)][maxlength] +=counts
		
	mapped2real_minconv[(mapped_type, read_type)][minconv] +=counts
	mapped2real_maxconv[(mapped_type, read_type)][maxconv] +=counts	
	
	
	
	
	
	
	
	
	
	
		
#print "\n".join(sorted(real2mapped.keys()));

#print
#for k, v in real2mapped['chimera(intergenic|intergenic)'].items():
	#print "%s\t%d" % (k, v)
	
			