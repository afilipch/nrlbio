#! usr/bin/python
'''script generates html reports and plots for given bed/gff file'''
import os;
import sys;
import argparse;
from collections import defaultdict

from pybedtools import BedTool

from nrlbio.statistics import bed
from nrlbio import html
from nrlbio import annotation 




parser = argparse.ArgumentParser(description='script generates html reports and plots for given bed/gff file');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to bed/gff file"); 
parser.add_argument('-o', '--output', nargs = '?', default = None, type = str, help = "path to the output report. If not set, report will be generated in current folder, deriving file name from bed/gff file");
parser.add_argument('--attributes', nargs = '+', default = [], type = str, help = "set gff attributes to produce statistics of");
parser.add_argument('--annotation', nargs = '?', default = False, const = True, type = bool, help = "if True, process given files as chiflex/annotate_bed.py annotated file");
parser.add_argument('--draw', nargs = '?', default = False, const = True, type = bool, help = "if True, script will also generate plots for html report");
parser.add_argument('--sparse_coefficient', nargs = '?', default = 1, type = int, help = "If sparse_coefficient set to 3, then only each 3rd record will be analyzed, saves time in a case of large files");
args = parser.parse_args();

#set for everything for annotation

def annotation_generator(bedtool):
	for interval in bedtool:
		annotation.assign_stat_attributes(interval)		
		yield interval;

		
		

if(args.annotation):
	attributes = list(annotation.STAT_ATTRS);
	for a in args.attributes:
		if(a not in attributes):
			attributes.append(a);
	generator = annotation_generator(BedTool(args.path));
else:
	attributes=args.attributes
	generator = BedTool(args.path)
				
############################################################################################################################

#parse options
if(args.output):
	report_name = args.output;
else:
	report_name = ".".join(os.path.basename(args.path).split(".")[:-1] + ['html'])
############################################################################################################################



#sys.stderr.write("%s\n" % str(attributes))





bedstat = bed.Stat(name=os.path.splitext(os.path.basename(report_name))[0])
bedstat.fill_stat_sheet(generator, attributes=attributes, sparse_coefficient = args.sparse_coefficient)
if(args.draw):
	dirname = os.path.splitext(report_name)[0]
	try:
		os.mkdir(dirname)
	except:
		sys.stderr.write("directory %s already exists, plots will be drawn there\n" % os.path.abspath(dirname));
	bedstat.generate_plots(configuration='bedstat', output=dirname)