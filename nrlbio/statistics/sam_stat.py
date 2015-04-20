#! usr/bin/python
'''this script generates statistics in html format for provided sam/bam files'''
import os;
import sys;
import argparse;

import pysam

import sam
from nrlbio import html




parser = argparse.ArgumentParser(description='script generates html reports and plots for each given sam/bam file separately');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "paths to sam/bam files"); 
parser.add_argument('-o', '--output', nargs = '+', default = [], type = str, help = "paths to the output reports. If not set, reports will be generated in current folder, deriving file names from sam/bam files");
parser.add_argument('--draw', nargs = '?', default = False, const = True, type = bool, help = "if True, scripts also generate plots for html report");
parser.add_argument('--ordered', nargs = '?', default = False, const = True, type = int, help = "If true sorts attributes by value, otherwise by key");
parser.add_argument('--top_entries', nargs = '?', default = 20, type = int, help = "Reports only top N entries for each statistics. Ignored if '--ordered'=False");
parser.add_argument('--sparse_coefficient', nargs = '?', default = 1, type = int, help = "If sparse_coefficient set to 3, then only each 3rd record will be analyzed, saves time in a case of large files");
parser.add_argument('--short_reference', nargs = '?', const = True, default = False, type = str, help = "If set, provides additional statistics. Useful for sam files derived from mapping to short reference(miRNAs, piRNAs, enhancers, etc.)")
args = parser.parse_args();

#parse options
if(args.output):
	if(len(args.output) == len(args.path)):
		report_names = args.output;
	else:
		sys.exit('number of output files have to be equal to the number of input ones, or \'--output\' should not be set\n')
else:
	report_names = [".".join(os.path.basename(x).split(".")[:-1] + ['html']) for x in args.path]
############################################################################################################################


#set attributes names and their respective headers
attributes = ["ascore", "query_start", "clipped_length_right", "query_end", "conv", "conv_weighted", "conv_number",	"clipped_seq_left", "clipped_seq_right"]
		
headers = ["Alignment Score", "Start position of the match in query", "number of nucleotides soft clipped downstream", "End position of the match in query", "Type of conversion", "Type of conversion weighted to a number of conversions in one read", "Number of conversion per read", "Soft clipped upstream sequence", "Soft clipped downstream sequence"]

if(args.short_reference):
	attributes.append("query_ref_start")
	attributes.append("ref_start")
	attributes.append("ref_end")
	headers.append("Start position of the match in the query and the reference")
	headers.append("Start position of the match in the reference")
	headers.append("End position of the match in the reference")

############################################################################################################################



for path, name in zip(args.path, report_names):
	samfile = pysam.Samfile(path)
	samstat = sam.Stat(name=os.path.splitext(os.path.basename(name))[0]);	
	samstat.fill_stat_sheet(samfile.fetch(until_eof=True), short_reference = args.short_reference, detailed = True, sparse_coefficient = args.sparse_coefficient)
	htmlstat = html.Stat.from_stat_object(name, samstat, attributes, headers=headers, ordered=args.ordered, top_entries=args.top_entries);
	htmlstat.report(output=name)
	if(args.draw):
		dirname = os.path.splitext(name)[0]
		try:
			os.mkdir(dirname)
		except:
			sys.stderr.write("directory %s already exists, plots will be drawn there" % os.path.abspath(dirname));
		samstat.generate_hist(configuration='samstat', output=dirname)