#! usr/bin/python
'''this script generates statistics in html format for provided sam/bam files'''
import os;
import sys;
import argparse;

import pysam

import sam
from nrlbio import html




parser = argparse.ArgumentParser(description='script generates html report and plots for given sam/bam files all together');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "paths to sam/bam files"); 
parser.add_argument('-o', '--output', nargs = '?', required = True, type = str, help = "paths to the output report");
parser.add_argument('--draw', nargs = '?', default = False, const = True, type = bool, help = "if True, scripts also generates plots for html report");
parser.add_argument('--ordered', nargs = '?', default = False, const = True, type = int, help = "If true sorts attributes by value, otherwise by key");
parser.add_argument('--top_entries', nargs = '?', default = 20, type = int, help = "Reports only top N entries for each statistics. Ignored if '--ordered'=False");
parser.add_argument('--sparse_coefficient', nargs = '?', default = 1, type = int, help = "If sparse_coefficient set to 3, then only each 3rd record will be analyzed, saves time in a case of large files");
parser.add_argument('--short_reference', nargs = '?', const = True, default = False, type = str, help = "If set, provides additional statistics. Useful for sam files derived from mapping to short reference(miRNAs, piRNAs, enhancers, etc.)")
args = parser.parse_args();


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


data = [];
for path in args.path:
	samfile = pysam.Samfile(path)
	samstat = sam.Stat(name=os.path.splitext(os.path.basename(path))[0]);	
	samstat.fill_stat_sheet(samfile.fetch(until_eof=True), short_reference = args.short_reference, detailed = True, sparse_coefficient = args.sparse_coefficient)
	data.append(samstat);
	
	
#htmlstat = html.Stat.from_stat_object(name, samstat, attributes, headers=headers, ordered=args.ordered, top_entries=args.top_entries);
#htmlstat.report(output=name)
if(args.draw):
	dirname = os.path.splitext(args.output)[0]
	try:
		os.mkdir(dirname)
	except:
		sys.stderr.write("directory %s already exists, plots will be drawn there" % os.path.abspath(dirname));
	sam.Stat.generate_multihist(data, configuration='samstat', output=dirname)