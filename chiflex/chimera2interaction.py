'''converts bed-like file of chimeric reads into interactions. That is merging chimeras with intersecting regions'''

import argparse;

import pybedtools

from nrlbio.interaction import bed2interactions

parser = argparse.ArgumentParser(description='converts bed-like file of chimeric reads into interactions. That is merging chimeras with intersecting regions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to chimeras bed-like file");
parser.add_argument('-d', '--distance', nargs = 2, default = [12, 12], type = int, help = "minimum overlaps in nucleotides to merge intervals");
parser.add_argument('-n', '--name', nargs = '?', required = True, type = str, help = "name for interactions")
args = parser.parse_args();


bed = pybedtools.BedTool(args.path);
for inter in bed2interactions(bed, args.distance, name=args.name):
	print inter