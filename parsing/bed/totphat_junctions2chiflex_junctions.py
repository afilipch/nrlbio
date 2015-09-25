#! /usr/bin/python
'''Converts tophat output (junctions.bed) to the format campatible with the format of chiflex linear splice junctions file(splicing/collapsed.lsj.gff). Also filters tophat output on basis of read support number. Format of the intervals coordinates after confersion should be:
start: left edge of junction(at splice signal)
end: right edge of junction(at splice signal)
'''

import argparse
import sys;

from pybedtools import BedTool;



parser = argparse.ArgumentParser(description='Converts tophat output (junctions.bed) to the format campatible with the format of chiflex linear splice junctions file(splicing/collapsed.lsj.gff). Also filters tophat output on basis of read support number. Format of the intervals coordinates after confersion should be:\nstart: left edge of junction(at splice signal)\nend: right edge of junction(at splice signal)', formatter_class = argparse.RawTextHelpFormatter);

parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the tophat junctions, bed12");
parser.add_argument('--support', nargs = '?', default = 2, type = int, help = "min read support for a splice junction to pass")
args = parser.parse_args();

for interval in BedTool(args.path):
	if(int(interval.score)>=args.support):
		block_sizes = [int(x) for x in interval[10].split(",")];
		interval.start = interval.start + block_sizes[0] 
		interval.end = interval.end - block_sizes[1]
		sys.stdout.write(str(interval));