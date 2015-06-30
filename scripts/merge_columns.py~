#! /usr/bin/python
import sys;
import argparse;


parser = argparse.ArgumentParser(description='merge collumns for different files. NOTE: order in the output is defined primary by order of files in path and secondary by order of comma delimeted numbers in [--columns](1-based)\nFor example: call {merge_columns.py file1 file2 --collumns 3,1 3,6} will create a file with a collumns: 3rd from the file1, 1st from the file1, 3rd from the file2, 6th from the file2');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to delimeted files");
parser.add_argument('-c', '--collumns', nargs = '+', required = True, help = "positions of collumns to output, for one file one set of comma delimeted numbers. NOTE: it is 1-based attribute");
parser.add_argument('-d', '--delimeter', nargs = '?', default = "\t", type = str, help = "delimeter of collumns in files");
args = parser.parse_args();

collumns = [map(lambda m: int(m)-1, x.split(",")) for x in args.collumns]; #list of tuples, first layer correspond to files, second to positions of collumns
min_col = min([min(x) for x in collumns])
if(min_col < 0):
	sys.exit("Collums numbers have to be positive integers. %d provided instead\n" % (min_col+1))

handlers = [open(path) for path in args.path]
delim = args.delimeter


for lines in zip(*handlers):
	a = [];
	for l,col in zip(lines, collumns):
		la=l.strip().split(delim)
		a.extend([la[x] for x in col]);
	print delim.join(a);
  

for handler in handlers:
	handler.close();