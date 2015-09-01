import parsing;
#! /usr/bin/python	
import sys;
import os;
import copy;
import argparse;

parser = argparse.ArgumentParser(description='create artificial genome on basis of profided fasta file');
parser.add_argument('path', metavar = 'N', nargs = 1, type = str, help = "Path to fasta file");
parser.add_argument('--From', nargs = '?', type = str, help = "conversion from");
parser.add_argument('--to',   nargs = '?', type = str, help = "conversion to");
parser.add_argument('-f', '--format', nargs = '?', type = str, required = True, help = "type of given file, fa or qfa");
args = parser.parse_args();

if(args.format == "fa"):
	f = open(args.path[0]);  
	for line in f:
		line = line.strip();
		if (line[0] == ">"):
			cur_key = line;
		else:  
			print cur_key + "_civ_-1";
			print line;
			for i in range(0, len(line)):
				if(line[i] == args.From):
					print cur_key + "_civ_"  + str(i);
					print line[:i] + args.to + line [i+1:];
	f.close();  
	  
elif(args.format == "qfa"):
	f = open(args.path[0], "r", 1);
	line_set = [];
	for line in f:
		line = line.strip();
		line_set.append(line);
		if (len(line_set) == 4):
			print line_set[0] + "_civ_-1";
			print line_set[1];
			print line_set[2];
			print line_set[3];
			for i in range(0, len(line_set[1])):
				if(line_set[1][i] == args.From):
					print line_set[0] + "_civ_" + str(i);
					print line_set[1][:i] + args.to + line_set[1][i+1:];
					print line_set[2];
					print line_set[3];
			line_set = [];
	f.close()	
else:
  raise Exception("please provide format of file as fa or qfa"); 
