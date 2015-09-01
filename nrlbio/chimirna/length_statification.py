#! /usr/bin/python	
'''Script produce different chipart.tsv files based on the length of the right part of chimeric read'''
from collections import *;
import parsing;
import chipart_lib;
import argparse;
import logging
import os;
import sys;
from time import strftime, gmtime

logger = logging.getLogger(__name__);
logger.setLevel(logging.DEBUG)
fh = logging.FileHandler(os.path.join("log", "length_stratification.txt"), 'a')
sh = logging.StreamHandler()
logger.addHandler(fh);
logger.addHandler(sh);

wf = open(os.path.join("log", "workflow.txt"), 'a');
wf.write("python " + " ".join(sys.argv) + "\n\n");
wf.close


parser = argparse.ArgumentParser(description='Script produce different chipart.tsv files based on the length of the right part of chimeric read');
# input files
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to anchors");
parser.add_argument('-l', '--length', nargs = '?', default = 12, type = int, help = "minimum length of the right chipart to produce right part further")
parser.add_argument('-p', '--postfix', nargs = '?', default = '', const = "_control", type = str, help = "postfix serves to tag resulted files")
args = parser.parse_args();
sboundary = 4;


l = open(os.path.join("left", 'long' + args.postfix + '.tsv'), 'w');
s = open(os.path.join("left", 'short' + args.postfix + '.tsv'), 'w');
n = open(os.path.join("left", 'null' + args.postfix + '.tsv'), 'w');

h = open(args.path[0]);

lc, sc, nc, uc = 0,0,0,0

for line in h:
	arr = line.split("\t")
	diff  = int(arr[4]) - int(arr[3]);
	if(diff >= args.length):
		lc += 1;
		l.write(line);
	elif(sboundary < diff < args.length):
		s.write(line);
		sc += 1;
	elif(diff == 0):
		n.write(line);
		nc += 1;
	else:
		uc += 1;
		pass;
		
h.close();
l.close();
s.close();
n.close();


logger.info(strftime("time of analysis %Y-%m-%d %H:%M:%S", gmtime()))	
logger.info("boundary between rightparts short rightparts: %d, boundary very short and short rightparts: %d" % (args.length, sboundary))	
logger.info("%d long rightparts, %d short rightparts, %d very short rightparts, %d none rightparts" % (lc, sc, nc, uc))	

