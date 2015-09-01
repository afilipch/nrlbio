from collections import *; 
import parsing;
import sys;
import re;
import os;
import time;
import csv;
import copy;
import argparse;
import multiprocessing;
from generators import *;


parser = argparse.ArgumentParser(description='remove adapter parts from reads');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to fastq file of reads");
parser.add_argument('-a', '--adapter', nargs = '?',  required = True, help = "path to adapters in fasta, or raw sequence");
parser.add_argument('-l', '--length', nargs = '?', default = 15, type = int, help = "min length of resulted reads");
parser.add_argument('-b', '--backward', nargs = '?', default = False, const = True, type = str, help = "search for 5'adapter");
args = parser.parse_args();

def remove_single(line_set, adapters, length, backward):
	found, removed = 0, 0
	for k, v in adapters.iteritems():
		pos =  line_set[1].find(v)
		if(pos != -1):
			found = 1;
			if(pos >= length):
				t = [line_set[0], line_set[1][:pos], line_set[2], line_set[3][:pos]];
			else:
				t = [];
				removed = 1;
			break;	
	else:	
		t = list(line_set);
	if(backward and t):
		t[1] = t[1][::-1];
		t[3] = t[3][::-1];
	return t, found, removed


def remove(arr):
	seq_list, adapters, length, backward = arr;
	found, removed = 0, 0
	ans = [];
	for line_set in seq_list:
		t,f,r = remove_single(line_set, adapters, length, backward)
		ans.append(t);
		found += f;
		removed += r;
	return ans, found, removed;	
				

try:
	adapters = parsing.fasta2dict(args.adapter, reverse = args.backward);
except:
	s = args.adapter.upper()
	if(args.backward):
		s = s[::-1]
	if (set(s) <= set('ACTGU')):
		adapters = {"cmd": s}
	else:	
		raise Exception("neither valid file name nor raw sequence")

#sys.stderr.write("%s\n" % adapters)		
  
#g =  grouper(generator_fastq(args.path, ["id", "seq", "sign", "qual"], reverse = args.backward), 10000)
#print g.next()


found, removed = 0, 0
for line_set in generator_fastq(args.path, ["id", "seq", "sign", "qual"], reverse = args.backward):
	t, f, r = remove_single(line_set, adapters, args.length, args.backward)
	found += f;
	removed += r;
	if(t):
		print "\n".join(t);

###multiprocessing versions works slower
#pool = multiprocessing.Pool(processes = 8)
#res = pool.imap(remove, ((seq_list, adapters, args.length, args.backward) for seq_list in grouper(generator_fastq(args.path, ["id", "seq", "sign", "qual"], reverse = args.backward), 100000)) )   

#found, removed = 0, 0
#for my_list, f, r in res:
	#found += f;
	#removed += r;
	#for el in my_list:
		#if(el):
			#print "\n".join(el);



sys.stderr.write("%s\t%d\n" % ("adapters found", found));
sys.stderr.write("%s\t%d\n" % ("short reads discarded", removed));