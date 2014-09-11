# /usr/bin/python
'''collections of classes and functions to deal with sam/bam files'''

import sys;
from collections import namedtuple;

from nrlbio import numerictools;


def key_alignment_score(arw):
	return arw.AS;


class ArWrapper(object):
	def __init__(self, aligned_read, rname):
		self.aligned_read = aligned_read;
		self.qname = aligned_read.qname;
		self.rname = rname
		self.AS = aligned_read.opt("AS")
		if(rname.split("_")[0] == "random"):
			self.control = True;
		else:
			self.control = False;


def demultiplex_read_hits(arwlist, key_function):
	real = filter(lambda x: not x.control, arwlist);
	control = filter(lambda x: x.control, arwlist);
	best_real, max_real = numerictools.maxes(real, key_function)
	best_control, max_control = numerictools.maxes(control, key_function)
	if(max_real > max_control):
		if(len(best_real) == 1):
			return best_real[0], [], None
		else:
			return None, best_real, None
	elif(max_control > max_real):
		return None, best_real, best_control[0];
	else:
		return None, best_real, None
	
