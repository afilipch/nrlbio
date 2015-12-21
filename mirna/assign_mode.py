#! /usr/lib/python
'''Assignes binding mode for mirna:target interaction'''

import argparse
import os
import sys
from collections import defaultdict

from nrlbio.mirna import fasta2mirnas, MODES_ORDER
from nrlbio.generators import generator_doublebed
from nrlbio.sequencetools import shuffle_string


parser = argparse.ArgumentParser(description='Assignes binding mode for mirna:target interaction');	
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the mirna:target interactions, doublegff format");
parser.add_argument('--mir', nargs = '?', required=True, type = str, help = "path to the miRNAs, fasta format");
parser.add_argument('-sc', '--set_control', nargs = '?', default = False, const=True, type = int, help = "If set, background(shuffled target sequence) seed match mode will be assigned");
args = parser.parse_args();

stat = defaultdict(lambda: defaultdict(int));

def get_mode(modes):
	for mode in MODES_ORDER:
		if(modes[mode]):
			return mode
	else:
		return 'none'

#get dictionary of mirna.Mirna objects
mirdict = fasta2mirnas(args.mir);

for i1, i2 in generator_doublebed(args.path):
	mode = get_mode(mirdict[i1.chrom].find_fixed_types(i2.attrs['seq']))
	i1.attrs['mode'] = mode;
	i2.attrs['mode'] = mode;
	stat[mode]['signal']+=1
	
	if(args.set_control):
		mode_shuffled = get_mode(mirdict[i1.chrom].find_fixed_types(shuffle_string(i2.attrs['seq'])))
		i1.attrs['mode_shuffled'] = mode_shuffled;
		i2.attrs['mode_shuffled'] = mode_shuffled;
		stat[mode_shuffled]['control']+=1
		
	sys.stdout.write(str(i1))
	sys.stdout.write(str(i2))
	
total = float(sum([x['signal'] for x in stat.values()]))
sys.stderr.write("mode\tsignal\tsignal_fraction\tcontrol\tcontrol_fraction\n")
for mode, td in stat.items():
	sys.stderr.write("%s\t%d\t%1.2f\t%d\t%1.2f\n" % (mode, td['signal'], td['signal']/total, td['control'], td['control']/total))
	