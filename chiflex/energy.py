#! /usr/lib/python
'''assignes hybridization energy to interactions'''

import argparse
import sys
import random 
from collections import defaultdict


from nrlbio.generators import generator_doublebed
from nrlbio.rnahybrid import get_rnahybrid
from nrlbio.numerictools import cdf
from nrlbio.sequencetools import shuffle_string




parser = argparse.ArgumentParser(description='assignes hybridization energy to interactions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to doublebed/gff file");
parser.add_argument('-r', '--report', nargs = '?', default = "reports", type = str, help = "path to the report folder");
args = parser.parse_args();

def interactions2energy(interaction_seqs):
	return [get_rnahybrid(x[0], x[1])[0] for x in interaction_seqs]

#real interactions	
interaction_seqs = [(x[0].attrs['seq'], x[1].attrs['seq']) for x in generator_doublebed(args.path)];
e_signal = cdf(interactions2energy(interaction_seqs))

#shuffled pairs
shuffled_pairs = zip(random.sample([x[0] for x in interaction_seqs], len(interaction_seqs)), random.sample([x[1] for x in interaction_seqs], len(interaction_seqs)))
e_shufpairs = cdf(interactions2energy(shuffled_pairs))

#shuffled sequences
shuffled_seqs = [(shuffle_string(x[0]), shuffle_string(x[1])) for x in interaction_seqs]
e_shufseq = cdf(interactions2energy(shuffled_seqs))

#for kv in e_signal:
	#print "%1.1f\t%d" % kv;
	
	
	
#e_shufseq = defaultdict(int)
#e_shufpairs = defaultdict(int)