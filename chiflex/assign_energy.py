#! /usr/lib/python
'''assignes hybridization energy to interactions'''

import argparse
import sys

from nrlbio.generators import generator_doublebed
from nrlbio.rnahybrid import get_rnahybrid, get_shuffled_rnahybrid


parser = argparse.ArgumentParser(description='assignes hybridization energy to interactions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to doublebed/gff file");
parser.add_argument('-st', '--shuffle_trials', nargs = '?', default = 0, type = int, help = "Controls how many times sequences will be shuffled to get background(control) hybridiztion energy. If set to 0, background energy will not be assigned");
parser.add_argument('-p', '--pattern', nargs = '?', default = False, const=True, type = bool, help = "If set, hybridization pattern will be assigned as well");
args = parser.parse_args();


for i1, i2 in generator_doublebed(args.path):
	energy, pattern = get_rnahybrid(i1.attrs['seq'], i2.attrs['seq']);
	i1.attrs['energy'] = str(energy);
	i2.attrs['energy'] = str(energy);
	if(args.pattern):
		i1.attrs['pattern'] = ",".join([str(x) for x in pattern]);
		i2.attrs['pattern'] = ",".join([str(x) for x in pattern]);
	
	if(args.shuffle_trials):
		energy_shuffled, pattern_shuffled = get_shuffled_rnahybrid(i1.attrs['seq'], i2.attrs['seq'], trials=args.shuffle_trials);
		i1.attrs['energy_shuffled'] = str(energy_shuffled);
		i2.attrs['energy_shuffled'] = str(energy_shuffled);
		if(args.pattern):
			i1.attrs['pattern_shuffled'] = ",".join(["%1.2f" % x for x in pattern_shuffled]);
			i2.attrs['pattern_shuffled'] = ",".join(["%1.2f" % x for x in pattern_shuffled]);
		
	sys.stdout.write(str(i1))
	sys.stdout.write(str(i2))