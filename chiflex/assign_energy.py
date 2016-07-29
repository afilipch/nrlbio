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
parser.add_argument('-si', '--shuffle_interactions', nargs = '?', default = False, const=True, type = bool, help = "If set, shuffling will be performed for miRNA tagrget pairs");
args = parser.parse_args();	


if(args.shuffle_interactions):
	import random
	mirnas = [];
	targets = [];
	
	for i1, i2 in generator_doublebed(args.path):
		mirnas.append(i1);
		targets.append(i2);
		
	for i1, i2 in zip(mirnas, targets):
		energy, pattern = get_rnahybrid(i2.attrs['seq'], i1.attrs['seq']);
		i1.attrs['energy'] = str(energy);
		i2.attrs['energy'] = str(energy);
		i1.attrs['pattern'] = ",".join([str(x) for x in pattern]);
		i2.attrs['pattern'] = ",".join([str(x) for x in pattern]);
		
		if(args.shuffle_trials):
			energy_shuffled, patterns  = 0, []
			for _ in range(args.shuffle_trials):
				
				for __ in range(100):
					index = random.randint(0, len(targets)-1)
					target = targets[index]
					if(mirnas[index].name != i1.name):
						break;
						
				energy, pattern = get_rnahybrid(target.attrs['seq'], i1.attrs['seq']);
				energy_shuffled += energy;
				patterns.append(pattern)
				
			energy_shuffled = energy_shuffled/float(args.shuffle_trials);
			pattern_shuffled = [];
			for p in zip(*patterns):
				pattern_shuffled.append(sum(p)/float(args.shuffle_trials))
				
			i1.attrs['energy_shuffled'] = str(energy_shuffled);
			i2.attrs['energy_shuffled'] = str(energy_shuffled);
			i1.attrs['pattern_shuffled'] = ",".join(["%1.2f" % x for x in pattern_shuffled]);
			i2.attrs['pattern_shuffled'] = ",".join(["%1.2f" % x for x in pattern_shuffled]);
		
		sys.stdout.write(str(i1))
		sys.stdout.write(str(i2))
		
		
else:
	for i1, i2 in generator_doublebed(args.path):
		energy, pattern = get_rnahybrid(i2.attrs['seq'], i1.attrs['seq']);
		i1.attrs['energy'] = str(energy);
		i2.attrs['energy'] = str(energy);
		if(args.pattern):
			i1.attrs['pattern'] = ",".join([str(x) for x in pattern]);
			i2.attrs['pattern'] = ",".join([str(x) for x in pattern]);
		
		if(args.shuffle_trials):
			energy_shuffled, pattern_shuffled = get_shuffled_rnahybrid(i2.attrs['seq'], i1.attrs['seq'], trials=args.shuffle_trials);
			i1.attrs['energy_shuffled'] = str(energy_shuffled);
			i2.attrs['energy_shuffled'] = str(energy_shuffled);
			if(args.pattern):
				i1.attrs['pattern_shuffled'] = ",".join(["%1.2f" % x for x in pattern_shuffled]);
				i2.attrs['pattern_shuffled'] = ",".join(["%1.2f" % x for x in pattern_shuffled]);
			
		sys.stdout.write(str(i1))
		sys.stdout.write(str(i2))