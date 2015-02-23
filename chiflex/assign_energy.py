#! /usr/lib/python
'''assignes hybridization energy to interactions'''

import argparse
import sys

from nrlbio.generators import generator_doublebed
from nrlbio.rnahybrid import get_rnahybrid


parser = argparse.ArgumentParser(description='assignes hybridization energy to interactions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to doublebed/gff file");
args = parser.parse_args();


for i1, i2 in generator_doublebed(args.path):
	energy, dummy = get_rnahybrid(i1.attrs['seq'], i2.attrs['seq']);
	i1.attrs['energy'] = str(energy);
	i2.attrs['energy'] = str(energy);
	sys.stdout.write(str(i1))
	sys.stdout.write(str(i2))