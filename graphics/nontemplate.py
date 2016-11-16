# /usr/bin/python
'''Draws a plot of nontemplate addiction to miRNAs'''
import sys
import argparse
from collections import defaultdict
import math

from pybedtools import BedTool
import matplotlib.pyplot as plt
import numpy as np

#from nrlbio.pyplot_extension import remove_top_left_boundaries


parser = argparse.ArgumentParser(description='Draws a plot of nontemplate addiction to miRNAs');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the interacting RNA, double gff file with mode and mode_shuffled assigned");
parser.add_argument('--norm', nargs = '?', default=False, const = True, type = bool, help = "If set, binding modes will be normalized");
parser.add_argument('--output', nargs = '?', type = str, help = "Path to the output");
parser.add_argument('--extlen', nargs = '?', required=True, type = int, help = "The length of downstream extensions");
#parser.add_argument('--title', nargs = '?', type = str, help = "Title for a plot");
#parser.add_argument('--control', nargs = '?', default = 'shuffled sequences', type = str, help = "Type of control: shuffled sequences or shuffled interactions");
args = parser.parse_args();





