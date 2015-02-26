#! /usr/lib/python
'''Creates makefile and directory structure for chiflex project'''
import argparse
import sys
import os

bowtie_settings = {'N': ('0','-'),
'L': ('16','-'),
'i': ('C,2', '-'),
'ignore-quals': ('True', '--'),
'norc': ('True', '--'),
'local': ('True', '--'), 
'mp': ('6,6', '-'),
'rfg': ('12,12', '--'),
'rdg': ('6,4', '--'),
'min-score': ('C,40', '--'), 
'k': ('4', '-'),
'D': ('30', '-'),
'R': ('2', '-'),
'no-unal': ('True', '--'), 
'p': ('8', '-') }
bs_string = "\n".join(["\t%s%s=%s" % (x[1][1], x[0], x[1][0]) for x in bowtie_settings.items()])

parser = argparse.ArgumentParser(description='assignes hybridization energy to interactions');
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the project folder");
parser.add_argument('-b', '--bowtie', nargs = '+', default = [], type = str, help = "Bowtie settings. For example, if one wants to set \'-p 4\', use \'--local\' alignment mode, but not \'--norc\' option then \'p=4 local=True norc=False\' should be provided. Given attributes replace default(for Chiflex, NOT for Bowtie) ones.\nDefault settings are:%s" % bs_string)
args = parser.parse_args();

#process bowtie options
for bo in args.bowtie:
	try:
		name, value = bo.split("=");
	except(ValueError):
		sys.exit("Bowtie options are provided in a malformatted way. Please see help for example\n") 
		
	if(name in bowtie_settings):
		bowtie_settings[name] = (value, bowtie_settings[name][1])
	else:
		sys.stderr.write("provided option \'%s\' is currently not supported by Chiflex and will be ignored\n" % name) 

#print bowtie_settings