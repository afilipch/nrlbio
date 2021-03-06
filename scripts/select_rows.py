#! /usr/bin/python
import sys;
import copy;
import argparse;


parser = argparse.ArgumentParser(description='output entry(line) of the file A only if the locker attribute value is present in the key field of of file B');

parser.add_argument('-a', '--locker_file', nargs = '?', required = True, type = str, help = "path to file A");
parser.add_argument('-b', '--key_file', nargs = '?', required = True, type = str, help = "path to file B");

parser.add_argument('-k', '--key', nargs = '?', required = True, type = int, help = "collumn in the second file used as key");
parser.add_argument('-l', '--locker', nargs = '?', required = True, type = int, help = "collumn in the first file used as locker");

parser.add_argument('-dl', '--delimeter_locker', nargs = '?', default = "\t", type = str, help = "delimeter of collums in locker file");
parser.add_argument('-dk', '--delimeter_key', nargs = '?', default = "\t", type = str, help = "delimeter of collums in key file");

parser.add_argument('-o', '--output', nargs = '?', default = "wa", choices = ["wa", "wab", "wi"],  type = str, help = "wa: outputs only A entries, keys for which were found in B\nwab: outputs A and corresponding B entries\nwi: outputs only A entries, keys for which were not found in B\n");
args = parser.parse_args();

key = args.key-1;
locker = args.locker-1;
kdelim = args.delimeter_key;
ldelim = args.delimeter_locker;


key2entry = {};
with open(args.key_file) as f:
	for l in f:
		a = l.strip().split(kdelim);
		if(args.output == "wab"):
			key2entry[a[key]] = copy.copy(a);
		else:
			key2entry[a[key]] = True;		

			
with open(args.locker_file) as f:
	for l in f:
		a = l.strip().split(ldelim);
		row = key2entry.get(a[locker], None);
		if(row): 
			if(args.output == "wa"):
				print "\t".join(a)
			if(args.output == "wab"):
				print "\t".join(a+row)
		elif(args.output == "wi"):
			print "\t".join(a)



	


