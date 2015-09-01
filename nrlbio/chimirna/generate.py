#! usr/bin/python
'''this script generates clash pipeline makefile and creates required folder tree'''
from project_lib import Project;
import gconst;
import os;
import sys;
import argparse




parser = argparse.ArgumentParser(description='script generates clash pipeline makefile and creates required folder tree');
parser.add_argument('--name', nargs = '?', type = str, required = True, help = "name of the project, will be used as name for root folder of the project created in current directory");
parser.add_argument('--reads', nargs = '?', type = str, required = True, help = "path to collapsed reads fastq file");
parser.add_argument('--ref', nargs = '?', type = str, required = True, help = "path to mapping reference bed file");
parser.add_argument('--index', nargs = '?', type = str, required = True, help = "path to mapping reference index");
parser.add_argument('--system', nargs = '?', type = str, choices = gconst.systems, required = True, help = "model system");
parser.add_argument('--mode', nargs = '?', type = str, choices = gconst.experiment.keys(), required = True, help = "type of experiment");
parser.add_argument('--gmode', nargs = '?', default = False, const = True, type = bool, help = "gmode enabled, we account for RNAse T1 preferences");
parser.add_argument('--explore', nargs = '?', default = False, const = True, type = bool, help = "slow exploration mode");
parser.add_argument('--change', nargs = '?', default = False, const = True, type = bool, help = "only change or create makefile");
args = parser.parse_args();


name = os.path.abspath(args.name);
if(not args.change):
	while True:
		try:
			os.makedirs(name);
			os.mkdir(os.path.join(name, "log"));
			os.mkdir(os.path.join(name, "left"));
			os.mkdir(os.path.join(name, "lstatistics"));
			os.mkdir(os.path.join(name, "right"));
			os.mkdir(os.path.join(name, "rstatistics"));
			os.mkdir(os.path.join(name, "sam"));
			os.mkdir(os.path.join(name, "output"));
			os.mkdir(os.path.join(name, "plots"));
			if(args.gmode):
				os.mkdir(os.path.join(name, "gmode"));
			break;
		except:
			name = os.path.abspath(raw_input("directory with this name is currently exists, please give another name: "))

p = Project(os.path.basename(name), args.reads, args.ref, args.index, args.system, args.mode, args.gmode, args.explore);		
wf = open(os.path.join(gconst.proot, p.name + ".txt"), 'w');
wf.write("python " + " ".join(sys.argv) + "\n\n");
wf.write(p.info())
wf.close


mh = open(os.path.join(name, "Makefile"), 'w');
mh.write(p.makefile());
mh.close();