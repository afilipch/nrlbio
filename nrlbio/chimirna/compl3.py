'''assess 3' complementarity for each miRNA'''
import interaction_lib;
import sys;
import os;
import argparse;
import copy
from collections import *;
from Bio.Seq import reverse_complement;


parser = argparse.ArgumentParser(description='assess 3\' complementarity for each miRNA');
parser.add_argument('path', metavar = 'N', nargs = '+', type = str, help = "path to interactions (output/interactions.bed)");
parser.add_argument('-w', '--window', nargs = '?', default = 4, type = int, help = "length of complementarity window")
parser.add_argument('-s', '--start', nargs = '?', default = 10, type = int, help = "position in miRNA to start from looking for complementarity")
args = parser.parse_args();

def sum_lists(*args):
	ml = max([len(x) for x in args]);
	r = [0]*ml;
	for l in args:
		for i, v in enumerate(l):
			r[i]+=v;
	return r;
	
def negate_list(l1, l2):
	r = [];
	for i, v in enumerate(l1):
		if(i<len(l2)):
			r.append(v-l2[i]);
		else:
			r.append(v);
	return r;
	
	
def max_lists(*args):
	ml = max([len(x) for x in args]);
	mymin = min([min(x) for x in args])
	r = [mymin]*ml;
	for l in args:
		for i, v in enumerate(l):
			r[i]=max([v, r[i]]);
	return r;	
	
def mean_lists(*args):
	ml = max([len(x) for x in args]);
	mymin = min([min(x) for x in args])
	r = [0]*ml;
	a = [0]*ml
	for l in args:
		for i, v in enumerate(l):
			r[i]+= v;
			a[i]+=1.0;
	for i in range(len(r)):
		r[i] = r[i]/a[i]
	return r;	
	
	
def generate_mismatch(match):	
	r = set()
	for i in range(6):
		for n in "ACTG":
			r.add(match[:i] + n + match[i:])
	return r;	
	
def decrease_tseq(tseq, mseq):
	seeds = generate_mismatch(reverse_complement(mseq[1:7]));
	for seed in seeds:
		pos = tseq.rfind(seed);
		if(pos>-1):
			return tseq[:pos]
	return tseq;		
		
	
	
	
def get_profile(mseq, tseq, start, lw):
	t = reverse_complement(decrease_tseq(tseq,mseq));
	m = mseq[start:]
	r = [0]*len(m);
	for i in range(len(m)-lw + 1):
		w = m[i:i+lw];
		#print w, t
		if(w in t):
			for j in range(i, i+lw):
				r[j]+=1;
	return r;			
	
def get_profiles(interactions, start, lw, support = 20):
	mir2profile = {};
	mir2profiles = defaultdict(list);
	for inter in interactions:
		mir2profiles[inter.mirid].append(get_profile(inter.mirseq, inter.tseq, start, lw))
		#print get_profile(inter.mirseq, inter.tseq, start, lw)
	for k, v in mir2profiles.iteritems():
		if(len(v)>support):
			mir2profile[k] = sum_lists(*v);
	return mir2profile;	
		
		
def get_profiles_control(interactions, start, lw, repeat = 10):
	mir2profile = {};
	for stub in range(repeat):
		cinteractions = interaction_lib.shufseq(interactions);
		temp = defaultdict(list)
		mir2profiles = defaultdict(list);
		for inter in cinteractions:
			mir2profiles[inter.mirid].append(get_profile(inter.mirseq, inter.tseq, start, lw))
			#print get_profile(inter.mirseq, inter.tseq, start, lw)
		for k, v in mir2profiles.iteritems():
			temp[k].append(sum_lists(*v));
		for k, v in temp.iteritems():
			mir2profile[k] = max_lists(*v);
	return mir2profile		
	
def clear_signal(ls, lc):
	for i, sv in enumerate(ls):
		if (sv-lc[i] > 10):
			return True;
	else:
		return False;

mir2compl = defaultdict(list);
interactions = interaction_lib.get(args.path[0]);


signal = get_profiles(interactions, args.start, args.window)
control = get_profiles_control(interactions, args.start, args.window, repeat = 10)

mirid2seq = dict([(x.mirid, x.mirseq) for x in interactions])

for mirid, compl in signal.iteritems():
	if(clear_signal(compl, control[mirid])):
		print "%s\t%s\n%s\n%s\n%s\n" % (mirid, mirid2seq[mirid], compl, control[mirid], negate_list(compl, control[mirid]));
		
dirname = 	"mirna_profiles"	
if not os.path.exists(dirname):
    os.makedirs(dirname)		

#for path in args.interactions:
	#f = open(path);
	#for l in f:
		#a = l.strip().split("\t");
		#mirids = a[6].split(",");
		#if(len(mirids) == 1):
			#mirid = mirids[0];
			#r = get_profile(a[7], a[8], start = 10, lw =4)
			#mir2compl.append(copy.copy(r))
			
			
		#int2mirseq[a[3]] = a[7]
	#f.close()	