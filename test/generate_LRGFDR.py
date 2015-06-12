# /usr/bin/python
'''creates grid to test LRGFDR'''

import argparse;
from collections import defaultdict
from random import randint, expovariate, gauss

from nrlbio import LRGFDR;
from visualize_LRGFDR import  visualize_grid;

parser = argparse.ArgumentParser(description='creates grid to test LRGFDR');
parser.add_argument('path', nargs = '?',  type = str, help = "name of the output file");
args = parser.parse_args();

c_elements = 10000;

#2D test
def add_elements(d, coordinates, elements):
	for _ in range(elements):
		key = [];
		for s, e in coordinates:
			key.append(randint(s, e-1));
		d[tuple(key)] += 1;	
	return True;


def exp(alpha, lower, upper):
	for _ in range(20):
		r = lower + int(expovariate(alpha)*(upper-lower))
		if(r<upper):
			return int(r);
	return lower;	

def gaussian(mu, var, lower, upper):
	cvar = var*((upper-lower)**0.5)
	for _ in range(20):
		r = gauss(mu, cvar);
		if(lower<=r<upper):
			return int(r);
	return lower;

#non-unfiorm distribution
def add_distribution(d, distributions, elements, arguments):
	for _ in range(elements):
		key = [];
		for f, args in zip(distributions, arguments):
			key.append(f(*args));
		d[tuple(key)] += 1;
	return True;

	
#we use 100*100 table;
control = defaultdict(int);
signal = defaultdict(int);
	

def overlap_problem(signal,control, c_elements):
	add_elements(control, [[20, 120], [50, 150]], c_elements);
	add_elements(signal, [[20, 120], [50, 150]], c_elements);
	add_elements(signal, [[20, 60], [50, 80]], 100000);
	add_elements(signal, [[40, 80], [70, 110]], 40000);

def negative_problem(signal, control, c_elements):
	add_elements(control, [[20, 120], [50, 150]], c_elements);
	add_elements(signal, [[20, 120], [50, 150]], c_elements);	
	add_elements(signal, [[50, 60], [60, 100]], 60000);
	add_elements(signal, [[30, 60], [90, 100]], 60000);
	
	
def distribution(signal, control, c_elements):
	add_distribution(signal, [gaussian, exp], 60000, [(30, 1.5, 0, 60), (3,0,40)])
	add_distribution(signal, [gaussian, randint], c_elements, [(50,1,0,100), (0,100)])
	add_distribution(control, [gaussian, randint], c_elements, [(50,1,0,100), (0,100)])
	#add_elements(control, [[0, 100], [0, 100]], c_elements);
	#add_elements(signal, [[0, 100], [0, 100]], c_elements);
	#d = defaultdict(int);
	#for _ in range(10000):
		#k = gaussian(30, 1, 10, 50);
		#d[k] += 1;
	#for k, v in sorted(d.items(), key = lambda x: x[1]):
		#print k, v
		

#negative_problem(signal, control, c_elements);	
distribution(signal, control, c_elements)	
grid = LRGFDR.Grid.from_dict(signal, control);
grid.to_csv(args.path);



