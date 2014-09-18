# /usr/bin/python
'''extends 'random' module functionality'''

import sys
import random
import bisect

import numpy as np

def weighted2interval(iterable):
	items = [x[0] for x in iterable]
	interval = np.array([x[1] for x in iterable], dtype=float)
	interval = interval/np.sum(interval)	
	interval = np.cumsum(interval)
	return items, interval;

	
def weighted_choice_fast(items, interval):
	r = random.random()
	index = bisect.bisect_left(interval, r)
	return items[index];
	
	
def weighted_choice(iterable):
	items, interval = weighted2interval(iterable);
	return weighted_choice_fast(items, interval)
	
	
	
	
#testing section
if(__name__ == "__main__"):
	from collections import defaultdict
	a = [('a', 10), ('b', 50), ('c', 30), ('d', 80), ('e', 30)]
	choices = defaultdict(int);
	for _ in range(200):
		choices[weighted_choice(a)] += 1
	print choices;	
	
	#print bisect.bisect_left(weighted2interval(a)[1], 0.9)
	