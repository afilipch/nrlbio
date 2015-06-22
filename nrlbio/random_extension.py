# /usr/bin/python
'''extends 'random' module functionality'''

import sys
import random
import bisect

import numpy as np


def list2interval(l):
	'''converts list of probabilities into interval, useable for random choice'''
	interval = np.array(l, dtype=float)
	interval = interval/np.sum(interval)	
	interval = np.cumsum(interval)	
	return interval
	

def weighted2interval(counter):
	'''creates list of items and probability intervals from counter. 
	For example, there are two items: A with probability 40 to be selected, B with probability 60. Then one should provide an counter [(A, 40), (B, 60)]. Outcome will be [A, B], [0.4, 1] 
	
		counter iterable: element is tuple. 1st element of the tuple is item to select, 2nd is a probability to select this item. Probabilities may be not normalized
	Returns:
		items list: list of items to select
		interval list: probability intervals. 
	'''
	items = [x[0] for x in counter]
	interval = list2interval([x[1] for x in counter])
	return items, interval;

	
def weighted_choice_fast(items, interval):
	'''randomly select item from items based on probability interval provided'''
	r = random.random()
	index = bisect.bisect_left(interval, r)
	return items[index];
	
	
def weighted_choice(counter):
	'''randomly select items from counter based on their probabilities
		
		counter iterable: element is tuple. 1st element of the tuple is item to select, 2nd is a probability to select this item. Probabilities may be not normalized
	
	Returns obj: randomly selected item
	'''
	items, interval = weighted2interval(counter);
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
	