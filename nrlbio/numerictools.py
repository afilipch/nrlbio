# /usr/bin/python
'''collections of classes and functions to solve diverse numerical problems'''
from math import log

import numpy as np;


def cdf(iterable):
	from collections import Counter
	'''Constructs cdf from given iterable with numeric entries
		iterable iterable: entries have to be of numeric type
	Returns list: each element is a tuple of two elements. 1st is any value from iterable, 2nd number of entries in iterable, which are less or equal to the 1st element.
	'''
	l = [];
	sum = 0;
	for k, v in sorted(Counter(iterable).items(), key = lambda x: x[0]):
		sum+=v;
		l.append((k, sum));
	return l;



def maxes(list_, key_function=lambda x: x):
	'''Returns all maximum(according to key function) elements in a given iterable.
	
		list_ sequence type: any sequence type to look for maximum elements in
		key_function function returns float: function to be called on the element of an iterable. The bigger the value function returns, the 'bigger' the elements
		
		Return tuple: 1st element is a list of maximum elements. 2nd element is an integer output of key_function value associated with maximum element.
	'''
	if(not list_):
		return [], None;
		
	m, max_list = key_function(list_[0]), []
	for s in list_:
		k = key_function(s)
		if k > m:
			m, max_list = k, [s]
		elif k == m:
			max_list.append(s)
	return max_list, m


def overlap(i1, i2):
	'''Returns overlap of two intervals, maybe negative
	
		i1 iterable: 1d interval. 2-element iterable. First element is start of interval(0-based inclusive). Second element is end of interval(0-based exclusive)
		i2 iterable: another 1d interval. 2-element iterable. First element is start of interval(0-based inclusive). Second element is end of interval(0-based exclusive)
		
		Returns iterable: Overlap of two intervals. 2-element iterable. First element is start of interval(0-based inclusive). Second element is end of interval(0-based exclusive).
	'''
	start = max(i1[0], i2[0])
	end = min(i1[1], i2[1])
	return start, end;
	
	
def merge_intervals(intervals, distance=0, assume_sorted=False):
	'''Yields lists of overlapping intervals
	
		intervals iterable: element is 2-element tuple(First element is start of interval(0-based inclusive). Second element is end of interval(0-based exclusive))
		distance int: minimum overlap(max gap if negative) required
		assume_sorted bool: if True, "intervals" argument is treated as sorted(according to start position) iterable
		
		Yields list: list of overlapping intervals (2-element tuple. First element is start of interval(0-based inclusive). Second element is end of interval(0-based exclusive))
	'''
	if(hasattr(intervals, 'next')):
		first = intervals.next();
		l = intervals
		
	else:	
		l = list(intervals)
		if(not assume_sorted):
			l.sort(key = lambda x: x[0]);	
		first = l.pop(0);
	
	merged = [first];
	start, end = first;
	for i in l:
		s, e = overlap(i, (start, end))
		if(e - s >= distance):
			merged.append(i);
			end = max(i[1], end)
		else:
			yield merged;
			start, end = i;
			merged = [i];
	yield merged		
		
		
def interval_intersection(i1, i2):
	'''Returns overlap of two intervals, only positive, or None
		
		i1 list|tuple: 1st element is the start of the interval(0-based), 2nd element is the end of the interval(exclusive)
		i2 list|tuple: 1st element is the start of the interval(0-based), 2nd element is the end of the interval(exclusive)
		
		Returns tuple|None: overlap of two intervals. 1st element start of the overlap(0-based), 2nd element end of the overlap(exclusive). None if there is no overlap
	'''	
	start, end = overlap(i1, i2);
	if(end>start):
		return start, end;
	else:
		return None;		
		
		
def overlap_hyperrectangles(hr1, hr2):
	'''Returns intersection of hyperrectangles
		
		hr1 list: hyperrectangle. List of intervals(2-element iterables: first element is start of interval(0-based inclusive), second element is end of interval(0-based exclusive))
		hr2 list: another hyperrectangle. List of intervals(2-element iterables: first element is start of interval(0-based inclusive), second element is end of interval(0-based exclusive))
		
		Returns list: overlap of two hyperrectangles. List of intervals(2-element iterables: first element is start of interval(0-based inclusive), second element is end of interval(0-based exclusive))
	'''
	c = [];
	for i1, i2 in zip(hr1, hr2):
		o = interval_intersection(i1, i2);
		if(o):
			c.append(list(o));
		else:
			return None
	return c;
	
	

		
		
def dict2entropy(d):
	'''Calculates entropy for a given dictionary
		
		d dict: Key may be everything, Value is a number of states(cases) assotiated with the Key;
		
		Returns float: entropy value;
	'''
	a = np.array(d.values(), dtype=float)
	probs = a/np.sum(a)
	
	return -1*sum([(p*log(p)) for p in probs])	
		
		
#testing section
if(__name__ == "__main__"):
	#print overlap([5,8], [4,10])
	#print overlap_hyperrectangles([[5,8], [4,10]],   [[3,513], [6,9]])
	
	a = [(1,9), (10,12), (11,16), (17,48),  (16,40), (100,125), (30, 35), (90,140), (9, 16), (22, 29), (130, 150)]
	for m in merge_intervals(a, distance = 1, assume_sorted=False):
		for i in m:
			print i;
		print	
	
	