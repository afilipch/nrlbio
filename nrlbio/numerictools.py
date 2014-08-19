# /usr/bin/python
'''collections of classes and functions to solve diverse numerical problems'''


def maxes(list_, key_function=lambda x: x):
	'''Returns all maximum(according to key function) elements in a given iterable.
	
		list_ sequence type: any sequence type to look for maximum elements in
		key_function function returns float: function to be called on the elemenent of an iterable. The bigger the value function returns, the 'bigger' the elements
		
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
	'''Returns overlap of two intervals
	
		i1 iterable: 1d interval. 2-element iterable. First element is start of interval(0-based inclusive). Second element is end of interval(0-based exclusive)
		i2 iterable: another 1d interval. 2-element iterable. First element is start of interval(0-based inclusive). Second element is end of interval(0-based exclusive)
		
		Returns iterable: Overlap of two intervals. 2-element iterable. First element is start of interval(0-based inclusive). Second element is end of interval(0-based exclusive). Returns None if there is no interval
	'''
	l = max(i1[0], i2[0])
	u = min(i1[1], i2[1])
	if(u>l):
		return l, u
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
		o = overlap(i1, i2);
		if(o):
			c.append(list(o));
		else:
			return None
	return c;	
		
#testing section
if(__name__ == "__main__"):
	print overlap([5,8], [4,10])
	print overlap_hyperrectangles([[5,8], [4,10]],   [[3,513], [6,9]])