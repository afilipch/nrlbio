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