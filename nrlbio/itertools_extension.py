# /usr/bin/python
'''library supporting advanced handling, processing, conversion of iterables'''

import sys
import os
import itertools
from collections import defaultdict

class LengthException(Exception):
	pass;

	
	
	
def pop_first(iterable):
	'''pops first element of any iterable, regardless does it has __getitem__ method or not'''
	if(hasattr(iterable, 'next')):
		return iterable.next();
	else:
		return iterable.pop(0);

		
def iterable_of_objects_to_counts_dict(iterable, attributes):
	'''Returns dictionary build on iterable. Key: tuple of certain values of attributes, Value: number of objects in iterable with values of attributes equal to the ones in Key.
	
		iterable iterable: any iterable of objects		
		attributes list: attributes used to pull oblects. For example if only attribute "age" is provided all people with the same age(35) will be gathered in one dictionary entry with Key=[(35)]
		if None all attributes will be used
		
		Returns dictionary: Key: tuple of certain values of attributes, Value: number of objects in iterable with values of attributes equal to the ones in Key.
	'''
	d = defaultdict(int);	
	first = pop_first(iterable);
	
	if(not attributes):
		attributes = vars(first).keys()
		
	k = tuple([getattr(first, attr) for attr in attributes])
	d[k]+=1;	
	
	for el in iterable:
		k = tuple([getattr(el, attr) for attr in attribute])
		d[k]+=1;
	return d;
	
	
def iterable_of_lists_to_counts_dict(iterable, indices=None):
	'''Returns dictionary build on iterable. Key: tuple of certain values of attributes(corresponding to indices), Value: number of objects in iterable with values of attributes equal to the ones in Key.
	
		iterable iterable: any iterable of objects		
		indices list: indices used to pull oblects. For example if only index 3 is provided all entries with the entry[3] will be gathered in one dictionary entry with Key=[(entry[3])]
		if None all attributes will be used
		
		Returns dictionary: Key: tuple of certain values of attributes(corresponding to indices), Value: number of objects in iterable with values of attributes equal to the ones in Key.
	'''
	d = defaultdict(int);
	
	first = pop_first(iterable);
	length = len(first);
	
	if(not indices):
		indices = range(length)
		
	k = tuple([first[i] for i in indices])
	d[k]+=1;	
		
	for el in iterable:
		if(len(el)!=length):
			raise LengthException('elements in iterable are of different length\n')
		else:		
			k = tuple([el[i] for i in indices])
			d[k]+=1;
	return d;	
	
	
def iterable_of_objects_to_list_dict(iterable, attributes=None):
	'''Returns dictionary build on iterable. Key: tuple of certain values of attributes, Value: list of objects in iterable with values of attributes equal to the ones in Key.
	
		iterable iterable: any iterable of objects		
		attributes list: attributes used to pull oblects. For example if only attribute "age" is provided all people with the same age(35) will be gathered in one dictionary entry with Key=[(35)]
		
		Returns dictionary: Key: tuple of certain values of attributes, Value: list of objects in iterable with values of attributes equal to the ones in Key.
	'''
	d = defaultdict(list);
	for el in iterable:
		k = tuple([getattr(el, key) for key in keys])
		d[k].append(el);
	return d;
	
	
def iterable_of_lists_to_list_dict(iterable, indices):
	'''Returns dictionary build on iterable. Key: tuple of certain values of attributes(corresponding to indices), Value: list of objects in iterable with values of attributes equal to the ones in Key.
	
		iterable iterable: any iterable of objects		
		indices list: indices used to pull oblects. For example if only index 3 is provided all entries with the entry[3] will be gathered in one dictionary entry with Key=[(entry[3])]
		if None all attributes will be used
		
		Returns dictionary: Key: tuple of certain values of attributes(corresponding to indices), Value: list of objects in iterable with values of attributes equal to the ones in Key.
	'''	
	d = defaultdict(list);
	
	first = pop_first(iterable);
	length = len(first);
	
	if(not indices):
		indices = range(length)
		
	k = tuple([first[i] for i in indices])
	d[k].append(first);	
		
	for el in iterable:
		if(len(el)!=length):
			raise LengthException('elements in iterable are of different length\n')
		else:		
			k = tuple([el[i] for i in indices])
			d[k].append(el);
	return d;