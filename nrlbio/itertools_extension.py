# /usr/bin/python
'''library supporting advanced handling, processing, conversion of iterables'''

import sys
import os
import itertools
from collections import defaultdict

class LengthException(Exception):
	pass;
	
	
class TypeException(Exception):
	pass;	

	
	
	
def pop_first(iterable):
	'''pops first element of any iterable, regardless does it has __getitem__ method or not'''
	if(hasattr(iterable, 'next')):	
		return iterable.next();
	else:
		return iterable.pop(0);

		
def iterable_of_objects_to_counts_dict(iterable, attributes=None):
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
		k = tuple([getattr(el, attr) for attr in attributes])
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
	
	
def iterable_to_counts_dict(iterable, attributes):
	'''Returns dictionary build on iterable. Key: tuple of certain values of attributes, Value: number of objects in iterable with values of attributes equal to the ones in Key.
	
		iterable iterable: any iterable of objects		
		attributes list: attributes used to pull oblects. For example if only attribute "age" is provided all people with the same age(35) will be gathered in one dictionary entry with Key=[(35)]
		OR indices used to pull oblects. For example if only index 3 is provided all entries with the entry[3] will be gathered in one dictionary entry with Key=[(entry[3])]
		if None all attributes will be used
		
		Returns dictionary: Key: tuple of certain values of attributes, Value: number of objects in iterable with values of attributes equal to the ones in Key.
	'''
	if(not attributes):
		try:
			return iterable_of_lists_to_counts_dict(iterable)
		except:
			return iterable_of_objects_to_counts_dict(iterable)
			
	elif(all([isinstance(x, str) for x in attributes])):
		return iterable_of_objects_to_counts_dict(iterable, attributes=attributes);
		
	elif(all([isinstance(x, int) for x in attributes])):
		return iterable_of_lists_to_counts_dict(iterable, indices=attributes)
		
	else:
		raise TypeException('atributes have to all integers or all strings\n')
	
	
def iterable_of_objects_to_list_dict(iterable, attributes=None):
	'''Returns dictionary build on iterable. Key: tuple of certain values of attributes, Value: list of objects in iterable with values of attributes equal to the ones in Key.
	
		iterable iterable: any iterable of objects		
		attributes list: attributes used to pull oblects. For example if only attribute "age" is provided all people with the same age(35) will be gathered in one dictionary entry with Key=[(35)]
		if None all attributes will be used
		
		Returns dictionary: Key: tuple of certain values of attributes, Value: list of objects in iterable with values of attributes equal to the ones in Key.
	'''
	d = defaultdict(list);
	first = pop_first(iterable);
	
	if(not attributes):
		attributes = vars(first).keys()
		
	k = tuple([getattr(first, attr) for attr in attributes])
	d[k].append(first);	
	
	for el in iterable:
		k = tuple([getattr(el, attr) for attr in attributes])
		d[k].append(el);
	return d;
	
	
def iterable_of_lists_to_list_dict(iterable, indices=None):
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
	
	
	
	
#testing section
if(__name__ == "__main__"):
	from random import randint
	
	class A(object):
		def __init__(self):
			self.a = randint(1,10);
			self.b = randint(1, 5);
			self.c = bool(randint(0,1));
		
	it_objects = (A() for _ in range(1000))
	it_lists = [(randint(1,4), randint(1,4), randint(1,4)) for _ in range(1000)]
	
	for kv in iterable_of_objects_to_list_dict(it_objects, ['a', 'c']).iteritems():
		print "%s|\t\t%s\n" % (kv[0], "\t".join([str(vars(x)) for x in kv[1]]));
	
	
	for kv in iterable_of_lists_to_list_dict(it_lists, [0, 2]).iteritems():
		print "%s|\t\t%s\n" % (kv[0], "\t".join([str((x)) for x in kv[1]]));
		
	#it_lists.append((1,7));
	
	#for kv in iterable_of_lists_to_list_dict(it_lists, [0, 2]).iteritems():
		#print "%s|\t\t%s\n" % (kv[0], "\t".join([str((x)) for x in kv[1]]));	
	
	
	
	
	
	