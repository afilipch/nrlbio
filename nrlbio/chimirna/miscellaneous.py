#! /usr/bin/python	
'''library of miscellaneous functions'''
import random;


def intersection(a1, a2, b1, b2):
	if( int(a1) <= int(b1) and int(a2) >= int(b1)):
		return True;
	elif( int(a1) <= int(b2) and int(a2) >= int(b2)):
		return True;
	elif( int(a1) >= int(b1) and int(a2) <= int(b2)):
		return True;
	return False;
  
  
def shuffle_str(string):
	my_list = list(string);
	random.shuffle(my_list)
	return ''.join(my_list);
  
def findsubstring(string, substring, intersect = False):
	positions = [];
	pointer = 0;
	while True:
		pos = string[pointer:].find(substring)
		if(pos == -1):
			return positions;
		else:
			positions.append(pos + pointer);
			if(intersect):
				pointer += pos + 1;
			else:
				pointer += pos + len(substring);
	return positions; 
      
      
def max_values(my_list, my_key = lambda x: x):
	maximum = my_key(max(my_list, key= my_key));
	return filter(lambda x: my_key(x) == maximum, my_list)
  
def min_values(my_list, my_key = lambda x: x):
	minimum = my_key(min(my_list, key= my_key));
	return filter(lambda x: my_key(x) == minimum, my_list) 