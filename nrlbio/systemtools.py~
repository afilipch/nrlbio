'''Collection of functions to interact with an operation system'''


import os

def onlyfiles(directory):
	'''Returns list of all files in a given directory. NOTE: files in the subfolders are NOT returned'''
	return [ f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory,f)) ]
	
	
def allfiles(directory)	:
	'''Returns list of all files in a given directory. NOTE: files in the subfolders are returned'''
	f = []
	for (dirpath, dirnames, filenames) in os.walk(directory):
		f.extend(filenames)	
	return f;	