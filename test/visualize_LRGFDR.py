import os

import numpy as np
import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import animation
from pylab import Rectangle

from nrlbio import LRGFDR;
from nrlbio.numpy_extension import get_slice

def array2axes(grid, array):
	ax = plt.axes()
	ax.pcolormesh(array, cmap=plt.cm.Blues, alpha=0.8)
	
	#adjust layout
	ax.set_frame_on(False)
	
	#set ticks and ticklabels

	xticks = ax.get_xticks();
	#xticks[-1] -= 1;
	ax.set_xticks(xticks+0.5, minor=False)
	ax.set_xticklabels([grid.decoding_table[1][int(x)] for x in xticks[:-1]], minor=False);
	
	yticks = ax.get_yticks();
	#yticks[-1] -= 1;
	#print yticks
	ax.set_yticks(yticks+0.5, minor=False)
	ax.set_yticklabels([grid.decoding_table[0][int(x)] for x in yticks[:-1]], minor=False);	

	#set lables:
	if(not grid.attribute_names):
		grid.set_attribute_names(["x1","x2"])
	ax.set_xlabel(grid.attribute_names[0]);
	ax.set_ylabel(grid.attribute_names[1]);	
	
	# Turn off all the ticks
	ax = plt.gca()
	for t in ax.xaxis.get_major_ticks():
		t.tick1On = False
		t.tick2On = False
	for t in ax.yaxis.get_major_ticks():
		t.tick1On = False
		t.tick2On = False
	return ax;	

def grid2heatmap(grid, f):
	a = f(grid.array)
	a = (a - a.mean()) / (a.max() - a.min())
	ax = array2axes(grid, a)	
	return ax

	
	
	
	
def animate_clustering(grid, function, clusters, nclusters, output):
	a = function(grid.array)
	a = (a - a.mean()) / (a.max() - a.min())
	_max = a.max()
	_min = a.min()
	coordinates = [None];
	for cl in clusters:
		coordinates += cl.history;
	lc = len(coordinates)
	for cl in nclusters:
		coordinates += cl.history;
	frames = len(coordinates);
	#return None
	
	#def init():
		#print "fuck you!"
		#return grid2heatmap(grid, function)

	# animation function.  This is called sequentially
	def animate(i):
		if(i==0):
			return grid2heatmap(grid, function)
				
		aslice = get_slice(a, coordinates[i])
		if(i<lc):
			aslice[:] = _max
		else:
			aslice[:] = _min
		print i;
		ax = plt.axes()
		ax.pcolormesh(a, cmap=plt.cm.Blues, alpha=0.8)
		return ax		
		#return grid2heatmap(grid, function)
		
	fig, ax = plt.subplots()
	fig = plt.gcf()
	fig.set_size_inches(10, 10)	
	# call the animator.  blit=True means only re-draw the parts that have changed.
	anim = animation.FuncAnimation(fig, animate,
								frames=frames, interval=500, blit=True)
	anim.save(output, fps=1)
	
	
