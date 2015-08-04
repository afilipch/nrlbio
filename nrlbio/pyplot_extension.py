# /usr/bin/python
'''Collection of classes and functions for plotting. Extends pyplot functionality'''
from collections import Counter, Iterable

import matplotlib
#matplotlib.use('pdf')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from nrlbio.itertools_extension import flatten

try:
	plt.style.use('ggplot')
except:
	pass;

def _set_bins(values, step = None, range_=None):
	if(step):
		if(range_):
			return int(round((range_[1]-range_[0])/step))
		else:
			return int(round((max(values)-min(values))/step+1))
	return len(values)
	
def _set_range(values, bins, step = None):
	if(step):
		return min(values), max(values)+step
	else:	
		ext = (max(values)-min(values))/(bins-1)
		return min(values), max(values)+ext
	
def _set_axis_extension(xvals, yvals, xext=0.1, yext=0.1):
	xext_ = (max(xvals) - min(xvals))*xext
	return min(xvals)-xext_, max(xvals)+xext_, 0, max(yvals)*(1+yext)
	
	
def _get_weights(data):
	if(isinstance(data, dict)):
		values = data.keys()
		weights = data.values()
	elif(isinstance(data, Iterable)):
		d = Counter(data)
		values = d.keys()
		weights = d.values()
	else:
		raise ValueError("{data} argument should be dictionary or iterable %s type is provided instead" % str(type(data)))
	return values, weights
	

def histogram(data, title=None, xlabel=None, ylabel=None, xticks=None, xticklabels=None, xticksrotation = 0, output=None, step = None, bins=None, range=None, normed=False, cumulative=False, bottom=None, histtype=u'bar', align=u'mid', orientation=u'vertical', rwidth=None, log=False, color=None, label=None, stacked=False, hold=None):
	'''Draws histogram with predifined style on basis of provided data
		
		data list|numpy.ndarray|dict: list or array of values to plot as histogram. Counter-like dict may be also provided.
		output str|None: path to the output pdf file. If None, plot will be shown in interactive mode 	
	'''
	values, weights = _get_weights(data);
	if(not bins):
		bins = _set_bins(values, step=step, range_=range);
	if(not range):
		range = _set_range(values, bins, step=step);
		
	fig = plt.figure()
	ax = plt.subplot(111)
	ax.patch.set_facecolor('white')
	
	yvals, xvals, npathes = plt.hist(values, bins=bins, range=range, normed=normed, weights=weights, cumulative=cumulative, bottom=bottom, histtype=histtype, align=align, orientation=orientation, rwidth=rwidth, log=log, color=color, label=label, stacked=stacked, hold=hold);
	
	
	plt.axis(_set_axis_extension(xvals, yvals))
	
	if(label):
		plt.legend(loc='upper right')
	if(xlabel):
		plt.xlabel(xlabel)
	if(ylabel):	
		plt.ylabel(ylabel)
	if(title):	
		plt.title(title)
	if(xticks):
		plt.xticks(xticks)
	if(xticklabels):
		ax.set_xticklabels(xticklabels, rotation=xticksrotation)	
	
	
	if(output):
		plt.savefig(output, bbox_inches='tight')
	else:
		plt.show(bbox_inches='tight')
		
	plt.close();
	return yvals, xvals, npathes


	
def multihistogram(data, title=None, xlabel=None, ylabel=None, xticks=None, xticklabels=None, xlabels=None, xticksrotation = 0, output=None, step = None, bins=None, range=None, normed=False, cumulative=False, bottom=None, histtype=u'bar', align=u'mid', orientation=u'vertical', rwidth=None, log=False, color=None, label=None, stacked=False, hold=None):
	'''Draws many histograms in one figure with predifined style on basis of provided data
		
		data list|numpy.ndarray|dict: list or array of values to plot as histogram. Counter-like dict may be also provided.
		output str|None: path to the output pdf file. If None, plot will be shown in interactive mode 
	'''
	vw = [];
	for datum in data:
		vw.append(_get_weights(datum));
		
	merged_values = set(flatten([x[0] for x in vw]))
		
	if(not bins):
		bins = _set_bins(merged_values, step=step);
	if(not range):
		range = _set_range(merged_values, bins, step=step);
		
	
	
	fig = plt.figure()
	ax = plt.subplot(111)
	#ax.patch.set_facecolor('white')
	yxp = [];
	for (values, weights), c, l in zip(vw, color, label):
		yxp.append(plt.hist(values, bins=bins, range=range, normed=normed, weights=weights, cumulative=cumulative, bottom=bottom, histtype=histtype, align=align, orientation=orientation, rwidth=rwidth, log=log, color=c, label=l, stacked=stacked, hold=hold, alpha=0.5));
	
	yvals = set(flatten([x[0] for x in yxp]))
	xvals = set(flatten([x[1] for x in yxp]))
	plt.axis(_set_axis_extension(xvals, yvals))
	plt.legend(loc='upper right')
	if(xlabel):
		plt.xlabel(xlabel)
	if(ylabel):	
		plt.ylabel(ylabel)
	if(title):	
		plt.title(title)
	if(xticks):
		plt.xticks(xticks)
	if(xticklabels):
		ax.set_xticklabels(xticklabels, rotation=xticksrotation)
		
	if(output):
		plt.savefig(output)
	else:
		plt.show()
		
	plt.close();	
	
	
def _set_pie_slices(data, labels, top=None, min_fraction=0):
	if(isinstance(data, dict)):
		ncounter = list(data.items());
	elif(isinstance(data, Iterable)):
		ncounter = list(zip(labels, data));
	else:
		raise ValueError("{data} argument should be dictionary or iterable %s type is provided instead" % str(type(data)));	
	ncounter.sort(key= lambda x: x[1], reverse=True)
	
	others = [];
	total = float(sum([x[1] for x in ncounter]))
	if(top):
		others = ncounter[top:]
		ncounter = ncounter[:top];
	if(min_fraction):
		tc = [];
		for l, v in ncounter:
			if(v/total>=min_fraction):
				tc.append((l,v));
			else:
				others.append((l, v));
		ncounter = tc;
		
	if(others):
		ncounter.append(('others', sum([x[1] for x in others])))

	
	slices = [x[1] for x in ncounter]
	labels = ['%s\nn=%d' % x for x in ncounter]
			
	return slices, labels, others
	

def _set_wedge_label(w, min_fraction=0.04):
	if(w/100.0<min_fraction):
		return ''
	else:	
		return "%1.1f%%" % w
	
_num_others = 5;	
	
def pie(data, top=10, min_fraction=0.04, title=None, labelsize='smaller', labelva='top', output=None, explode=None, labels=None, colors=('yellowgreen', 'gold', 'lightskyblue', 'lightcoral', 'azure', 'seashell', 'darkorchid', 'chartreuse'), autopct=_set_wedge_label, pctdistance=0.8, shadow=False, labeldistance=1.15, startangle=0, radius=None, counterclock=True, wedgeprops=None, textprops={'fontsize': 'smaller'}):
      
	slices, labels, others= _set_pie_slices(data, labels, top=top, min_fraction=min_fraction);
	
	#fix for others	
	if(others):
		startangle = 271 + (float(slices[-1])/sum(slices))*180
		others.sort(key= lambda x: x[1], reverse=True)
		labels[-1] = "\n".join(['%s(n=%d)' % x for x in others[:_num_others]])
	 
	plt.figure(1, figsize=(8,8))
	ax = plt.axes([0.2, 0.2, 0.6, 0.6])
	patches, texts, autotexts = plt.pie(slices, explode=explode, labels=labels, colors=colors, autopct=autopct, pctdistance=pctdistance, shadow=shadow, labeldistance=labeldistance, startangle=startangle, radius=radius, counterclock=counterclock, wedgeprops=wedgeprops, textprops=textprops)
	
	#plt.axis('equal')
	if(title):	
		plt.title(title)
	if(labelsize):
		for text in texts:
			text.set_fontsize(labelsize)
	if(labelva):
		for text in texts:
			text.set_verticalalignment(labelva)
	
	if(output):
		plt.savefig(output)
	else:
		plt.show()
		
	plt.close();	
	
	
	
	
	
	
if(__name__ == '__main__'):
	data = [1,1,1,1,4,5,6,5,5,5,6,6,6,7,7,2,2,2,3,3,3,3,1,1,1,1,4,4,4, 12, 12]
	c = Counter(data);
	di = {6: 50, 6.1: 150, 6.2: 250, 6.3: 150, 6.4: 50, 6.8: 70};
	dj = {5.7: 40, 5.8: 80, 5.7: 100, 6: 120, 6.1: 180, 6.2: 50};
	#multihistogram([di, dj], title = 'konfetka', xlabel='day time', ylabel='bu bu bu level', xticks=list(np.arange(5.5, 6.9, 0.1)), xticklabels= ['%dh' % x for x in range(5,20)], step=0.1, color=('lightgreen', '0.15'), labels=('signal', 'control'));
	labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
	sizes = [15, 30, 45, 10]
	pie(dj, labels=labels, title="be", top = 10, min_fraction= 0.1)
	
		
		
	
