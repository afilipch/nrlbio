# /usr/bin/python
'''collection of basic statistics functions'''

import scipy.stats as st
from math import sqrt


def wilson_score_interval_binomial(k, n, pval):
	'''Calculates confidence interval for binomial proportion(probability of success) based on one sample observation
	
		k int: number of successes in a sample
		n int: size of a sample
		pval float: significance level (proabilty that actual binomial proportion lies outside reported interval)
		
	Returns tuple: lower and upper bounfary of the confidence interval for binomial proportion(probability of success) based on one sample observation
	'''
	
	if(n==0):
		return 0
	
	z = st.norm.ppf(1-pval/2)
	p = float(k)/n
	
	d = z*sqrt((p*(1-p)+z*z/(4*n))/n)
	
	a = p + z**2/(2*n)
	c = 1+z**2/n

	return (a-d)/c, (a+d)/c
	
	
	
	
	
if(__name__=='__main__'):
	print wilson_score_interval_binomial(1, 1000, 0.01);