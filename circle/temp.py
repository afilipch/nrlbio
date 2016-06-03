from itertools import combinations
from collections import defaultdict
from nrlbio.itertools_extension import powerset


routes = defaultdict(list)


for p1, p2 in combinations(range(4), 2):
	if(p1==1 and p2==2):
		pass;
		##pairs[p1].append(p2);
	else:	
		routes[p2].append([p1])
		for el in routes[p1]:
			routes[p2].append(el+[p1])
	
	
for k, l  in routes.items():
	for v in l:
		print "\t".join([str(x) for x in v+[k]])
		 
		 
print 

for ss in powerset(range(4)):
	if(len(ss)>1):
		print "\t".join([str(x) for x in ss])