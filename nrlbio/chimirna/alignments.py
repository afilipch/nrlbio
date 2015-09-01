#! /usr/bin/python	
def get_best_alignments(refseq_dict, sid, seq, m=2, mm=-5, g=-6, e=-4):
	'''based on Bio.pairwise.align.local outputs list of alignments to refseq_dict with best scores'''
	import pairwise2;
	mymax = 0;
	ans = [];
	for rid, refseq in refseq_dict.items():
		for a in pairwise2.align.localms(refseq, seq, m, mm, g, e):
			if(a[2] > mymax):
				b = [sid, rid] + list(a);
				b.append(len(refseq))
				b.append(len(seq));
				ans = [b]
				mymax = a[2]
			elif(a[2] == mymax):
				b = [sid, rid] + list(a);
				b.append(len(refseq))
				b.append(len(seq))
				ans.append(b);
	for el in ans:
		el.append(1.0/len(ans))
	return ans;
	
def analyse_alignment(alignment, represent = False, add_conversion = ''):
	r_a, s_a, score, start, end, r_length, s_length, weight = alignment[2:10];
	s_start, r_start, s_end, r_end = 0,0,0,0 ## start of match part in read first nt 0based and is included in match, start of match part in reference first nt 0based and is included in match, end of match part in read first nt 0based and is excluded in match, end of match part in reference first nt 0based and is excluded in match;
	conv = [] ## conversions in form (position relative to reference, nucl/gap in reference, nucl/gep in read)
	if(add_conversion):
		conv.append(add_conversion)
	cut_left, cut_right = "NN", "NN" ## cut is represented by two letters first is left part of cut, second is right;
	
	# search for start positions
	for i in range(start):
		if(r_a[i] != "-"):
			r_start += 1;
			r_end += 1;
		if(s_a[i] != "-"):
			s_start += 1;
			s_end += 1;
	#print r_start, s_start
	
	for i in range(start, end):
		if(r_a[i] != s_a[i]):
			conv.append((r_end, r_a[i], s_a[i]))		
		if(r_a[i] != "-"):
			r_end += 1;
		if(s_a[i] != "-"):
			s_end += 1;

		
	if(r_start and r_a[start-1] != "-"):
		cut_left = r_a[start-1] + r_a[start];
	if(len(r_a) > end and r_a[end] != "-"):
		cut_right = r_a[end-1] + r_a[end];	
		
	#print r_start, s_start			
	#print r_end, s_end
	#for el in conv:
		#print el;
	#print cut_left, cut_right			
	if(represent):
		if(conv):
			conv_str = "|".join([",".join([str(y) for y in x]) for x in conv])
		else:
			conv_str = "-1,N,N"
		return 	"\t".join([str(x) for x in [alignment[0], alignment[1], s_start, s_end, s_length, r_start, r_end, r_length, cut_left, cut_right, conv_str, weight, score]]);
	return alignment[0], alignment[1], s_start, s_end, s_length, r_start, r_end, r_length, cut_left, cut_right, conv, weight, score;		
	


	

#refseq_dict = {"hsa-let-7i" : "TGAGGTGAGTAGTTTGTGCTGTT", "hsa-let-7c" : "TGAGGTAGTAGGTTGTATGGTT", "stub" : "TGACGACGGAATA"}	
#seq = "TGAGGTAGTAGTTTGTACTGTTGCGC"

#g = get_best_alignments(refseq_dict, "@test", seq)
#for el in g:
	#print el
#print analyse_alignment(g[0])	
#print analyse_alignment(g[0], True)		
	