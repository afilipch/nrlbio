#! /usr/bin/python	
import commands;

def paired(sl):
	'''input: list of strings representing hybrid. output: list of values 0-> unbound nucleotide,1-> bound nucleotide along miRNA'''
	bound = [];
	micro_match = sl[2][::-1]; # 5->3 direction
	micro_unmatch = sl[3][::-1]; # 5->3 direction
	for i in range(len(micro_match)):
		if(micro_match[i] != " "):
			bound.append(1);
		elif(micro_unmatch[i] != " "):
			bound.append(0);
	return bound;        
  
def run(target, micro, arguments = ""):
	call = "RNAhybrid -s 3utr_human \"" + target + "\"  \"" + micro + "\" -b 1 -c" + arguments + " ";
	status, output = commands.getstatusoutput(call)
	e = float(output.split(":")[4]);  
	sl = output.split(":")[7:11];  
	p = paired(sl);
	return p, e;
 
def represent(target, micro, arguments = ""):
	call = "RNAhybrid -s 3utr_human \"" + target + "\"  \"" + micro + "\" -b 1 " + arguments + " ";
	status, output = commands.getstatusoutput(call)
	t = output.split("\n")
	e = float(t[5].split()[1]);
	return "\n".join(t[3:]), e; 