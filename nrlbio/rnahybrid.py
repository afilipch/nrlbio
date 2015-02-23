'''Library contains some wrappers/parsers for RNA hybridization tools'''
import commands;

### RNAhybrid wrappers

def _paired_rnahybrid(basepairing):
	'''Extracts pairing information from raw RNAhybrid text output 
		basepairing list: strings representing raw RNAhybrid text output 
	Returns list: list of values 0-> unbound nucleotide,1-> bound nucleotide along miRNA
	'''
	bound = [];
	micro_match = basepairing[2][::-1]; # 5->3 direction
	micro_unmatch = basepairing[3][::-1]; # 5->3 direction
	for i in range(len(micro_match)):
		if(micro_match[i] != " "):
			bound.append(1);
		elif(micro_unmatch[i] != " "):
			bound.append(0);
	return bound;  


def get_rnahybrid(target, miRNA, arguments = ""):
	"""Calculates hybridization energy of intermolecular interaction between given sequences via RNAhybrid tool. NOTE, RNAhybridis is designed to hybridize miRNA and respective target (not folding themselves). Therefore input sequences are called 'target' and 'miRNA' respectively, and basepairing of miRNA is also reported.
		target str: miRNAs target sequence
		miRNA str: miRNA sequence
		arguments str: optional arguments for RNAhybid. For detailed information see RNAhybid documentation
		
	Returns: 
		float: hybridization energy
		list: basepairing, 0-> unbound nucleotide,1-> bound nucleotide along miRNA
	"""
	call = "RNAhybrid -s 3utr_human \"" + target + "\"  \"" + miRNA + "\" -b 1 -c" + arguments + " ";
	status, output = commands.getstatusoutput(call)
	e = float(output.split(":")[4]);  
	basepairing = output.split(":")[7:11];  
	p = _paired_rnahybrid(basepairing);
	return e, p;
