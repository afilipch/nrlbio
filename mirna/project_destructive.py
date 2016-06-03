#! /usr/bin/python
'''Creates makefile and directory structure for destructive miRNA sites search'''
import argparse
import sys
import os


from nrlbio.makefiles import dependence, get_header, get_bowtie_call, get_script




parser = argparse.ArgumentParser(description='Creates makefile and directory structure for destructive miRNA sites search')

#Required files
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "Path to the project folder. If folder does not exist, it will be created");
seqerr = parser.add_argument('--sequences', nargs = '?', type = os.path.abspath, help = "Path to the sequences to look for destructive sites in, fasta format");
parser.add_argument('--interactions', nargs = '?', type = os.path.abspath, help = "Path to the chimera-based mirna:target interactions, double gff format. Is ignored if [--sequences] argument set");
parser.add_argument('--mir', nargs = '?', type = os.path.abspath, required = True, help = "Path to the miRNA sequences, fasta format");

#Constants for the project run
parser.add_argument("--system", nargs = '?', required=True, type=str, help="Model system/reference specie (hg19|dm6|...)")
parser.add_argument('--nrlbio', nargs = '?', type = os.path.abspath, required = True, help = "Path to the nrlbio folder")

#Conservation options
conserr = parser.add_argument("--maf", nargs = '?', type=os.path.abspath, help="Path to indexed MAF files. If set along with [--maf, --genomes, --mircons, --mir2ucsc], destructive conservation will be evaluated")
parser.add_argument("--genomes", nargs = '?', type=os.path.abspath, help="Path to genomes corresponding to the multiple specie alignment. If set along with [--maf, --genomes, --mircons, --mir2ucsc], destructive conservation will be evaluated")
parser.add_argument('--mircons', nargs = '?', type = os.path.abspath, help = "Path to the miRNAs conservation file, fasta format. If set along with [--maf, --genomes, --mircons, --mir2ucsc], destructive conservation will be evaluated");
parser.add_argument('--mir2ucsc', nargs = '?', type = os.path.abspath, help = "Path to the table which connects miRNA specie names with those from MAF files, tsv file. If set along with [--maf, --genomes, --mircons, --mir2ucsc], destructive conservation will be evaluated");

#Miscellaneous options
parser.add_argument('--annotation', nargs = '?', default =None, type = os.path.abspath, help = "Path to an annotation file in gff format. If provided, found genomic loci will be annotated");
parser.add_argument('--precursors', nargs = '?', default =None, type = os.path.abspath, help = "Path to the miRNA precursors, if set precursors will be removed from predictions, bed/gff format");
#parser.add_argument('--liftover', nargs = 2, default =None, type = os.path.abspath, help = "Path to the liftover files (1st file for conversion before conservation analyses, 2nd file for conservation after analyses), If set liftovers coordinates before conservation analyses and liftovers back after");
parser.add_argument('--threads', nargs = '?', default = 8, type = int, help = "Number of threads to use")
parser.add_argument('--minscore', nargs = '?', default = 20.0, type = float, help = "Only the regions with destructive score greater or equal to [--minscore] will be selected")
parser.add_argument('--reassign', nargs = '?', default = False, const=True, type = bool, help = "If set, regions position on a reference are reassigned to the genomic ones. Usefull in the case of nongenomic references(transcriptome, rRNAs, etc.). NOTE: reference headers have to be in [chrom]|[strand]|[start]|[stop] format")
parser.add_argument("--refgenome", nargs = '?', type=os.path.abspath, help="Path to genome. The option has to be set along with [--reassign]")

#Interface options
parser.add_argument('--only_makefile', nargs = '?', default = False, const = True, type = bool, help = "if set, a new makefile is created, but not folder structure");
args = parser.parse_args();

#######################################################################################################################
#Check input arguments
if(not (args.interactions or args.sequences)):
	raise argparse.ArgumentError(seqerr, "Value for the [--sequences] or [--interactions] argument has to be provided\n")

consargs = [args.maf, args.genomes, args.mircons, args.mir2ucsc]
if(any(consargs)):
	if(all(consargs)):
		checkcons = True
	else:	
		raise argparse.ArgumentError(conserr, "Values for the [--maf, --genomes, --mircons, --mir2ucsc] arguments have to be provided simultaneously or not provided at all\n")
else:
	checkcons = False
	
	
if(args.reassign):
	if(args.refgenome):
		genseq = args.refgenome
	else:	
		sys.exit('[--reassign] has to be accompanied with [--refgenome] option')
else:
	genseq = args.sequences;


#######################################################################################################################
#Set paths' constants

chiflex_package = os.path.abspath(os.path.join(args.nrlbio, 'chiflex'));
advanced_package = os.path.abspath(os.path.join(args.nrlbio, 'advanced'));
mirna_package = os.path.abspath(os.path.join(args.nrlbio, 'mirna'));
parsing_package = os.path.abspath(os.path.join(args.nrlbio, 'parsing', 'bed'));
scripts_package = os.path.abspath(os.path.join(args.nrlbio, 'scripts'));






#######################################################################################################################
#Main function to create top level Makefile
def makefile_main():
	mlist=[];
	clean = []
	output_files = args.interactions
	
	
	#Chimera independent prediction candidates preselection
	if(args.sequences):
		#Look for 2-8 seed matches with up to one mismatch in the provided sequences
		input_files = args.sequences
		output_files = 'anchors.bed'
		script = get_script('search_destructive_anchors.py', arguments={'--mir': args.mir}, inp = input_files, out = output_files,  package=mirna_package)
		mlist.append(dependence(input_files, output_files, script))
		clean.append(output_files)
	
		
		#Reassign to the genomic coordinates
		if(args.reassign):
			input_files = output_files
			output_files = 'coordassigned.bed'
			script = get_script('assign_coordinates.py', arguments={}, inp = input_files, out = output_files, package=chiflex_package)
			mlist.append(dependence(input_files, output_files, script))
			clean.append(output_files)
		
		
		#Extend regions to get sequences long enough for hybridization with miRNAs
		input_files = output_files
		output_files = 'extended.bed'
		script = get_script('extend_regions.py', arguments={'-l': 18, '-r': 2, '--strand': True}, inp = input_files, out = output_files,  package=parsing_package)
		mlist.append(dependence(input_files, output_files, script))
		clean.append(output_files)
		
	
		#Change format to double gff for further comaptibility
		input_files = output_files
		output_files = 'extended.gff'
		script = get_script('destructive2interactions.py', arguments={}, inp = input_files, out = output_files,  package=mirna_package)
		mlist.append(dependence(input_files, output_files, script))
		clean.append(output_files)
		
		
		#Assign sequences for miRNAs and their potential destructive sites
		input_files = output_files
		output_files = 'seqassigned.gff'
		script = get_script('assign_seq.py', arguments={'--fasta' :  [genseq, args.mir]}, inp = input_files, out = output_files,  package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
		clean.append(output_files)
		
	else:
		pass;
		
	#Assign destructive score for each miRNA:target interaction, filter out candidates  with a low score
	input_files = output_files
	output_files = 'dscored.gff'
	script = get_script('assign_destructive_score.py', arguments={'--mir' : args.mir, '--ftype': 'double', '--system': args.system, '--threads': args.threads, '--minscore': args.minscore}, inp = input_files, out = output_files,  package=mirna_package)
	mlist.append(dependence(input_files, output_files, script))
	
	if(args.precursors):
		#Assign destructive score for each miRNA:target interaction, filter out candidates  with a low score
		input_files = output_files
		output_files = 'dscored_wo_precursors.gff'
		script = ('bedtools', 'intersect', '-a', input_files, '-b', args.precursors, '-s', '-v', '>', output_files)
		mlist.append(dependence(input_files, output_files, script))
	
	
	if(checkcons):
		#Extract MAF blocks for potential destructive sites
		input_files = output_files
		output_files = 'conserved.fa'
		script = get_script('extract_maf.py', arguments={'--maf': args.maf, '--genomes': args.genomes, '--table': args.mir2ucsc, '--system': args.system, '--left': 5, '--right': 5, '--min-len': 18, '--muscle': True}, inp = input_files, out = output_files, package = scripts_package)
		mlist.append(dependence(input_files, output_files, script))
		clean.append(output_files)
		
		#Assign conservation destructive score based on previously extracted MAF blocks
		input_files = 'dscored.gff', 'conserved.fa'
		output_files = 'cons_assigned.gff'
		script = get_script('destructive_conservation.py', arguments={'--mir' : args.mircons, '--mir2ucsc': args.mir2ucsc, '--refspecie': args.system, '--mafasta': input_files[1]}, inp = input_files[0], out = output_files,  package=mirna_package)
		mlist.append(dependence(input_files, output_files, script))
		
	if(args.annotation):
		#Create final HTML report on destructive sites found
		input_files = output_files
		output_files = 'annotated.gff'
		script = get_script('annotate_bed_with_gff3.py', arguments={'--gff3' : args.annotation}, inp = input_files, out = output_files,  package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
		
		
	#Create final HTML report on destructive sites found
	input_files = output_files
	output_files = 'destructive.gff', 'destructive.html'
	script = get_script('destructive_report.py', arguments={'--output' : output_files[1], '--system': args.system, '--best_only': 1000}, inp = input_files, out = output_files[0],  package=mirna_package)
	mlist.append(dependence(input_files, output_files, script))


		
	#Get header and cleaner for the makefile
	mlist.insert(0, get_header(output_files))
	mlist.append('.PHONY: clean\nclean:\n\trm %s\n' % " ".join(clean));
	
	return "\n\n".join(mlist)


#######################################################################################################################
#Create folder

project_path = os.path.abspath(args.path)

while (not args.only_makefile):
	try:
		os.makedirs(project_path);
		break;
	except:
		answer = raw_input("\nProject directory \'%s\' is currently exists, please type 'N' if you don't want to create a new project, 'MO' if you want to change/create only the Makefile, [another project name] if you want to create a new folder structure and makefile: " % project_path)
		if(answer=='N'):
			sys.exit('Project was not created')
		elif(answer=='MO'):
			sys.stderr.write('Makefile was changed/added to the existing project %s\n' % project_path)
			break
		else:
			project_path = os.path.abspath(answer)
			
			
#######################################################################################################################
#Create Makefile
with open(os.path.join(project_path, 'Makefile'), 'w') as mf:
	mf.write(makefile_main());



#######################################################################################################################
#Create a report

def multipath(l):
	return "\t".join([os.path.abspath(x) for x in l])


arguments_report = (
('path', ('Project folder', os.path.abspath)),
('sequences', ('Sequences to search for destructive sites in', os.path.abspath)),
('interactions', ('miRNA:target Interactions to search for destructive sites', os.path.abspath)), 
('mir', ('miRNAs sequences', os.path.abspath)), 

('maf', ('Maf blocks used for conservation analysis', os.path.abspath)),
('genomes', ('Genomes used for conservation analysis', os.path.abspath)), 
('mircons', ('miRNAs\' conservation file used for conservation analysis', os.path.abspath)), 
('mir2ucsc', ('MirBase to UCSC translation table used for conservation analysis', os.path.abspath)), 

('system', ('Genomic system used', str)), 
('annotation', ('Annotation system used for interactions annotation', os.path.abspath)), 
('only_makefile', ('New Makefile was generated', str)),  
)

with open(os.path.join(project_path, 'call.txt'), 'w') as rf:
	rf.write("Project call:\npython %s\n\n" % " ".join(sys.argv))
	
	for arg, (description, fun) in arguments_report: 
		av = getattr(args, arg);
		if(av):
			rf.write("%s:\t%s\n" % (description, fun(av)))
		else:
			rf.write("%s:\tnot set\n" % description)
			
		
		
		