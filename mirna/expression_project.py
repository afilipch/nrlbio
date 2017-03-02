#! /usr/bin/python
'''Creates makefile and directory structure for mirna expression estimation from smal RNA sequencing experiment'''
import argparse
import sys
import os

from nrlbio.config.config import load_config;
from nrlbio.makefiles import dependence, get_header, get_bowtie_call, get_script, get_bowtie_help

#Load default parameters for bowtie2
conf = load_config('doublechiflex')

#Bowtie options preliminary parsing
bowtie_settings = conf['first']['bowtie']
bowtie_help_str = "[%s]" % " ".join(["%s%s=%s" % (x[1][1], x[0], x[1][0]) for x in bowtie_settings.items()])



parser = argparse.ArgumentParser(description='Creates makefile and directory structure for mirna expression estimation from smal RNA sequencing experiment');
#Required options
parser.add_argument('path', metavar = 'N', nargs = '?', type = os.path.abspath, help = "Path to the project folder. If folder does not exist, it will be created");
parser.add_argument('--reads', nargs = '?', type = os.path.abspath, required = True, help = "Path to sequencing reads. fastq/fasta file");
parser.add_argument('--index', nargs = '?', type = os.path.abspath, required = True, help = "Path to the mapping reference(miRNA sequences) bowtie2 index");
parser.add_argument('--mirbase_precursors', nargs = '?', required = True, type = os.path.abspath, help = "Path to the mirbase mirna precursors gff file");
parser.add_argument('--chiflex', nargs = '?', type = os.path.abspath, required = True, help = "Path to the Chiflex folder")
parser.add_argument('--name', nargs = '?', type = str, required = True, help = "Name of the project, will be used as a name for interactions ids")

parser.add_argument('--only_makefile', nargs = '?', default = False, const = True, type = bool, help = "If set, a new makefile is created, but not folders");
parser.add_argument('--nonunique', nargs = '?', default = False, const=True, type = bool, help = "If set, nonunique mappings with high alignment score are kept")
parser.add_argument('--collapsed', nargs = '?', default = False, const = True, type = str, help = "If set, reads are considered to be sequence-collapsed");
#bowtie2 options
parser.add_argument('--bowtie', nargs = '+', default = [], type = str, help = "Bowtie settings for the first round of mapping. For example, if one wants to set \'-p 4\', use \'--local\' alignment mode, but not \'--norc\' option then \'p=4 local=True norc=False\' should be provided. Given attributes replace default(for Chiflex, NOT for Bowtie) ones. Default settings for the modes are: %s" % bowtie_help_str)

args = parser.parse_args();

chiflex_package = os.path.join(args.chiflex, 'chiflex')
mirna_package = os.path.join(args.chiflex, 'mirna')







#######################################################################################################################
#Main function to create top level Makefile
def makefile_main():
	mlist=[];
	
	#Processing of the left chimeric part bowite2 settings
	bs_list = get_bowtie_call(bowtie_settings, args.bowtie, args.index, args.reads, args.name)
	
	#Map reads with bowtie2
	input_files = args.reads
	output_files = os.path.join('sam', '%s.mapped.sam' % args.name)
	script = bs_list
	mlist.append(dependence(input_files, output_files, script))
	

	
	#Collapse confident, but nonunique mappings into single sam entry, to protect them from further filtering out
	if(args.nonunique):
		input_files = output_files
		output_files = os.path.join('sam', 'collapsed.first.bam'), os.path.join('auxillary', 'collapsed.first.bed')
		script = get_script('collapse_nonunique_sam.py', arguments={'-s': output_files[0], '-b': output_files[1], '--minscore': 26}, inp = input_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
		#following reassignment is done for furthe consistency
		output_files = output_files[0]
		
	#Demultiplex mapping hits into control and true mapped reads
	input_files = output_files;
	output_files = [os.path.join('sam', '%s.%s.bam' % (args.name, x)) for x in ['unique', 'control']]
	script = get_script('demultiplex_sam.py', arguments={'--output': 'sam', '--name': args.name, '--score': conf['first']['demultiplex_sam']['score'], '--bestdistance': conf['first']['demultiplex_sam']['bestdistance']}, inp = input_files, package=chiflex_package)
	mlist.append(dependence(input_files, output_files, script));
	
	#Filter mapped reads using control ones
	input_files = output_files;
	output_files = os.path.join('sam', '%s.filtered.bam' % args.name)
	script = get_script('filter_alignment.py', arguments={'--output': output_files, '--features': " ".join(conf['first']['filter_alignment']['features']), '--fdr': conf['first']['filter_alignment']['fdr'], '--signal': input_files[0], '--control': input_files[1]}, package=chiflex_package)
	mlist.append(dependence(input_files, output_files, script));

	#Calculate raw expression
	input_files = output_files;
	output_files = os.path.join('sam', '%s.sorted.bam' % args.name)
	script = ('samtools', 'sort', '-f', input_files, output_files)
	mlist.append(dependence(input_files, output_files, script));
	
	input_files = output_files;
	output_files = os.path.join('expression', '%s.raw.gff' % args.name)
	if(args.collapsed):
		script = get_script('sam2expression.py', arguments={'--mirbase_precursors': args.mirbase_precursors, '--collapsed': True}, inp = input_files, out = output_files, package=mirna_package)
	else:
		script = get_script('sam2expression.py', arguments={'--mirbase_precursors': args.mirbase_precursors}, inp = input_files, out = output_files, package=mirna_package)
	mlist.append(dependence(input_files, output_files, script));	
	
	#Fix expression values for identical nonunique mappings  
	if(args.nonunique):
		input_files = output_files, os.path.join('auxillary', 'collapsed.first.bed');
		output_files = os.path.join('expression', '%s.adjusted.gff' % args.name)
		script = get_script('fix_expression_with_collapsed.py', arguments={'--collapsed': input_files[1]}, inp = input_files[0], out = output_files, package=mirna_package)
		mlist.append(dependence(input_files, output_files, script))

	
	
	#Get header and cleaner for the makefile
	mlist.insert(0, get_header(output_files))
	mlist.append('clean:\n\techo "nothing to clean."\n');
	
	return "\n\n".join(mlist)


#######################################################################################################################
#Create folders

project_path = args.path
folders = ['sam', 'reports', 'expression', 'log', 'auxillary']


while (not args.only_makefile):
	try:
		os.makedirs(project_path);
		for folder in folders:
			os.mkdir(os.path.join(project_path, folder));
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
('name', ('Project name', str)),
('path', ('Project folder', str)),
('reads', ('Sequencing reads', str)),
('chiflex', ('Chiflex module used in the project', str)), 
('index', ('Mapping reference index (genome, transcriptome) used for mapping', str)), 
 
('mirbase_precursors', ('Path to miRNAs precursors used for expression analysis', os.path.abspath)), 

('nonunique', ('Nonunique mappings with high alignment score are kept', str)),

('collapsed', ('Reads were supposed to be collapsed', str)),
)

with open(os.path.join(project_path, 'log/project.txt'), 'w') as rf:
	rf.write("Project call:\npython %s\n\n" % " ".join(sys.argv))
	
	for arg, (description, fun) in arguments_report: 
		av = getattr(args, arg);
		if(av):
			rf.write("%s:\t%s\n" % (description, fun(av)))
		else:
			rf.write("%s:\tnot set\n" % description)
			
	rf.write("bowtie2 settings:\n")
	for arg, (value, dashes) in bowtie_settings.items():
		rf.write("\t%s%s:\t%s\n" % (dashes, arg, value))


		
		
		