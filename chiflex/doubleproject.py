#! /usr/bin/python
'''Creates makefile and directory structure for chiflex double reference project'''
import argparse
import sys
import os

from nrlbio.config.config import load_config;
from nrlbio.makefiles import dependence, get_header, get_bowtie_call, get_script, get_bowtie_help

conf = load_config('doublechiflex')

#Bowtie options preliminary parsing
bowtie_settings1 = conf['first']['bowtie'];
bowtie_settings2 = conf['second']['bowtie'];

bowtie_help_str1 = "[%s]" % " ".join(["%s%s=%s" % (x[1][1], x[0], x[1][0]) for x in bowtie_settings1.items()])
bowtie_help_str2 = "[%s]" % " ".join(["%s%s=%s" % (x[1][1], x[0], x[1][0]) for x in bowtie_settings2.items()])


parser = argparse.ArgumentParser(description='Creates makefile and directory structure for chiflex project')#, formatter_class = argparse.RawTextHelpFormatter);
#Required options
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the project folder. If folder does not exist, it will be created");
parser.add_argument('--reads', nargs = '?', type = str, required = True, help = "path to sequencing reads. fastq/fasta file");
parser.add_argument('--indices', nargs = 2, type = str, required = True, help = "path to the mapping reference bowtie2 indices. NOTE: At the first step, mapping will be done for the first index provided");
parser.add_argument('--chiflex', nargs = '?', type = str, required = True, help = "path to the Chiflex folder")
parser.add_argument('--name', nargs = '?', type = str, required = True, help = "name of the project, will be used as name for interactions ids")


#Paths to the files for additional annotation;
parser.add_argument('--references', nargs = '+', type = str, help = "path to the references in fasta format. If set, sequences will be assigned to resulting interactions as well as energy and binding pattern(pattern of paired nucleotides on a left chimeric part)");
parser.add_argument('--annotation', nargs = '?', type = str, help = "path to an annotation file in gff format. If provided, found genomic loci will be annotated");
parser.add_argument('--mirna', nargs = '?', type = str, help = "path miRNA sequences in fasta format. If set interactions will be further analized as miRNA:target chimeras. NOTE: This mode is applicable for any small RNA");

#Options for the mapping result postprocessing
parser.add_argument('--repetitive', nargs = '?', default = False, const=True, type = bool, help = "if set, repetitive mapped sequences are removed")
parser.add_argument('--nonunique', nargs = '?', default = False, const=True, type = bool, help = "if set, nonunique mappings with high alignment score are kept")
parser.add_argument('--reassign', nargs = '?', default = False, const=True, type = bool, help = "if set, hits position on a reference will be reassigned to the genomic ones. Usefull in the case of nongenomic references(transcriptome, rRNAs, etc.). NOTE: reference headers has to be in [chrom]|[strand]|[start]|[stop] format")

#Options for the output control
parser.add_argument('--only_makefile', nargs = '?', default = False, const = True, type = bool, help = "if set, a new makefile is created, but not folder structure");

#bowtie2 options
parser.add_argument('--bowtie_args1', nargs = '+', default = [], type = str, help = "Bowtie settings for the first round of mapping. For example, if one wants to set \'-p 4\', use \'--local\' alignment mode, but not \'--norc\' option then \'p=4 local=True norc=False\' should be provided. Given attributes replace default(for Chiflex, NOT for Bowtie) ones. Default settings for the modes are: %s" % bowtie_help_str1)
parser.add_argument('--bowtie_args2', nargs = '+', default = [], type = str, help = "Bowtie settings for the second round of mapping. For example, if one wants to set \'-p 4\', use \'--local\' alignment mode, but not \'--norc\' option then \'p=4 local=True norc=False\' should be provided. Given attributes replace default(for Chiflex, NOT for Bowtie) ones. Default settings for the modes are: %s" % bowtie_help_str2)
args = parser.parse_args();

#######################################################################################################################
##Set some constants
firstname = "%s.%s" % (args.name, "first")
secondname = "%s.%s" % (args.name, "second")

chiflex_package = os.path.abspath(os.path.join(args.chiflex, 'chiflex'));
advanced_package = os.path.abspath(os.path.join(args.chiflex, 'advanced'));
mirna_package = os.path.abspath(os.path.join(args.chiflex, 'mirna'));







#######################################################################################################################
#Function to create top level Makefile
def makefile_main():
	mlist=[];
	
	#processing of the right chimeric part
	bs_list1 = get_bowtie_call(bowtie_settings1, args.bowtie_args1, args.indices[0], args.reads, firstname)
	#Map reads with bowtie2
	input_files = os.path.abspath(args.reads)
	output_files = os.path.join('sam', '%s.mapped.sam' % firstname)
	script = bs_list1
	mlist.append(dependence(input_files, output_files, script))
	
	#Remove repetitive mappings if option 'repetitive' set
	if(args.repetitive):
		input_files = output_files
		output_files = os.path.join('sam', 'nonrep.first.bam')
		script = get_script('filter_sam.py', arguments={'--output': output_files, '--filters': 'repetitive', '--arguments': 'min_entropy=1.6'}, inp = input_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
	
	#Collapse confident, but nonunique mappings into single sam entry, to protect them from further filtering out
	if(args.nonunique):
		input_files = output_files
		output_files = os.path.join('sam', 'collapsed.first.bam'), os.path.join('auxillary', 'collapsed.first.bed')
		script = get_script('collapse_nonunique_sam.py', arguments={'-s': output_files[0], '-b': output_files[1], '--minscore': 26}, inp = input_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
		#following reassignment is done for furthe consistency
		output_files = output_files[0]
		
		
	input_files = output_files;
	output_files = [os.path.join('sam', '%s.%s.bam' % (firstname, x)) for x in ['unique', 'control']]
	script = get_script('demultiplex_sam.py', arguments={'--output': 'sam', '--name': firstname, '--score': conf['first']['demultiplex_sam']['score'], '--bestdistance': conf['first']['demultiplex_sam']['bestdistance']}, inp = input_files, package=chiflex_package)
	mlist.append(dependence(input_files, output_files, script));
	
	input_files = output_files;
	output_files = os.path.join('sam', '%s.filtered.bam' % firstname)
	script = get_script('filter_alignment.py', arguments={'--output': output_files, '--features': " ".join(conf['first']['filter_alignment']['features']), '--fdr': conf['first']['filter_alignment']['fdr'], '--signal': input_files[0], '--control': input_files[1]}, package=chiflex_package)
	mlist.append(dependence(input_files, output_files, script));
	
	input_files = output_files;
	output_files = os.path.join('sam', '%s.length%d.bam' % (firstname, conf['min_right_length'])) 
	script = get_script('filter_sam.py', arguments={'--output': output_files, '--filters': 'chimera_right', '--arguments': "minright=%d" % conf['min_right_length']}, inp = input_files, package=chiflex_package)
	mlist.append(dependence(input_files, output_files, script));

	input_files = output_files;
	output_files = os.path.join('fastq', '%s.length%d.fastq' % (firstname, conf['min_right_length'])) 
	script = get_script('sam2fastq_masked.py', arguments={'--pad': '\'\''}, inp = input_files, out=output_files, package=advanced_package)
	mlist.append(dependence(input_files, output_files, script));
	
	
	#processing of the right chimeric part
	bs_list2 = get_bowtie_call(bowtie_settings2, args.bowtie_args2, args.indices[1], os.path.join(args.path, output_files), secondname)#reads will come further	
	#Map reads with bowtie2
	input_files = output_files
	output_files = os.path.join('sam', '%s.mapped.sam' % secondname)
	script = bs_list2
	mlist.append(dependence(input_files, output_files, script))
	
	#Remove repetitive mappings if option 'repetitive' set
	if(args.repetitive):
		input_files = output_files
		output_files = os.path.join('sam', 'nonrep.second.bam')
		script = get_script('filter_sam.py', arguments={'--output': output_files, '--filters': 'repetitive', '--arguments': 'min_entropy=1.4'}, inp = input_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
	
	#Collapse confident, but nonunique mappings into single sam entry, to protect them from further filtering out
	if(args.nonunique):
		input_files = output_files
		output_files = os.path.join('sam', 'collapsed.second.bam'), os.path.join('auxillary', 'collapsed.second.bed')
		script = get_script('collapse_nonunique_sam.py', arguments={'-s': output_files[0], '-b': output_files[1], '--minscore': 36}, inp = input_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
		#following reassignment is done for furthe consistency
		output_files = output_files[0]
		
	input_files = output_files;
	output_files = [os.path.join('sam', '%s.%s.bam' % (secondname, x)) for x in ['unique', 'control']]
	script = get_script('demultiplex_sam.py', arguments={'--output': 'sam', '--name': secondname, '--score': conf['second']['demultiplex_sam']['score'], '--bestdistance': conf['second']['demultiplex_sam']['bestdistance']}, inp = input_files, package=chiflex_package)
	mlist.append(dependence(input_files, output_files, script));
	
	input_files = output_files;
	output_files = os.path.join('sam', '%s.filtered.bam' % secondname)
	script = get_script('filter_alignment.py', arguments={'--output': output_files, '--features': " ".join(conf['second']['filter_alignment']['features']), '--fdr': conf['second']['filter_alignment']['fdr'], '--signal': input_files[0], '--control': input_files[1]}, package=chiflex_package)
	mlist.append(dependence(input_files, output_files, script));
	
	
	input_files = os.path.join('sam', '%s.length%d.bam' % (firstname, conf['min_right_length'])), output_files
	output_files = os.path.join('chimeras', '%s.merged.gff' % args.name)
	script = get_script('merged2chimeras.py', arguments={'--oformat': 'gff', '--fixgap': ""}, inp=input_files, out=output_files, package=chiflex_package)
	mlist.append(dependence(input_files, output_files, script));
	
	if(args.reassign):
		input_files = output_files 
		output_files = os.path.join('chimeras', '%s.coord_assigned.gff' % args.name)
		script = get_script('assign_coordinates.py', arguments={'--quite': ''}, inp = input_files, out = output_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
		
	input_files =  output_files
	output_files = os.path.join('interactions', '%s.sorted.gff' % args.name) 
	script = get_script('sort.py', inp=input_files, out = output_files, package=chiflex_package)
	mlist.append(dependence(input_files, output_files, script));

	input_files =  output_files
	output_files = os.path.join('interactions', '%s.interactions.gff' % args.name),  os.path.join('interactions', '%s.rid2iid.tsv' % args.name)
	script = get_script('collapse2interaction.py', arguments={'-od': output_files[1], '--name': args.name, '--distance': conf['collapse2interaction']['distance']}, inp=input_files, out = output_files[0], package=chiflex_package)
	mlist.append(dependence(input_files, output_files, script));
	output_files = output_files[0];
	
	if(args.references):
		input_files = output_files
		output_files = os.path.join('interactions', '%s.seqassigned.gff' % args.name)
		script = get_script('assign_seq.py', arguments={'--fasta': " ".join([os.path.abspath(x) for x in args.references])}, inp=input_files, out = output_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script));
		
		input_files =  output_files
		output_files = os.path.join('interactions', '%s.energyassigned.gff' % args.name)
		if(args.mirna):
			script = get_script('assign_energy.py', arguments={'--shuffle_trials': conf['shuffle_trials'], '--pattern': ''}, inp=input_files, out=output_files, package=chiflex_package)
		else:
			script = get_script('assign_energy.py', arguments={'--shuffle_trials': conf['shuffle_trials']}, inp=input_files, out=output_files,  package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script));
		
		if(args.mirna):
			input_files =  output_files
			output_files = os.path.join('interactions', '%s.modeassigned.gff' % args.name)
			script = get_script('assign_mode.py', arguments={'--mir': os.path.abspath(args.mirna), '--set_control': ''}, inp=input_files, out = output_files, package=mirna_package)
			mlist.append(dependence(input_files, output_files, script));
		
	mlist.insert(0, get_header(output_files))
	# makefie cleaner
	mlist.append('clean:\n\techo "nothing to clean."\n');
	return "\n\n".join(mlist)


#######################################################################################################################
#Create folder structure
project_path = os.path.abspath(args.path)
folders = ['sam', 'reports', 'chimeras', 'log', 'auxillary', 'interactions', 'fastq']


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
#Create Makefiles
with open(os.path.join(project_path, 'Makefile'), 'w') as mf:
	mf.write(makefile_main());



def multipath(l):
	return "\t".join([os.path.abspath(x) for x in l])

#report project call:
arguments_report = (
('name', ('Project name, assigned to the generated interactions', str)),
('path', ('Project folder', os.path.abspath)),
('reads', ('Sequencing reads', os.path.abspath)),
('chiflex', ('Chiflex module used in the project', os.path.abspath)), 
('indices', ('Mapping reference indices (genome, transcriptome) used for mapping', multipath)), 
('references', ('Mapping reference sequences (genome, transcriptome) used for mapping', multipath)), 
('mirna', ('Path to miRNAs used for downstream analysis', os.path.abspath)), 

('annotation', ('Annotation system used for interactions annotation', os.path.abspath)), 
('repetitive', ('Repetitive mapped sequences are removed', str)),
('nonunique', ('Nonunique mappings with high alignment score are kept', str)),
('only_makefile', ('New Makefile was generated', str)),  
)

with open(os.path.join(project_path, 'log/project.txt'), 'w') as rf:
	rf.write("Project call:\npython %s\n\n" % " ".join(sys.argv))
	
	for arg, (description, fun) in arguments_report: 
		av = getattr(args, arg);
		if(av):
			rf.write("%s:\t%s\n" % (description, fun(av)))
		else:
			rf.write("%s:\tnot set\n" % description)
			
	rf.write("bowtie2 settings for first round of mapping:\n")
	for arg, (value, dashes) in bowtie_settings1.items():
		rf.write("\t%s%s:\t%s\n" % (dashes, arg, value))

	rf.write("bowtie2 settings for second round of mapping:\n")
	for arg, (value, dashes) in bowtie_settings2.items():
		rf.write("\t%s%s:\t%s\n" % (dashes, arg, value))