#! /usr/bin/python
'''Creates makefile and directory structure for chiflex project'''
import argparse
import sys
import os

from nrlbio.config.config import load_config;
from nrlbio.makefiles import dependence, get_header, get_bowtie_call, get_script, get_bowtie_help

conf = load_config('chiflex')

#Modes order preliminary parsing
modes_order = conf['modes_order'];

#set possible chimera types
interaction_types = ['inter', 'intra', 'csj', 'lsj']


#Bowtie options preliminary parsing
bowtie_configurations = conf['bowtie'];
bowtie_help_str = get_bowtie_help(bowtie_configurations)
	

parser = argparse.ArgumentParser(description='Creates makefile and directory structure for chiflex project')#, formatter_class = argparse.RawTextHelpFormatter);

#Required options
parser.add_argument('path', metavar = 'N', nargs = '?', type = str, help = "path to the project folder. If folder does not exist, it will be created");
parser.add_argument('--reads', nargs = '+', type = os.path.abspath, required = True, help = "path to sequencing reads. fastq/fasta file. If more than one file provided, reads will be collapsed");
parser.add_argument('--reference', nargs = '?', type = os.path.abspath, required = True, help = "path to the mapping reference");
parser.add_argument('--chiflex', nargs = '?', type = os.path.abspath, required = True, help = "path to the Chiflex folder")
parser.add_argument('--name', nargs = '?', type = str, required = True, help = "name of the project, will be used as name for interactions ids")
arg_modes = parser.add_argument('--modes', nargs = '+', type = str, required = True, choices = modes_order.keys(),  help = "Mode of the Chiflex. If you are looking for circular or linear splice junctions - choose mode \'splicing\'. If you want to find read clusters(binding sites) - choose mode \'clustering\'. If you explore RNA-RNA intra- or intermolecular interactions - choose mode \'interaction\'. NOTE: it is possible to run more than one mode together. However mapping options will be selected for the mode with a highest(1st is higher than 2nd) priority. Priorities: %s. NOTE: Without setting [--genome] or [--exon] option interaction and splicing modes cannot be used. NOTE: Without setting [--genome_index] option clustering mode cannot be used" % ["%s=%d" % x for x in modes_order.items()])

#Paths to the files for additional annotation;
parser.add_argument('--genome', nargs = '?', type = os.path.abspath, default = None, help = "path to a reference sequences(genome, transcriptome, etc.) in fasta format. Is required for de novo splice sites search and downstream sequence analysis. NOTE: if [--genome] and [--exons] are provided simultaneously, then splice sites will be de novo annotated based on [--genome]");
parser.add_argument('--exons', nargs = '?', type = os.path.abspath, default = None, help = "path to a file of exonic regions in bed format. If provided, specific type(circular or linear splice junction, intra- or intermolecular interaction) will be assigned to each interaction. That is, only known splice sites will be found. NOTE: if [--genome] and [--exons] are provided simultaneously, then splice sites will be de novo annotated based on [--genome]");
parser.add_argument('--annotation', nargs = '?', type = os.path.abspath, default = None, help = "path to an annotation file in gff format. If provided, found genomic loci will be annotated");
parser.add_argument('--genome_index', nargs = '?', type = os.path.abspath, default = None, help = "reference(genome) index file (.fai)");

#Options for the mapping result postprocessing
parser.add_argument('--collapsed', nargs = '?', default = False, const=True, type = bool, help = "If set, reads are assumed collapsed with collpse.pl script. Read count appendix of the read id will be used to calculate read support of the interactions")
parser.add_argument('--repetitive', nargs = '?', default = False, const=True, type = bool, help = "if set, repetitive mapped sequences are removed")
parser.add_argument('--nonunique', nargs = '?', default = False, const=True, type = bool, help = "if set, nonunique mappings with high alignment score are kept")
parser.add_argument('--reassign', nargs = '?', default = False, const=True, type = bool, help = "if set, hits position on a reference will be reassigned to the genomic ones. Usefull in the case of nongenomic references(transcriptome, rRNAs, etc.). NOTE: reference headers has to be in [chrom]|[strand]|[start]|[stop] format")
parser.add_argument('--stranded', nargs = '?', default = False, const=True, type = bool, help = "Should be set if the sequencing data are stranded")

#Options for the output control
parser.add_argument('--only_makefile', nargs = '?', default = False, const = True, type = bool, help = "if set, a new makefile is created, but not folder structure");
parser.add_argument('--reports', nargs = '?', default = False, const = True, type = bool, help = "if set, html reports will be produced");

#bowtie2 options
parser.add_argument('--bowtie', nargs = '+', default = [], type = str, help = "Bowtie settings. For example, if one wants to set \'-p 4\', use \'--local\' alignment mode, but not \'--norc\' option then \'p=4 local=True norc=False\' should be provided. Given attributes replace default(for Chiflex, NOT for Bowtie) ones. Default settings for the modes are:%s" % bowtie_help_str)
args = parser.parse_args();

print args.reads

#######################################################################################################################
#Set priotity for the modes
main_mode = min(args.modes, key= lambda x: modes_order[x]);


#######################################################################################################################
##Parse bowtie options
bowtie_settings = bowtie_configurations[main_mode];


#######################################################################################################################
#Configure input arguments
if('interaction' in args.modes or 'splicing' in args.modes):
	if(args.genome):
		annotate_with_genome = True;
		genome_path = args.genome
	elif(args.exons):
		annotate_with_genome = False;
	else:
		raise argparse.ArgumentError(arg_modes, "Without setting [--genome] or [--exon] option interaction and splicing modes cannot be used\n")
		
if(('clustering' in args.modes) and (not args.genome_index)):
	raise argparse.ArgumentError(arg_modes, "Without setting [--genome_index] option clustering mode cannot be used\n")
	
	
chiflex_package = os.path.join(args.chiflex, 'chiflex')
splicing_package = os.path.join(args.chiflex, 'splicing')
interaction_package = os.path.join(args.chiflex, 'interaction')




#######################################################################################################################
#Function to create top level Makefile
def makefile_main():
	clean = []
	mlist=[];
	generated_makefiles = [];
	modes_makefiles = [];
	
	#Collapse reads, if more than one fasta/fastq file is provided
	if(len(args.reads) > 1):
		input_files = args.reads
		output_files = 'merged.fastq'
		script = 'cat', " ".join(input_files), '>', output_files
		mlist.append(dependence(input_files, output_files, script))
		clean.append(output_files)
		
		
		input_files = output_files
		output_files = 'collapsed.fastq'
		script = 'collapse.pl', input_files, '>', output_files
		mlist.append(dependence(input_files, output_files, script))
		clean.append(output_files)
		
		collapsed = True
		
	else:
		input_files = args.reads
		collapsed = args.collapsed
	
	#Map reads with bowtie2
	input_files = output_files
	output_files = os.path.join('sam', '%s.mapped.sam' % args.name)
	script = get_bowtie_call(bowtie_settings, args.bowtie, args.reference, input_files, args.name)
	mlist.append(dependence(input_files, output_files, script))
	clean.append(output_files)
	
	#Remove repetitive mappings if option 'repetitive' set
	if(args.repetitive):
		input_files = output_files
		output_files = os.path.join('sam', 'nonrep.bam')
		script = get_script('filter_sam.py', arguments={'--output': output_files, '--filters': 'repetitive', '--arguments': 'min_entropy=1.6'}, inp = input_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
		clean.append(output_files)
	
	#Collapse confident, but nonunique mappings into single sam entry, to protect them from further filtering out
	if(args.nonunique):
		input_files = output_files
		output_files = os.path.join('sam', 'collapsed.bam'), os.path.join('auxillary', 'collapsed.bed')
		script = get_script('collapse_nonunique_sam.py', arguments={'-s': output_files[0], '-b': output_files[1], '--minscore': 62}, inp = input_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
		#following reassignment is done for furthe consistency
		output_files = output_files[0]
		clean.append(output_files)
		
		
		
	#If only the clustering mode has been set, then script switches to the makefile_clustering at this point. For other cases switch happens later
	if(main_mode=='clustering'):
		input_files = output_files;
		output_files = [os.path.join('sam', '%s.%s.bam' % (args.name, x)) for x in ['unique', 'control']]
		script = get_script('demultiplex_sam.py', arguments={'--output': 'sam', '--name': args.name, '--score': conf['demultiplex_sam']['score'], '--bestdistance': conf['demultiplex_sam']['bestdistance']}, inp = input_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script));
		
		input_files = ['makefile_clustering'] + output_files
		output_files = 'clustering'
		modes_makefiles.append('clustering')
		script = ["$(MAKE)", '-f', '$<']
		mlist.append(dependence(input_files, output_files, script))		
		generated_makefiles.append('makefile_clustering');
		
		
	else:
		#demultiplex mapping hits into single and chimeric reads
		input_files = output_files
		output_files = [os.path.join('sam', '%s.%s.bam' % (args.name, x)) for x in ['unique', 'unique_chimera', 'control_chimera', 'control']]
		script = get_script('demultiplex_chimera.py', arguments={'--output': 'sam', '--name': args.name, '--score': conf['demultiplex_chimera']['score'], '--score_chimera': conf['demultiplex_chimera']['score_chimera'], '--maxgap': conf['demultiplex_chimera']['maxgap'], '--s_distance': conf['demultiplex_chimera']['s_distance'], '--ch_distance': conf['demultiplex_chimera']['ch_distance']}, inp = input_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
			
		#Merge sam hits into chimeras in doublebed format
		input_files = os.path.join('sam', '%s.unique_chimera.bam' % args.name) 
		output_files = os.path.join('chimeras', 'unique.bed') 
		script = get_script('merged2chimeras.py', arguments={}, inp = input_files, out = output_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
		
		#Merge sam hits into chimeras in doublebed format for decoy(control) reference
		input_files = os.path.join('sam', '%s.control_chimera.bam' % args.name) 
		output_files = os.path.join('chimeras', 'control.bed') 
		script = get_script('merged2chimeras.py', arguments={}, inp = input_files, out = output_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
		
		if(args.reassign):
			input_files = os.path.join('chimeras', 'control.bed') 
			output_files = os.path.join('chimeras', 'control.assigned.bed') 
			script = get_script('assign_coordinates.py', arguments={}, inp = input_files, out = output_files, package=chiflex_package)
			mlist.append(dependence(input_files, output_files, script))
			
			input_files = os.path.join('chimeras', 'unique.bed') 
			output_files = os.path.join('chimeras', 'unique.assigned.bed') 
			script = get_script('assign_coordinates.py', arguments={}, inp = input_files, out = output_files, package=chiflex_package)
			mlist.append(dependence(input_files, output_files, script));
			
			controlf = os.path.join('chimeras', 'control.assigned.bed')
			uniquef = os.path.join('chimeras', 'unique.assigned.bed') 
		else:
			controlf = os.path.join('chimeras', 'control.bed');
			uniquef = os.path.join('chimeras', 'unique.bed');
			
			
		if(annotate_with_genome):
			input_files = uniquef
			output_files = os.path.join('chimeras', 'unique.annotated.gff')
			script = get_script('annotate_novel.py', arguments={'--reference': genome_path}, inp = input_files, out=output_files, package=chiflex_package)
			mlist.append(dependence(input_files, output_files, script));
		
			input_files = output_files
			output_files = [os.path.join('chimeras', 'unique.annotated.%s.gff' % x) for x in interaction_types]
			script = get_script('stratify_gff.py', arguments={'--attributes': 'ntype', '--output': 'chimeras', '--rtypes': " ".join(interaction_types)}, inp = input_files, package=chiflex_package)
			mlist.append(dependence(input_files, output_files, script));
			
			
			input_files = controlf
			output_files = os.path.join('chimeras', 'control.annotated.gff')
			script = get_script('annotate_novel.py', arguments={'--reference': genome_path}, inp = input_files, out=output_files, package=chiflex_package)
			mlist.append(dependence(input_files, output_files, script));
		
			input_files = output_files
			output_files = [os.path.join('chimeras', 'control.annotated.%s.gff' % x) for x in interaction_types]
			script = get_script('stratify_gff.py', arguments={'--attributes': 'ntype', '--output': 'chimeras', '--rtypes': " ".join(interaction_types)}, inp = input_files, package=chiflex_package)
			mlist.append(dependence(input_files, output_files, script));
			
		else:
			input_files = uniquef, args.exons
			output_files = os.path.join('chimeras', 'unique.annotated.gff')
			if(args.stranded):
				cargs = arguments={'--exons': input_files[1], '--stranded': ''}
			else:
				cargs = arguments={'--exons': input_files[1]}
			script = get_script('annotate_chimera.py', arguments=cargs, inp = input_files[0], out=output_files, package=chiflex_package)
			mlist.append(dependence(input_files, output_files, script));
			
			input_files = output_files
			output_files = [os.path.join('chimeras', 'unique.annotated.%s.gff' % x) for x in interaction_types]
			script = get_script('stratify_gff.py', arguments={'--attributes': 'ktype', '--output': 'chimeras', '--rtypes': " ".join(interaction_types)}, inp = input_files, package=chiflex_package)
			mlist.append(dependence(input_files, output_files, script));			
			
	
			input_files = controlf, args.exons
			output_files = os.path.join('chimeras', 'control.annotated.gff')
			if(args.stranded):
				cargs = arguments={'--exons': input_files[1], '--stranded': ''}
			else:
				cargs = arguments={'--exons': input_files[1]}
			script = get_script('annotate_chimera.py', arguments=cargs, inp = input_files[0], out=output_files, package=chiflex_package)
			mlist.append(dependence(input_files, output_files, script));
			
			input_files = output_files
			output_files = [os.path.join('chimeras', 'control.annotated.%s.gff' % x) for x in interaction_types]
			script = get_script('stratify_gff.py', arguments={'--attributes': 'ktype', '--output': 'chimeras', '--rtypes': " ".join(interaction_types)}, inp = input_files, package=chiflex_package)
			mlist.append(dependence(input_files, output_files, script));	
			
			
		
		if('clustering' in args.modes):
			input_files = ['makefile_clustering'] + [os.path.join('sam', '%s.%s.bam' % (args.name, x)) for x in ['unique', 'control']]
			modes_makefiles.append('clustering')
			output_files = 'clustering'
			script = ["$(MAKE)", '-f', '$<']
			mlist.append(dependence(input_files, output_files, script))
			generated_makefiles.append('makefile_clustering');
			
		if('interaction' in args.modes):
			input_files = ['makefile_interaction'] + [os.path.join('chimeras', 'control.annotated.%s.gff' % x) for x in ['inter', 'intra']] + [os.path.join('chimeras', 'unique.annotated.%s.gff' % x) for x in ['inter', 'intra']]
			modes_makefiles.append('interaction')
			output_files = 'interaction'
			script = ["$(MAKE)", '-f', '$<']
			mlist.append(dependence(input_files, output_files, script))
			generated_makefiles.append('makefile_interaction');
			
		if('splicing' in args.modes):
			input_files = ['makefile_splicing'] + [os.path.join('chimeras', 'control.annotated.%s.gff' % x) for x in ['lsj', 'csj']] + [os.path.join('chimeras', 'unique.annotated.%s.gff' % x) for x in ['lsj', 'csj']]
			modes_makefiles.append('splicing');
			output_files = 'splicing'
			script = ["$(MAKE)", '-f', '$<']
			mlist.append(dependence(input_files, output_files, script))
			generated_makefiles.append('makefile_splicing');		
		
		
	#makefile header
	mlist.insert(0, get_header(modes_makefiles, phony=True))
	# makefie cleaner
	mlist.append( 'clean:\n%s\n\trm %s' %  ("\n".join(["\t$(MAKE) -f %s clean" % x for x in generated_makefiles]), ' '.join(clean)) );

	return "\n\n".join(mlist)


		
		
#######################################################################################################################
#Function to create clustering Makefile
def makefile_clustering():
	mlist=[]; 
	input_files = [os.path.join('sam', '%s.%s.bam' % (args.name, x)) for x in ['unique', 'control']]
	output_files = os.path.join('sam', '%s.filtered.bam' % args.name)
	script = get_script('filter_alignment.py', arguments={'--signal': input_files[0], '--control': input_files[1], '--output': output_files, '--features': " ".join(conf['filter_alignment']['features']), '--fdr': conf['filter_alignment']['fdr']}, package=chiflex_package)
	mlist.append(dependence(input_files, output_files, script));
	
	
	#sort filtered bam files for following clustering
	input_files = output_files
	output_files = os.path.join('sam', '%s.sorted.bam' % args.name)
	script = ['samtools', 'sort', '-o', output_files, '-T', 'temp', input_files]
	mlist.append(dependence(input_files, output_files, script));
	
	input_files = output_files
	output_files = os.path.join('sam', '%s.sorted.bam.bai' % args.name)
	script = ['samtools', 'index', input_files]
	mlist.append(dependence(input_files, output_files, script));
	
	
	#do preliminary clustering
	if(args.stranded):
		input_files = os.path.join('sam', '%s.sorted.bam' % args.name), os.path.join('sam', '%s.sorted.bam.bai' % args.name)
		output_files = os.path.join('clusters', 'minus.bedgraph')
		script = ['bedtools', 'genomecov', '-ibam', input_files[0], '-g', args.genome_index, '-bg', '-strand', '-', '>', output_files]
		mlist.append(dependence(input_files, output_files, script));
		
		input_files = os.path.join('clusters', 'minus.bedgraph'), os.path.join('sam', '%s.sorted.bam' % args.name)
		output_files = os.path.join('clusters', 'clusters.minus.gff')
		script = get_script('cluster_hits.py', arguments={'-s': '-', '--sbam': input_files[1]}, inp = input_files[0], out=output_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))


		input_files = os.path.join('sam', '%s.sorted.bam' % args.name), os.path.join('sam', '%s.sorted.bam.bai' % args.name)
		output_files = os.path.join('clusters', 'plus.bedgraph')
		script = ['bedtools', 'genomecov', '-ibam', input_files[0], '-g', args.genome_index, '-bg', '-strand', '+', '>', output_files]
		mlist.append(dependence(input_files, output_files, script));
		
		input_files = os.path.join('clusters', 'plus.bedgraph'), os.path.join('sam', '%s.sorted.bam' % args.name)
		output_files = os.path.join('clusters', 'clusters.plus.gff')
		script = get_script('cluster_hits.py', arguments={'-s': '+', '--sbam': input_files[1]}, inp = input_files[0], out=output_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))	


		input_files = os.path.join('clusters', 'clusters.plus.gff'), os.path.join('clusters', 'clusters.minus.gff')
		output_files = os.path.join('clusters', 'clusters.gff')
		script = ['cat', input_files[0], input_files[1], '>', output_files]
		mlist.append(dependence(input_files, output_files, script));

	else:
		input_files = os.path.join('sam', '%s.sorted.bam' % args.name), os.path.join('sam', '%s.sorted.bam.bai' % args.name)
		output_files = os.path.join('clusters', 'all.bedgraph')
		script = ['bedtools', 'genomecov', '-ibam', input_files[0], '-g', args.genome_index, '-bg', '>', output_files]
		mlist.append(dependence(input_files, output_files, script));
		
		input_files = os.path.join('clusters', 'all.bedgraph'), os.path.join('sam', '%s.sorted.bam' % args.name)
		output_files = os.path.join('clusters', 'clusters.gff')
		script = get_script('cluster_hits.py', arguments={'-s': '.', '--sbam': input_files[1]}, inp = input_files[0], out=output_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
		
		
	if(args.reassign):
		input_files = output_files
		output_files = os.path.join('clusters', 'clusters.assigned.gff') 
		script = get_script('assign_coordinates.py', arguments={}, inp = input_files, out = output_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
		
		
	if(args.annotation):
		input_files = output_files, args.annotation
		output_files = os.path.join('clusters', 'clusters.annotated.gff')
		script = get_script('annotate_bed.py', arguments={'--annotation': input_files[1]}, inp = input_files[0], out=output_files, package=chiflex_package)
		m.append(dependence(input_files, output_files, script))
		
	
	
	#makefile header
	mlist.insert(0, get_header(output_files))
	# makefie cleaner
	mlist.append('clean:\n\techo "nothing to clean."\n');
	return "\n\n".join(mlist)



#######################################################################################################################
#Function to create interaction Makefile		
def makefile_interaction():
	mlist=[];
	final_files = [];
	
	for itype in ['inter', 'intra']:
		input_files =  [os.path.join('chimeras', '%s.annotated.%s.gff' % (x, itype)) for x in ['unique', 'control']]
		output_files = os.path.join('interactions', 'filtered.%s.gff' % itype) 
		script = get_script('filter_chimera.py', arguments={'-s': input_files[0], '-c' : input_files[1], '--features': " ".join(conf['filter_chimera']['features']), '--fdr': conf['filter_chimera']['fdr']}, out = output_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
		
		input_files =  output_files
		output_files = os.path.join('interactions', 'sorted.%s.gff' % itype) 
		script = get_script('sort.py', inp=input_files, out = output_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script));
		
		input_files =  output_files
		output_files = os.path.join('interactions', 'interactions.%s.gff' % itype),  os.path.join('interactions', 'rid2iid.%s.bed' % itype)
		script = get_script('collapse2interaction.py', arguments={'-od': output_files[1], '--name': "%s_%s" % (args.name, itype), '--distance': conf['collapse2interaction']['distance']}, inp=input_files, out = output_files[0], package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script));
		output_files = os.path.join('interactions', 'interactions.%s.gff' % itype)
		
		input_files = output_files
		output_files = os.path.join('interactions', 'interactions.%s.itype.gff' % itype)
		script = get_script('spilt2subtypes.py', inp = input_files, out = output_files, package=interaction_package)
		mlist.append(dependence(input_files, output_files, script));
				
		if(args.exons):
			input_files = output_files, args.exons
			output_files = os.path.join('interactions', 'interactions.%s.itype.ktype.gff' % itype)
			if(args.stranded):
				cargs = arguments={'--exons': input_files[1], '--stranded': ''}
			else:
				cargs = arguments={'--exons': input_files[1]}
			script = get_script('annotate_chimera.py', arguments=cargs, inp = input_files[0], out=output_files, package=chiflex_package)
			mlist.append(dependence(input_files, output_files, script));
		
		if(args.annotation):
			input_files = output_files, args.annotation
			output_files = os.path.join('interactions', 'annotated.%s.gff' % itype)
			script = get_script('annotate_bed.py', arguments={'--annotation': input_files[1]}, inp = input_files[0], out=output_files, package=chiflex_package)
			mlist.append(dependence(input_files, output_files, script));
			
		final_files.append(output_files)
	
	#makefile header
	mlist.insert(0, get_header(final_files))
	# makefie cleaner
	mlist.append('clean:\n\techo "nothing to clean."\n');	
	return "\n\n".join(mlist)



#######################################################################################################################
#Function to create splicing Makefile
def makefile_splicing():
	mlist=[];
	final_files = [];
	
	for itype in ['lsj', 'csj']:
		input_files =  [os.path.join('chimeras', '%s.annotated.%s.gff' % (x, itype)) for x in ['unique', 'control']]
		output_files = os.path.join('splicing', 'filtered.%s.gff' % itype) 
		script = get_script('filter_chimera.py', arguments={'-s': input_files[0], '-c' : input_files[1], '--features': " ".join(conf['filter_chimera']['features']), '--fdr': conf['filter_chimera']['fdr']}, out = output_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script))
		
		if(args.exons):
			input_files = output_files, args.exons
			output_files = os.path.join('interactions', 'filtered.%s.ktype.gff' % itype)
			if(args.stranded):
				cargs = arguments={'--exons': input_files[1], '--stranded': ''}
			else:
				cargs = arguments={'--exons': input_files[1]}
			script = get_script('annotate_chimera.py', arguments=cargs, inp = input_files[0], out=output_files, package=chiflex_package)
			mlist.append(dependence(input_files, output_files, script));		
		
		input_files = output_files
		output_files = os.path.join('splicing', 'single.%s.gff' % itype)
		script = get_script('double2single.py', arguments={'--jtype': itype}, inp = input_files, out=output_files, package=splicing_package)
		mlist.append(dependence(input_files, output_files, script));
		
		input_files = output_files
		output_files = os.path.join('splicing', 'sorted.%s.gff' % itype)
		script = get_script('sort.py', inp=input_files, out = output_files, package=chiflex_package)
		mlist.append(dependence(input_files, output_files, script));		
			
		input_files = output_files
		output_files = os.path.join('splicing', 'collapsed.%s.gff' % itype), os.path.join('splicing', 'read2id.%s.bed' % itype)
		script = get_script('merge_junctions.py', arguments={'--jtype': itype, '--dictionary': output_files[1]}, inp = input_files, out=output_files[0], package=splicing_package)
		mlist.append(dependence(input_files, output_files, script));
		
		#output_files = os.path.join('splicing', 'collapsed.%s.gff' % itype)
		#final_files.append(os.path.join('splicing', 'read2id.%s.gff' % itype));
				
		final_files.extend(output_files);
		

	#makefile header
	mlist.insert(0, get_header(final_files))
	# makefie cleaner
	mlist.append('clean:\n\techo "nothing to clean."\n');	
	return "\n\n".join(mlist)



#######################################################################################################################
#Create folder structure
project_path = os.path.abspath(args.path)
folders = ['sam', 'reports', 'auxillary', 'log']
if('clustering' in args.modes):
	folders.append('clusters')
if('interaction' in args.modes):
	folders.append('interactions');
if('splicing' in args.modes):
	folders.append('splicing')
if('splicing' in args.modes or 'interaction' in args.modes):
	folders.append('chimeras');	


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
	
modes2functions = {'clustering': (makefile_clustering, 'makefile_clustering'), 'interaction': (makefile_interaction, 'makefile_interaction'), 'splicing': (makefile_splicing, 'makefile_splicing')}
for mode in args.modes:
	with open(os.path.join(project_path, modes2functions[mode][1]), 'w') as mf:
		mf.write(modes2functions[mode][0]());
	
	
def multipath(l):
	return "\t".join([os.path.abspath(x) for x in l])	
	
	
#report project call:
arguments_report = (
('name', ('Project name, assigned to the generated interactions', str)),
('modes', ('Modes of the project', str)),
('path', ('Project folder', os.path.abspath)),
('reads', ('Sequencing reads', multipath)),
('chiflex', ('Chiflex module used in the project', str)), 
('reference', ('Mapping reference index(genome, transcriptome) used for mapping', str)), 
('genome', ('Reference sequences(genome)', str)),
('exons', ('Exon system used for interactions annotation', str)), 
('annotation', ('Annotation system used for interactions annotation', str)), 
('repetitive', ('Repetitive mapped sequences are removed', str)),
('nonunique', ('Nonunique mappings with high alignment score are kept', str)),
('reports', ('html reports are produced', str)), 
('only_makefile', ('New Makefile was generated', str)),  
('reassign', ('Positions were reassigned to the genomic ones', str)),
('stranded', ('Sequencing data are implied to be stranded', str)), 
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



