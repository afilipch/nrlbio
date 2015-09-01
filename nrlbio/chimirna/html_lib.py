'''library to produce html reports'''
import gconst;
import copy;
import sys;

class HTML_read(object):
	def __init__(self):
		self.left_mapped = "";
		self.right_mapped = "";
		self.unmapped = ["","",""]

class HTML_part(object):
	def __init__(self, padding, lines, conversion_lines):
		self.padding = padding;
		self.lines = copy.copy(lines);
		if(conversion_lines):
			self.conversion_lines = conversion_lines
			if(len(conversion_lines) > 1 and "".join(conversion_lines[-2:]) == "-"):
				#sys.stderr.write("%s\n" % conversion_lines);
				conversion_lines[-2] = "";			
		else:
			pass;
			#self.conversion_lines = [""]
			#sys.stderr.write("%s\n" % "BF");
		
		
class HTML_interaction(object):
	def __init__(self, interaction, lparts, rparts, seqs, rids, system, refseq, depth):
		self.rids = rids;
		self.interaction = interaction
		#sys.stderr.write(str(self.interaction) + "\n")
		self.refseq = refseq
		self.depth = depth		
		self.link = self.get_link(system, inner=False)
		self.left_parts = self.get_lparts(lparts);
		self.right_parts = self.get_rparts(rparts);
		self.reads = self.get_reads(seqs, lparts, rparts)
		

		#sys.stderr.write(self.link + "\n")		
		#self.test();
		
	def get_link(self, system, inner = True):	
		if(inner):
			source = "http://141.80.186.52/"
		else:
			source = "http://genome.ucsc.edu/"	
	
		if(system == "ce6"):
			system_string = "cgi-bin/hgTracks?hgsid=7901&org=C.elegans&db=ce6&position=";
		elif(system == "hg19"):
			system_string = "cgi-bin/hgTracks?hgsid=7901&org=H.sapiens&db=hg19&position=";
		elif(system == "mm9"):
			system_string = "cgi-bin/hgTracks?hgsid=7901&org=M.musculus&db=mm9&position="; 
		else:
			pass;
			
		return "%s%s%s:%d..%d" % (source, system_string, self.interaction.chromosome , self.interaction.start, self.interaction.end)	
		
	def get_lparts(self, lparts):
		#if( self.interaction.iid == 'ago2_kishore-pc01495'):
			#sys.stderr.write("%d\t%s\n" % (lparts[0].read_end, self.rids[0]));		
		ans = [];
		for chipart in lparts:
			padding  = chipart.ref_start
			lines = [];
			conversion_lines = [];
			
			aligned = list(self.interaction.mirseq[chipart.ref_start: chipart.ref_end]);
			conversions = [];
			for c in chipart.conversions:
				conversions.append(c.split(","))
				conversions[-1][0] = int(conversions[-1][0]) - chipart.ref_start
			conversions.sort(key= lambda x: x[0])	
			start = 0;	
			for c in conversions:
				if(chipart.ref_start <= c[0] < chipart.ref_end):
					lines.append("".join(aligned[start:c[0]]));
					conversion_lines.append(c[2]);
					start = c[0] + 1;
				else:
					pass;
			lines.append("".join(aligned[start:]));
			conversion_lines.append("");
			if(not lines):
				lines.append("".join(aligned));
			ans.append(HTML_part(padding, lines, conversion_lines))
		return ans
		
	def get_rparts(self, rparts):
		ans = []
		for rp in rparts:
			if(self.interaction.strand == "+"):
				pad = rp.genome_start - self.interaction.start;
			else:
				pad = self.interaction.end - rp.genome_end ;

				
			padding = self.depth+pad;
			lines = [];
			conversion_lines = [];
			line = ""
			conv_line = "";
			refseq = self.refseq[padding:]	
			
			#if(len(rp.aligned) > len(refseq)):
				#ans.append(HTML_part(padding, [rp.aligned], [""]))		
				#continue;
			#sys.stderr.write("%s\t%s\t%d\n" % (self.refseq, rp.aligned, padding))
			#sys.stderr.write(str(self.interaction) + "\n")
			#sys.stderr.write("%s\n%s\n%s\n%d\t%d\n" % (self.refseq, refseq, rp.aligned, rp.genome_start, self.interaction.start))
			for i in range(len(rp.aligned)):
				if(rp.aligned[i] == refseq[i]):
					line += rp.aligned[i];
				else:
					conv_line = rp.aligned[i];
					lines.append(line);
					conversion_lines.append(conv_line)
					conv_line = ""
					line = ""
			lines.append(line);
			conversion_lines.append(conv_line)					
			ans.append(HTML_part(padding, lines, conversion_lines))		
		return ans;	
	
	def get_reads(self, seqs, lparts, rparts):
		#if( self.interaction.iid == 'ago2_kishore-pc01495-2'):
			#sys.stderr.write("%d\t%s\t%s\n" % (lparts[0].read_end, self.rids[0], seqs[0]));
		ans = [];
		
		#for i in range(len(seqs)):
			#if(self.rids[i] == '@ago2_kishore-pc_9202443_x1'):
				#html_read = HTML_read();
				#html_read.unmapped[0] = seqs[i][:lparts[i].read_start];
				#html_read.left_mapped = seqs[i][lparts[i].read_start:lparts[i].read_end]
				#rseq = rparts[i].aligned.replace("-","");
				##sys.stderr.write('%s\t%s\n' % (rseq, rparts[i].aligned));
				##conv = rparts[i].conversion - rparts[i].read_start;
				##if(conv >= 0 and conv < len(rseq)):
					##rseq[conv] = "C";
				##rseq = "".join(rseq);
				#sys.stderr.write('%s\t%s\n' % (rseq, rparts[i].aligned));
				#rstart = seqs[i].rfind(rseq);
				#html_read.unmapped[1] = seqs[i][lparts[i].read_end:rstart];
				#html_read.right_mapped = rseq;
				#html_read.unmapped[2] = seqs[i][len(rseq)+rstart:];	
		
		for i in range(len(seqs)):
			html_read = HTML_read();
			html_read.unmapped[0] = seqs[i][:lparts[i].read_start];
			html_read.left_mapped = seqs[i][lparts[i].read_start:lparts[i].read_end]
			rseq = list(rparts[i].aligned.replace("-",""));
			#conv = rparts[i].conversion - rparts[i].read_start;
			#if(conv >= 0 and conv < len(rseq)):
				#rseq[conv] = "C";
			rseq = "".join(rseq);
			rstart = seqs[i].rfind(rseq);
			html_read.unmapped[1] = seqs[i][lparts[i].read_end:rstart];
			html_read.right_mapped = rseq;
			html_read.unmapped[2] = seqs[i][len(rseq)+rstart:];
			ans.append(html_read);
		return ans;	
			
			
			
			
			
				
	def test(self):
		#print self.interaction.mirseq
		#print
		#for hr in self.left_parts:
			#line = hr.padding;
			#for i in range(len(hr.lines)):
				#line += hr.lines[i] + hr.conversion_lines[i].lower();
			#print line
		#print "\n"*2	
		
		print self.interaction.tseq, "target"
		print
		for hr in self.right_parts:
			line = hr.padding;
			for i in range(len(hr.lines)):
				line += hr.lines[i] + hr.conversion_lines[i].lower();
			print line
		print "\n"*2
		
		
		
			
			


def interaction2link(interaction, system):
	if(system == "ce6"):
		system_string = "<a href=\"http://141.80.186.52/cgi-bin/hgTracks?hgsid=7901&org=C.elegans&db=ce6&position=";
	elif(system == "hg19"):
		system_string = "<a href=\"http://141.80.186.52/cgi-bin/hgTracks?hgsid=7901&org=H.sapiens&db=hg19&position=";
	elif(system == "mm9"):
		system_string = "<a href=\"http://141.80.186.52/cgi-bin/hgTracks?hgsid=7901&org=M.musculus&db=mm9&position="; 
	link = "%s%s:%d..%d\" target=\"_blank\">ucsc link</a>" % (system_string, interaction.chromosome , interaction.start, interaction.end)	
	print r"<pre>"
	print "\t".join([str(x) for x in [interaction.chromosome, interaction.start - 1, interaction.end, interaction.strand, interaction.mirid, interaction.score, link]])
	#print r"</pre>" 
	#print r"<pre>"	
	#print
	#print r"</pre>"	
	
def full_representation(interaction, mydict, system):
	link = "%s%s:%d..%d\" target=\"_blank\">ucsc link</a>" % (gconst.links[system], interaction.chromosome , interaction.start, interaction.end)
	#mirna section
	print interaction.mirseq;
	print
	for rid, d in mydict.items():
		chipart = d[0]
		line  = " "*chipart.ref_start
		aligned = list(interaction.mirseq[chipart.ref_start: chipart.ref_end]);
		conversions = [];
		for c in chipart.conversions:
			conversions.append(c.split(","))
			conversions[-1][0] = int(conversions[-1][0]) - chipart.ref_start
		for c in conversions:
			if(chipart.ref_start <= c[0] < chipart.ref_end):
				aligned[c[0]] = c[2].lower();
		line = line + "".join(aligned) + " "*50
		print line	
			
	print "\n"*4		
    


