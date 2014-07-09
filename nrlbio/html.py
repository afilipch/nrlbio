def get_link(interval, system, internal = False):	
	if(internal):
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
		
	return "%s%s%s:%d..%d" % (source, system_string, interval.chrom , interval.start+1, interval.stop)	