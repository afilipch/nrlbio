'''removes all gaps in provided file'''
import sys;

nucl = set("ACTGactg")

with open(sys.argv[1]) as f:
	for l in f:
		t = [];
		if(l.startswith('>')):
			print l.strip()
		else:
			for s in l.strip("\n"):
				if s in nucl:
					t.append(s)
			print "".join(t);
				