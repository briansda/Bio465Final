#!/usr/bin/python

import re
import sys
print("Fixing " + sys.argv[1])
with open( sys.argv[ 1 ], "r" ) as fh:
	lines = fh.readlines()
	header = lines[0]
	#print header
	header = header.split('\t')
	#print header
	p = re.compile(".*(GSM\d{4,8}).*")
	for i in range(1, len(header)):
		m = p.match(header[i])
		header[i] =  m.group(1)
	lines[0] = "\t".join(header)
	lines[0] = lines[0]+'\n'
	fh.close()
	fh = open(sys.argv[1], 'w')	
	fh.writelines(lines)

	

fh.close()
sys.exit()
	
