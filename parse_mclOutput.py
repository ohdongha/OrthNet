#!/usr/bin/env python

synopsis = "\n\n usage: parse_mclOutput.py <mclOutput> <clusterID_header, default = 'c'>\n\
 \n\
 1. <mclOutput>\n\
	- an output of 'mcl', a tab-delimited text file with all members listed, one cluster per line.\n\
 2. Output\n\
	- '<mclOutput>.parsed.txt' will be created.\n\
#	- a header 'memberID\tclusterID' in the first line.\n\
	- afterwards, <memberID> and <clusterID>, tab-delimited, one member per line.\n\
	- for <clusterID>, clusters are numbered numerically, from 'c00001' to 'cXXXXX'.\n\
	- number of '0's will be adjusted based on the number of entire clusters.\n\
	- with <clusterID_header> specified, <clusterID> will be formatted as '<clusterID_header>_XXXXX'.\n\
\n\
 by ohdongha@gmail.com \n\
 20160511 ver 0.1\n"

import sys
import math
	
try: 
	fin = open(sys.argv[1], "rU")
	fout = open(sys.argv[1].strip() + ".parsed.txt", "w")
	if len(sys.argv) > 1:
		clusterID_header = sys.argv[2]
	else:
		clusterID_header = 'default'
except (ValueError, IndexError) :
	print synopsis
	sys.exit(0)

#assign OG group number for genes in <mclOutput> and write to <output.txt>
#from contextlib import nested

clusterID = 1
clusterID_dict = dict()
digit_for_clusterID = 0

for line in fin:
	tok = line.split()
	for i in range(0, len(tok)):
		clusterID_dict[tok[i].strip()] = clusterID
	clusterID = clusterID + 1
fin.close()

digit_for_clusterID = int(math.log(clusterID,10)) + 1

#fout.write("memberID\tclusterID\n")
for key, value in sorted(clusterID_dict.iteritems(), key=lambda (k, v): (v,k), reverse=False):
	if clusterID_header == 'default':
		fout.write(key + '\tc' + str(clusterID_dict[key]).rjust(digit_for_clusterID, '0') + '\n')
	else:
		fout.write(key + '\t' + clusterID_header + '_' + str(clusterID_dict[key]).rjust(digit_for_clusterID, '0') + '\n')
fout.close()

print "done parsing %d clusters, for %s" % (clusterID, sys.argv[1])
