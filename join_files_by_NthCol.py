#!/usr/bin/env python

synopsis = "\njoin_files_by_NthCol.py <file1.txt> <N> <header_switch> <file2.txt> <output.txt>\n\
# <file1.txt> and <file2.txt> are tab-delimited.\n\
# add <file2.txt> to <file1.txt>, by matching values in the 1st column of <file2.txt> ...\n\
# ... with <N>th column of the <file1.txt>.\n\
# if <header_switch> = 1, assumes there are headers in both files. <header_switch> = 0, when files have no header.\n\
# all lines in <file1.txt> will be kept, with fields with no matches filled with 'na'.\n\
# if multiple lines with the same value in the 1st column exist in <File2.txt>, only the first occurance will be used.\n\
# the 1st column of <File2.txt> will be omitted. \n\
# 20170715 ver 1.1 reporting more details when a line has an error,\n\
# 20150904 ver 1.0\n\
# copyleft by ohdongha@gmail.com\n\n"

import sys

try: 
	fin_file1 = open(sys.argv[1], "rU")
	col_to_match = int(sys.argv[2]) - 1
	header_switch = int(sys.argv[3])
	fin_file2 = open(sys.argv[4], "rU")
	fout = open(sys.argv[5], "w")
except ValueError :
	print synopsis
	sys.exit(0)
except IndexError :
	print synopsis
	sys.exit(0)

line_in_file2_dict = {}
header_file2 = ""
header_or_not = 1
number_of_fields_in_file2 = 0

## reading the file2
print "\n### reading " + sys.argv[4]
for line in fin_file2:
	tok = line.split('\t')
	number_of_fields_in_file2 = len(tok)
	if header_switch == 1 and header_or_not == 1 :
		for i in range(1, number_of_fields_in_file2):
			header_file2 = header_file2 + '\t' + tok[i].strip()
		header_or_not = 0	
	elif tok[0].strip() not in line_in_file2_dict:
		line_in_file2_dict[tok[0].strip()] = ""
		for i in range(1, number_of_fields_in_file2):
			line_in_file2_dict[tok[0].strip()] = line_in_file2_dict[tok[0].strip()] + '\t' + tok[i].strip()
fin_file2.close()

## adding file2 fileds to the file1
print "\n### adding to " + sys.argv[1] + " and writing to " + sys.argv[5]

file1_ID = ""
new_line = ""
header_or_not = 1
num_line = 0

for line in fin_file1:
	num_line += 1
	try:
		new_line = line.strip('\n')
		if header_switch == 1 and header_or_not == 1 :
			new_line = new_line + header_file2 
			header_or_not = 0
		else:
			file1_ID = new_line.split('\t')[col_to_match].strip()
			if file1_ID not in line_in_file2_dict:
				for i in range(1, number_of_fields_in_file2):
					new_line = new_line + '\t' + 'na'
			else:
				new_line = new_line + line_in_file2_dict[file1_ID]
		fout.write(new_line + '\n')
	except ValueError:
		print "Line %d has a non-valid value in %s," % (num_line, sys.argv[1])
	except IndexError :
		print "Line %d has a non-valid column in %s," % (num_line, sys.argv[1])
		
print "\n### done\n "
fin_file1.close()
fout.close()

