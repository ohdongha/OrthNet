#!/usr/bin/env python
import sys

synopsis = "\n\
join_files_by_NthCol.py <file1> <N> <header_switch> <file2> <output.txt> [e]\n\
 - <file1> and <file2> are tab-delimited text files,\n\
 - add <file2> to <file1>, by matching values in the 1st column of <file2>,\n\
     with the <N>st/nd/th column of the <file1>,\n\
 - <header_switch> == 1, if both files have headers; 0, if both have no header,\n\
 - all lines in <file1> are kept, with fields with no matches filled with 'na',\n\
 - if [e] is specified, fields with no matches are [e]mpty, instead of 'na',\n\
 - uses only the first occurrence of the same value in the 1st column in <file2>\n\n\
# copyleft by ohdongha@gmail.com 20180317 ver 1.1.1\n\n"

# 20180317 ver 1.1.1 options to add empty cells, instead of 'na,' for non-matches,\n\
# 20170715 ver 1.1 reporting more details when a line has an error,\n\
# 20150904 ver 1.0\n\

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

empty_nonMatch = False
if len(sys.argv) > 6:
	if sys.argv[6] == 'e':
		empty_nonMatch = True
	
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
					if empty_nonMatch:
						new_line = new_line + '\t'					
					else:
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