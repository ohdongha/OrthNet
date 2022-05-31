#!/usr/bin/env python
import sys, argparse
from argparse import RawTextHelpFormatter

#######################################################
### 0. script description, arguments, and functions ###
#######################################################
## 0.0 synopsis 
synopsis1 = "\
 - add <file2> to <file1>, by matching values in the 1st column of <file2>,\n\
    with the <N>st/nd/th column of the <file1>;\n\
 - similar to the old 'join_files_by_NthCol.py' but with more flexibility;\n"
  
synopsis2 = "detailed description:\n\
1. Input files, argument, and options\n\
 - <file1> and <file2> are tab-delimited text files;\n\
2. Merging options:\n\
 - '-m Mode' [1]: dealing with multiple matches in <file2>;\n\
    -m 1 : only the first occurence in <file2> is merged (default)\n\
    -m 2 : merge all matching <file2> lines, repeat the <file1> line\n\
 - '-n N' [1]: the <file1> column number that acts as keys to merge <file2>;\n\
    default is to use the 1st column of <file1>;\n\
 - '-H'|'--header' [False]: first lines of <file1> and <file2> are column headers\n\
    and merged regardless of matching;\n\
 - '-k'|'--keep_keys' [False]: include the 1st column of <file2> (keys) in the\n\
    merged file; default is not to include, since it will be redundant;\n\
*** options related to when a line in <file1> does not have any match in <file2>;\n\
 - '-x' [False]: do not print lines in <file1> without a match in <file2>;\n\
 - '-e filler_when_empty' ['na']: filler string to use when there is no match\n\
    in <file2>; -e '' will fill in with empty tabs;\n\
 - '-N totalCol' [0]: set the total number of columns in <file1> - lines in\n\
    <file1> with different column numbers will be printed as they are, without\n\
    looking for a match in <file2>; totalCol must be the same or larger than N;\n\
    if not given (default), all lines in <file1> will be matched to <file2>;\n\
3. Output:\n\
 - the resulting merged file is printed to STDOUT;\n\
 - the 1st column of <file2>, used as the key, is omitted in the merged file;\n\
by ohdongha@gmail.com 20220609 ver 0.3.2\n\n"
#version_history
# 20220609 ver 0.3.2 fix a bug where empty columns at the end of a records "stripped" + single column <file2>
# 20220531 ver 0.3.1 fix a bug where certain warnings are printed to stdout rather than stderr 
# 20210429 ver 0.3 add '-x' to remove <file1> lines without a match in <file2> 
# 20210124 ver 0.2 add '-N' to set the total column number and skip lines in <file1> with a different number of columns. 
# 20210102 ver 0.1 bug fix and add '-k' option
# 20201224 ver 0.0 modified from "join_files_by_NthCol.py"

## 0.1 parsing arguments 
parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)
# positional parameters
parser.add_argument('file1', type=argparse.FileType('r'))
parser.add_argument('file2', type=argparse.FileType('r'))
# options
parser.add_argument('-m', dest="mode", type=int, default= 1)
parser.add_argument('-n', dest="N", type=int, default= 1)
parser.add_argument('-H', '--header', action="store_true", default=False)
parser.add_argument('-k', '--keep_keys', action="store_true", default=False)
parser.add_argument('-x', dest="exclude_if_no_match", action="store_true", default=False)
parser.add_argument('-e', dest="filler_when_empty", type=str, default= "na")
parser.add_argument('-N', dest="totalCol", type=int, default= 0)
args = parser.parse_args()

# parse options
if args.N > 0:
	N = args.N
else:
	sys.stderr.write("Warning: '-n N' need a positive interger; set to default (n=1)")
	N = 1

if args.mode in {1,2}: 
	mode = args.mode
else:
	sys.stderr.write("Warning: '-m mode' parameter out of range; set to default (mode=1)")
	mode = 1

exclude_if_no_match = args.exclude_if_no_match
filler = args.filler_when_empty
number_of_fields_f1 = args.totalCol


##########################
### 1. reading <file2> ###
##########################
sys.stderr.write(" reading %s ...\n" % args.file2.name)

lines_in_f2_dict = dict() # key = string in the 1st column; value = list of lines
header_f2 = ""
number_of_fields_f2 = 0
header = True

for line in args.file2:
	line_to_add = ""
	tok = line.strip(" \n").split('\t') # v.0.3.2
	if args.keep_keys:
		line_to_add = line.strip(" \n")
	else:
		line_to_add = '\t'.join( tok[1:] )
	
	if header and args.header:
		header_f2 = line_to_add
	elif tok[0].strip() not in lines_in_f2_dict:
		lines_in_f2_dict[ tok[0].strip() ] = [ line_to_add ]
	else:
		lines_in_f2_dict[ tok[0].strip() ].append( line_to_add )

	if header: # counting fields in <file2> to "fill" with "filler" when there is no match,
		number_of_fields_in_f2 = len(tok)
		sys.stderr.write("Reading %d fields from %s \n" % (number_of_fields_in_f2, args.file2.name ) )
		header = False
	elif number_of_fields_in_f2 != len(tok):
		sys.stderr.write("Warning: a line in %s has a different number of fields:\n" % args.file2.name )
		sys.stderr.write(line)
		
args.file2.close()

if number_of_fields_in_f2 ==1: 
	if exclude_if_no_match:
		sys.stderr.write("%s has only one column and '-x' is on == extracting lines rather than merging ;)\n" % args.file2.name )
	else:
		sys.stderr.write("%s has only one column and '-x' is off - exiting since there is nothing to do ;)\n" % args.file2.name )
		sys.exit(0)	


#################################################
### 2. merging to <file1> and write to STDOUT ###
#################################################
sys.stderr.write(" merging to %s using column #%d ...\n" % (args.file1.name, N) )

header = True
num_line = 0

if args.keep_keys:
	number_of_fields_to_fill = number_of_fields_in_f2
else:
	number_of_fields_to_fill = number_of_fields_in_f2 - 1

for line in args.file1:
	line_merged = line.strip(" \n"); key_f1 = ""
	num_line += 1
	try:
		if args.header and header :
			if number_of_fields_in_f2 > 1:
				print( line_merged + '\t' + header_f2 ) 
			else:
				print( line_merged  )
			header = False
		elif number_of_fields_f1 == 0 or number_of_fields_f1 == len(line.split('\t')):
			key_f1 = line.split('\t')[N-1].strip()
			if key_f1 in lines_in_f2_dict:
				if number_of_fields_in_f2 == 1: # v.0.3.2
					print ( line.strip(" \n") )
				else:
					if mode == 1:
						print ( line.strip(" \n") + '\t' +  lines_in_f2_dict[key_f1][0] )
					else:
						for n in range(0, len( lines_in_f2_dict[key_f1] ) ):
							print ( line.strip(" \n") + '\t' +  lines_in_f2_dict[key_f1][n] )
			elif not exclude_if_no_match: #v.0.3
				for i in range(0, number_of_fields_to_fill):
					line_merged = line_merged + '\t' + filler
				print( line_merged )
		elif number_of_fields_f1 != len(line.split('\t')):
			print( line_merged ) # v0.2 skip <file1> lines not matching the set column number and print as they are, 
		
	except IndexError :
		sys.stderr.write( "Warning: line %d has a non-valid column in %s\n" % (num_line, args.file1.name) )

args.file1.close()