#!/usr/bin/env python
import sys, os, argparse
from argparse import RawTextHelpFormatter


###################################################
### 0. script description and parsing arguments ###
###################################################
synopsis1 = "\
  select representative gene models from a .gtfParsed.txt file"
synopsis2 = "detailed description:\n\
 - from <gtfParsed_input>, select gene models in <representative_list> and\n\
    print to <gtfParsed_rep_output>.\n\
 - unless otherwise indicated with an option, print the first geneID that\n\
    includes an item in <represenative_list>, i.e. no exact match needed,\n\
 - '-x'|'--exact_match': print only gene models with the first exact match;\n\
    default=False,\n\
 - '-e'|'--ends_with': print only gene models with the first match at the end\n\
    of geneID; default=False, ignored if '-x' is on,\n\
 - '-r'|'--replace': replace the geneID with the one in <representative_list>,\n\
    default=False, ignored if '-x' is on.\n\
by ohdongha@gmail.com 20171223 ver 0.0\n\n"

#version_history
#20171223 ver 0.0 

parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)

# positional parameters
parser.add_argument('representative_list', type=argparse.FileType('r'))
parser.add_argument('gtfParsed_input', type=argparse.FileType('r'))
parser.add_argument('gtfParsed_rep_output', type=argparse.FileType('w'))

# options
parser.add_argument('-x', '--exact_match', action="store_true", default=False)
parser.add_argument('-e', '--ends_with', action="store_true", default=False)
parser.add_argument('-r', '--replace', action="store_true", default=False)

args = parser.parse_args()


###########################################
### 1. reading in <representative_list> ###
###########################################
representatives = set()

print "reading %s as the <representative_list>:" % args.representative_list.name

for line in args.representative_list:
	tok = line.split()
	if tok[0].strip() in representatives:
		print "duplicated item: %s" % tok[0].strip()
	else:
		representatives.add( tok[0].strip() )

print "%d gene models accepted from %s," % ( len(representatives), args.representative_list.name)
args.representative_list.close()


########################################################################
### 2. process <gtfParsed_input> and write to <gtfParsed_rep_output> ###
########################################################################
geneID =""
num_geneID_found = 0
num_line = 0
new_line = ""
header = True

for line in args.gtfParsed_input:
	tok = line.split('\t')
	geneID = tok[0].strip()
	num_line += 1
	if header == True:
		args.gtfParsed_rep_output.write(line)
		header = False
	else:
		if args.exact_match == True:
			if geneID in representatives:
				num_geneID_found += 1
				representatives.remove(geneID)
				args.gtfParsed_rep_output.write(line)
		else:
			for item in representatives:
				if (args.ends_with == False and item in geneID) or geneID.endswith(item):
					num_geneID_found += 1
					representatives.remove(item)
					if args.replace == True:
						tok[0] = item
						new_line = '\t'.join(tok)
						args.gtfParsed_rep_output.write(new_line)						
					else:
						args.gtfParsed_rep_output.write(line)
					break
	
	# counter display
	if ( num_geneID_found % 500 == 0):
		sys.stdout.write("\r   found %d representatives among %d gene models." % ( num_geneID_found, num_line ) )
		sys.stdout.flush()

print "\ntotal %d representatives found and printed to %s." % (num_geneID_found, args.gtfParsed_rep_output.name)
if len(representatives) > 0:
	print "the following representatives couldn't be found:"
	print representatives

args.gtfParsed_input.close()
args.gtfParsed_rep_output.close()



