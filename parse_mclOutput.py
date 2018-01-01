#!/usr/bin/env python
import sys, math, os, subprocess, datetime, argparse
from argparse import RawTextHelpFormatter


###################################################
### 0. script description and parsing arguments ###
###################################################
synopsis1 = "\
  parse mclOutput to a tab-delimited list of geneID and clusterID"
synopsis2 = "detailed description:\n\
 1. Input arguments and options:\n\
  - the input file, <mclOutput>, is an output of 'mcl', a tab-delimited text\n\
     with all members (geneIDs) in a cluster listed, one cluster per line.\n\
  - <clusterID_header> is a short string.  each cluster is named as\n\
     '<clusterID_header>_XXXXX',\n\
  - '-H'|'--Header': add a header line 'geneID\t<clusterID_header>'\n\
     ;default=False,\n\
  - '-r'|'--remove_speciesID': when geneID is formatted as 'spcsID|geneID',\n\
     remove 'spcsID|' part; default=False,\n\
 2. Output:\n\
  - write results to '${<mclOutput>%%.txt}.parsed.txt'\n\
  - clusters are numbered numerically, from '00001' to 'XXXXX'; number of '0's\n\
     is based on the number of entire clusters (i.e. rjust).\n\
by ohdongha@gmail.com 20171225 ver 0.2\n\n"

#version_history
#20171225 ver 0.2 modified to work with orthoMCL results to prepare CLfinder input
#20160511 ver 0.1 

parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)

# positional parameters
parser.add_argument('mclOutput', type=argparse.FileType('r'), help="See below")
parser.add_argument('clusterID_header', type=str, help="See below")

# options
parser.add_argument('-H', '--Header', action="store_true", default=False)
parser.add_argument('-r', '--remove_speciesID', action="store_true", default=False)

args = parser.parse_args()

if args.mclOutput.name.endswith(".txt"):
	output_filename = args.mclOutput.name[:-4] + ".parsed.txt"
else:	
	output_filename = args.mclOutput.name + ".parsed.txt"

fout = open(output_filename, 'w')

	
###############################
### 1. read input and parse ###
###############################

clusterID = 1
clusterID_dict = dict()
digit_for_clusterID = 0
geneID = ""

for line in args.mclOutput:
	tok = line.split()
	for i in range(0, len(tok)):
		if args.remove_speciesID == True:
			geneID = tok[i].strip().split('|')[1]
		else:
			geneID = tok[i].strip()
		clusterID_dict[geneID] = clusterID
	clusterID = clusterID + 1

digit_for_clusterID = int(math.log(clusterID,10)) + 1


#######################################
### 2. read input, parse, and write ###
#######################################
clusterID_header = args.clusterID_header + '_'

if args.Header == True:
	fout.write("geneID\t" + args.clusterID_header + '\n')

for key, value in sorted(clusterID_dict.iteritems(), key=lambda (k, v): (v,k), reverse=False):
	fout.write(key + '\t' + clusterID_header + str(clusterID_dict[key]).rjust(digit_for_clusterID, '0') + '\n')

print "done parsing %d clusters, for %s" % (clusterID, args.mclOutput.name)

args.mclOutput.close()
fout.close()