#!/usr/bin/env python

import re, os, sys, subprocess
import argparse
from argparse import RawTextHelpFormatter

###################################################
### 0. script description and parsing arguments ###
###################################################

synopsis1 = "\
  a script to create all possible pairs of blast commands given a list\n\
  of speciesIDs. This script creates only command lines and does not actually\n\
  run blast,"
synopsis2 = "\
 1. Input files and parameters:\n\
  - './<Project>.list' includes all species IDs <spcsIDs>, one per each line.\n\
  - '-q'|'--Path2Query': path to input .fa files in blast commands; default='.'\n\
  - '-d'|'--Path2DBfiles': path to blast blast DB files; default='.'\n\
  - '-o'|'--Path2Output': path to blast output files; default='.'\n\
  - '-n'|'--options': runs blastp instead of blastn; default=False\n\
  - '-p'|'--blastp': runs blastp instead of blastn; default=False\n\
  - '-I'|'--InputFormat': format string after <spcsID> in .fasta files for both\n\
     query and DB; default='.cds.rep.fa'\n\
 2. Output:\n\
  - print command lines for blastn (default) with tabular format to STDOUT,\n\
  - blast output files are formatted as 'out__<spcsID1>__vs__<spcsID2>.txt'\n\
 3. Example options (for blastn):\n\
  -n \"-task blastn -num_threads 16 -evalue 1e-5 -max_target_seqs 10 \n\
  -outfmt ' 6 std qlen slen stitle'\"\n\
 by ohdongha@gmail.com ver0.0 20171224\n"

#version_history
#20171224 ver 0.0

# an example blastn command: 
# blastn -task blastn -num_threads 16 -evalue 1e-5 -max_target_seqs 10 -outfmt ' 6 std qlen slen stitle' -query spcsID1.cds.fa -db spcsID2.cds.fa -out out_spcsID1__vs__spcsID2.txt

parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)

## positional arguments
parser.add_argument('Project', type=str, help="'./<Project>.list' includes spcsIDs being compared")

## options
parser.add_argument('-q', '--Path2Query', dest="Path2Query", type=str, default=".", help="see below")
parser.add_argument('-d', '--Path2DBfiles', dest="Path2DBfiles", type=str, default=".", help="see below")
parser.add_argument('-o', '--Path2Output', dest="Path2Output", type=str, default=".", help="see below") 
parser.add_argument('-n', dest="options", type=str, default="", help="see below") 
parser.add_argument('-p', '--blastp', action="store_true", default=False, help="see below")
parser.add_argument('-I', dest="InputFormat", type=str, default=".cds.rep.fa", help="see below") 

args = parser.parse_args()

# default formats
blastn_backbone = "blastn -query %s -db %s -out %s "
blastp_backbone = "blastp -query %s -db %s -out %s "
input_backbone = "%s" + args.InputFormat
output_backbone_blastn = "out__%s__vs__%s.bln.txt"
output_backbone_blastp = "out__%s__vs__%s.blp.txt"

# defining PATHs
path_query = args.Path2Query
if path_query[-1] != "/": path_query = path_query + "/"
path_db = args.Path2DBfiles
if path_db[-1] != "/": path_db = path_db + "/"
path_output = args.Path2Output
if path_output[-1] != "/": path_output = path_output + "/"

# variables
file_query = ""
file_db = ""
file_output = ""

spcsID1 = ""
spcsID2 = ""
command_string = ""


################################
### 1. reading the list file ###
################################
try:
	fin_SpcsList = open(args.Project + '.list')
except IOError:
	fin_SpcsList = open(args.Project)
	
spcsID_list = []

for line in fin_SpcsList:
	spcsID_list.append(line.strip())

fin_SpcsList.close()


########################################
### 2. printing blast command lines  ###
########################################
for spcsID1 in spcsID_list:
	for spcsID2 in spcsID_list:
		if spcsID1 != spcsID2:
		
			file_query = path_query + input_backbone % spcsID1
			file_db = path_db + input_backbone % spcsID2
			
			if args.blastp == True:
				file_output = path_output + output_backbone_blastp % (spcsID1, spcsID2)
				command_string = blastp_backbone % (file_query, file_db, file_output)
			else:
				file_output = path_output + output_backbone_blastn % (spcsID1, spcsID2)
				command_string = blastn_backbone % (file_query, file_db, file_output)

			command_string = command_string + args.options
			print command_string