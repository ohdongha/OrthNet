#!/usr/bin/env python
import re, os, sys, subprocess, argparse
from argparse import RawTextHelpFormatter

###################################################
### 0. script description and parsing arguments ###
###################################################

synopsis1 = "\
create all possible pairs of blast (or mmseqs2) commands given a list of \n\
speciesIDs. This script creates only command lines and does not actually \n\
run blast (or mmseqs2);"
synopsis2 = "\
1. Input files and parameters:\n\
 - './<Project>.list' includes all species IDs <spcsIDs>, one per each line.\n\
 - '-q'|'--Path2Query': path to input .fa files in blast commands; default='.'\n\
 - '-d'|'--Path2DBfiles': path to blast blast DB files; default='.'\n\
 - '-o'|'--Path2Output': path to blast output files; default='.'\n\
 - '-n'|'--options': options to add to base command lines; default=''\n\
 - '-P'|'--blastp': create blastp commands instead of blastn; default=False\n\
 - '-M'|'--MMSeqs2': create MMSeqs2 commands, instead of blast; default=False\n\
 - '-I'|'--InputFormat': format string after <spcsID> in .fasta files for both\n\
    query and DB; ignored if '-M' is on; default='.cds.rep.fa'\n\
2. Output:\n\
 - print command lines for blastn (default) with tabular format to STDOUT,\n\
 - blast outputs are formatted as 'out__<spcsID1>_vs_<spcsID2>.[bln,blp].txt'\n\
3. Misc:\n\
 - assume blast DBs (or mmseqs DBs and indexes) are already created;\n\
 - suggested option for blastn: -n \"-task blastn -num_threads 16 -evalue 1e-5\n\
    -max_target_seqs 10 -outfmt ' 6 std qlen slen stitle'\"\n\
 - suggested option for mmseqs2: -n \"--max-seqs 10\"\n\
 - with '-M' option, assume there is a folder '.\tmp' for temporary files;\n\
 - with '-M' option, print also 'mmseqs convertalis' commands in the format of:\n\
    mmseqs convertalis <spcsID1>_DB <spcsID2>_DB out__<spcsID1>_vs_<spcsID2> \\\n\
    out__<spcs1>_vs_<spcs2>.mmseqs2.txt --format-mode 2\n\n\
by ohdongha@gmail.com ver0.1 20180507\n"

#version_history
#20180507 ver 0.1 # added option to print mmseqs2 commands, insead of blastn
#20171224 ver 0.0

## example blastn command: 
# blastn -task blastn -num_threads 16 -evalue 1e-5 -max_target_seqs 10 -outfmt ' 6 std qlen slen stitle' -query spcsID1.cds.fa -db spcsID2.cds.fa -out out_spcsID1__vs__spcsID2.txt
## example MMSeqs2 command:
# mmseqs search spcsID1_DB spcsID2_DB out__spcsID1_vs_spcsID2 tmp --max-seqs 10 
# mmseqs convertalis spcsID1_DB spcsID2_DB out__spcsID1_vs_spcsID2 out__spcsID1_vs_spcsID2.mmseqs2.txt --format-mode 2

parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)

## positional arguments
parser.add_argument('Project', type=str, help="'./<Project>.list' includes spcsIDs being compared")

## options
parser.add_argument('-q', '--Path2Query', dest="Path2Query", type=str, default=".", help="see below")
parser.add_argument('-d', '--Path2DBfiles', dest="Path2DBfiles", type=str, default=".", help="see below")
parser.add_argument('-o', '--Path2Output', dest="Path2Output", type=str, default=".", help="see below") 
parser.add_argument('-n', dest="options", type=str, default="", help="see below") 
parser.add_argument('-P', '--blastp', action="store_true", default=False, help="see below")
parser.add_argument('-M', '--MMSeqs2', action="store_true", default=False, help="see below")
parser.add_argument('-I', dest="InputFormat", type=str, default=".cds.rep.fa", help="see below") 

args = parser.parse_args()

# default formats
cmd_fmt_blastn = "blastn -query %s -db %s -out %s "
cmd_fmt_blastp = "blastp -query %s -db %s -out %s "
input_fmt = "%s" + args.InputFormat
out_fmt_blastn = "out__%s_vs_%s.bln.txt"
out_fmt_blastp = "out__%s_vs_%s.blp.txt"

cmd_fmt_mmseqs_search = "mmseqs search %s %s %s ./tmp "
cmd_fmt_mmseqs_convertalis = "mmseqs convertalis %s %s %s %s --format-mode 2 "
out_fmt_mmseqs_search = "out__%s_vs_%s "
out_fmt_mmseqs_convertalis = "out__%s_vs_%s.mmseqs2.txt "

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
file_output_m = ""

spcsID1 = ""
spcsID2 = ""
cmd_str = ""


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
		
			if args.MMSeqs2:
				file_query = path_query + "%s_DB" % spcsID1
				file_db = path_db + "%s_DB" % spcsID2				
			else:
				file_query = path_query + input_fmt % spcsID1
				file_db = path_db + input_fmt % spcsID2
			
			if args.MMSeqs2:
				file_output = path_output + out_fmt_mmseqs_search % (spcsID1, spcsID2)
				file_output_m = path_output + out_fmt_mmseqs_convertalis % (spcsID1, spcsID2)
				cmd_str = cmd_fmt_mmseqs_search % (file_query, file_db, file_output) + args.options
				cmd_str = cmd_str + '\n' + cmd_fmt_mmseqs_convertalis % (file_query, file_db, file_output, file_output_m)
			elif args.blastp:
				file_output = path_output + out_fmt_blastp % (spcsID1, spcsID2)
				cmd_str = cmd_fmt_blastp % (file_query, file_db, file_output) + args.options
			else:
				file_output = path_output + out_fmt_blastn % (spcsID1, spcsID2)
				cmd_str = cmd_fmt_blastn % (file_query, file_db, file_output) + args.options

			print cmd_str
