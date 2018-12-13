#!/usr/bin/env python
import os, sys, subprocess, datetime
import argparse
from argparse import RawTextHelpFormatter

###################################################
### 0. script description and parsing arguments ###
###################################################

synopsis1 = "\
  create a tabulated summary of a 'CL_finder_muli.py' run on multiple genomes.\n"

synopsis2 = "detailed description:\n\
 0. Pre-requisite: a 'CL_finder_multi.py' run on multiple genomes,\n\
 1. Input files and parameters:\n\
  - './<Project>.list' includes all species IDs (spcsIDs), one per each line.\n\
  - '-p'|'--Path2CLfm_output': path to the 'CL_finder_multi.py' output folder\n\
     that contains 'spcsID.CL_compared2*.rep1.txt' files for all species,\n\
     created by a 'CL_finder_multi.py -r' run; defaults = './'\n\
 2. Output and options\n\
  - for all query-target genome pairs, print the numbers of cl, ls, nd, tr\n\
     best-hit pairs in a tabulated format and print to <Output>.\n\
  - '-e'|'--cl_end': report 'cl_end' best-hit pairs separately from 'cl' pairs\n\
     ;default=False,\n\
 by ohdongha@gmail.com ver0.0.1 20181204\n"
 
#version_history
#20181204 ver 0.0.1 remove temporary file after running
#20180101 ver 0.0

parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)

## positional arguments
parser.add_argument('Project', type=str, help="'./<Project>.list' includes spcsIDs being compared")
parser.add_argument('Output', type=argparse.FileType('w'))

## options to receive PATH to input files
parser.add_argument('-p', '--Path2CLfm_output', dest="Path2CLfm_output", type=str, default=".", help="see below")

## other options
parser.add_argument('-e', '--cl_end', action="store_true", default=False, help="see below")

args = parser.parse_args()

# defining PATH and the temporary file name
path_CLfm_output = args.Path2CLfm_output
if path_CLfm_output[-1] != "/": path_CLfm_output = path_CLfm_output + "/"

tempFile_name = path_CLfm_output.split('/')[-2] + '_' + datetime.datetime.now().strftime("%y%m%d_%H%M%S") + '.rep1Combined.temp'


################################
### 1. reading the list file ###
################################
projectID = args.Project
try:
	fin_SpcsList = open(projectID + '.list', 'r')
except IOError:
	fin_SpcsList = open(projectID, 'r')
	projectID = projectID[:-5]
spcsID_list = []

print "\nreading the list file:" + fin_SpcsList.name
for line in fin_SpcsList:
	spcsID_list.append(line.strip())
	print line.strip()
fin_SpcsList.close()
print "Total %d species IDs detected." % (len(spcsID_list))


##############################################
### 2. create the temp file and reading it ###
##############################################
# create dictionaries: num_BestHitPairs[query_spcsID][target_spcsID][CL_type] = number_of_BestHitPairs
num_BestHitPairs = dict()

for query_spcsID in spcsID_list:
	num_BestHitPairs[query_spcsID] = dict()
	for target_spcsID in spcsID_list:
		if query_spcsID != target_spcsID:
			num_BestHitPairs[query_spcsID][target_spcsID] = dict() # dictionary with key = CL_type, value = number_of_BestHitPairs

# create and open the Rep1Combined file
subprocess.call("cat " + path_CLfm_output + "*.rep1.txt > ./" + tempFile_name, shell=True)
fin_Rep1Combined = open("./" + tempFile_name, 'r')

CL_type_set = set(['cl', 'ls', 'nd', 'tr'])
# take care of '-e' option here:
if args.cl_end == True:
	CL_type_set.add('cl_end')

# start parsing the Rep1Combined file
query_spcsID = ""
target_spcsID = ""
target_list = list()
CL_type = ""

for line in fin_Rep1Combined:
	tok = line.split('\t')
	if tok[0].strip().endswith('-'):
		if tok[0].strip()[:-1] in spcsID_list:
			query_spcsID = tok[0].strip()[:-1]
			print "processing the query species %s" % query_spcsID
		else:
			print "something is wrong, %s is not in the list of species IDs!" %  tok[0].strip()[:-1]
			print line
		target_list = tok[1:]
	elif tok[0].strip().split('_')[0] in CL_type_set:
		# take care of '-e' option here:
		if args.cl_end == True and tok[0].strip() == 'cl_end':
			CL_type = 'cl_end'
		else:
			CL_type = tok[0].strip().split('_')[0]
		# reading in num_BestHitPairs values	
		for i in range(0, len(target_list)):
			target_spcsID = target_list[i].split('_')[0]
			try:
				num_BestHitPairs[query_spcsID][target_spcsID][CL_type] = num_BestHitPairs[query_spcsID][target_spcsID].get(CL_type, 0) + int(tok[i+1])
			except ValueError:
				print "something is wrong, a non-numerical value detected among numbers of best-Hit pairs ..."
				print line
	else:
		print "something is wrong, there is a line that is neither a header nor a list of num_BestHitPairs ..."
		print line

fin_Rep1Combined.close()
print "completed reading .rep1 files.  now printing the summary to %s" % args.Output.name

#######################
### 3. print Output ###
#######################
# print header
args.Output.write("query_genome\tCL_type\t" + '\t'.join(spcsID_list) + '\n')
row_to_print = ""

# start printing the matrix
for query_spcsID in spcsID_list:
	for CL_type in sorted(list(CL_type_set)):
		row_to_print = query_spcsID + '\t' + CL_type
		for target_spcsID in spcsID_list:
			if query_spcsID == target_spcsID:
				row_to_print += '\t' # empty cells when query and target are the same
			else:
				row_to_print += ('\t' + str( num_BestHitPairs[query_spcsID][target_spcsID].get(CL_type, 0) ) )
				
		args.Output.write(row_to_print + '\n')

args.Output.close()	

## clean the .temp file
try:
    os.remove(tempFile_name)
except OSError:
    pass
	
print "\ndone\n"
