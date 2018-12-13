#!/usr/bin/env python

import re, sys, os
import argparse
from argparse import RawTextHelpFormatter

###################################################
### 0. script description and parsing arguments ###
###################################################

synopsis1 = "\
  - updates 'BestHitPairs' files so that reciprocal 'cl' relationship of\n\
     'non-bestHit' is preferred than non-reciprocal relationship of 'bestHit'.\n\
  - accepts '.4OrthNet.input' file generated with 'CL_finder_multi.py -n'.\n\
  - updated files can be fed back to 'CL_finder_multi.py -u'.\n"

synopsis2 = "detailed description:\n\
 1. Input files and parameters\n\
  - './<Project>.list' includes all species IDs (spcsIDs), one per line.\n\
  - <4OrthNet_input> is an output file of 'CL_finder_multi.py -n'; typical\n\
     filename == '<Project>.4OrthNet.input'; see 'CL_finder_mulyi.py -h'.\n\
  - '-b'|'--Path2BHPs': path to 'BestHitPairs(BHPs)' files for all possible\n\
     species pairs with file names 'BestHits__spcsID1__vs__spcsID2.list'.\n\
  - 'BHPs' files can be tabulated blast+ results with the '-max_target_seqs N'\n\
     (N>=10) option; i.e. multiple subjects for each query, with the most \n\
     similar 'BestHit', followed by next best subjects.\n\
 2. Updating 'BestHitPairs' file\n\
  - if CL_type(q->s)!='cl' and CL_type(s->q) not in <4OrthNet_input>, search\n\
     for q among the 'BestHits' of s; if q was found as the 'n'th best\n\
     subject loci for s, then 'BHn' (n>1) is added as the 3rd column to the\n\
     line that contains s and q in the 1st and 2nd columns.\n\
  - if q was not found among the list of 'BestHits' for s, a line containing\n\
     s, the 1st 'BestHit' of s, 'BHx', tab-delimited, is added.\n\
  - if both CL_type(q->s) and CL_type(s->q) are in <4OrthNet_input>, add 'BH1'\n\
     as the 3rd column to the line containing s and the 1st BestHit of s (==q).\n\
  - for all other lines of 'BestHitPairs' file, the first 2 columns are printed\n\
     with duplicated entries removed.\n\
 3. Output files\n\
  - '-o'|'--Path2BHPs_fxd' defines the output path (dafault='./BHPairs_updtd.1').\n\
  - File names for the updated 'BestHitPairs' are unchanged. \n\
 4. Misc. details:\n\
  - if genes in <4OrthNet_input> and 'BestHitPairs' are in the form of\n\
     '<spcsID>|<geneID>', only the <geneID> portion is used.\n\
  - assumes no header in 'BestHitPairs' files.\n\n\
 by ohdongha@gmail.com 20171225 ver 0.1.3\n"
 
#version_history
#20171225 ver 0.1.3 # if <Project> argument ends with '.list', just ignore it
#20170316 ver 0.1.2 # 'translocated (tlc)' is now called 'transposed (tr)' 
#20160815 ver 0.1.1 # modified to work with 'update_OrthNet_after_mcl.py'.
#20160801 ver 0.1 # print all lines from 'BestHitPairs', with updated preference in the third column.
#20160623 ver 0.0

parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)

parser.add_argument('Project', type=str, help="'./<Project>.list' includes spcsIDs being compared")
parser.add_argument('forOrthNetInput', type=argparse.FileType('r'), help="an output of 'CL_finder_multi_report.py -n'; see below")

parser.add_argument('-b', '--Path2BHPs', dest="Path2BHPs", type=str, default=".", help="PATH to 'BestHitPairs' files; default='.'")
parser.add_argument('-o', '--Path2BHPs_fxd', dest="Path2BHPs_fxd", type=str, default="./BHpairs_updtd.1", help="PATH for updated 'BestHitPairs' files") 

args = parser.parse_args()

path_BestHitPairs = args.Path2BHPs
if path_BestHitPairs[-1] != "/": path_BestHitPairs = path_BestHitPairs + "/"
path_BestHitPairs_fixed = args.Path2BHPs_fxd
if path_BestHitPairs_fixed[-1] != "/": path_BestHitPairs_fixed = path_BestHitPairs_fixed + "/"

try: 
	os.makedirs(path_BestHitPairs_fixed)
except OSError:
	if not os.path.isdir(path_BestHitPairs_fixed): raise

## defining expected format of input and output file names.  modify as needed: 
input_BestHitPairs_filename_format = path_BestHitPairs + "BestHits__%s__vs__%s.list" 
output_BestHitPairs_filename_format = path_BestHitPairs_fixed + "BestHits__%s__vs__%s.list" 


########################################
### 1. reading the species list file ###
########################################
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


########################################
### 2. reading the OrthNetInput file ###
########################################

queryID_gene = ""
subjectID_gene = ""
CL_type = ""

# create bestHit_CLtype_dict; bestHit_CLtype_dict[queryID_gene___subjectID_gene] = CL_type
bestHit_CLtype_dict = dict()

print "\n\nreading %s:" % args.forOrthNetInput.name

for line in args.forOrthNetInput:
	tok = line.split('\t')
	queryID_gene = tok[0].split('|')[ len(tok[0].split('|')) - 1 ].strip()
	subjectID_gene = tok[1].split('|')[ len(tok[1].split('|')) - 1 ].strip()
	CL_type = tok[2].strip()
	bestHit_CLtype_dict[queryID_gene + '___' + subjectID_gene] = CL_type
args.forOrthNetInput.close()


##############################
### 3. starting the update ###
##############################
query_spcsID = ""
subject_spcsID = ""

queryID_gene = ""
queryID_gene_prev = ""
subjectID_gene = ""
subjectID_gene_prev = ""
subjectID_gene_marked = ""
BestHit_rank = 0

BestHitPair_processed = False
BestHitPair_written = False
FirstEntry = True

for i in range( len(spcsID_list) ):
	query_spcsID = spcsID_list[i]
	subject_spcsID = ""
	
	print "\n\nprocessing %s as query_spcsID, comparing with all other species:" % query_spcsID

	for j in range( len(spcsID_list) ):
		if spcsID_list[j] != query_spcsID:
			subject_spcsID = spcsID_list[j]
			
			## opening the BestHitPairs files for reading and writing
			fin_BestHitPairs = open(input_BestHitPairs_filename_format % (query_spcsID, subject_spcsID), "rU")
			fout_BestHitPairs = open(output_BestHitPairs_filename_format % (query_spcsID, subject_spcsID), "w")
			print "updating %s and writing to %s," % (fin_BestHitPairs.name, fout_BestHitPairs.name)
			
			for line in fin_BestHitPairs:
				tok = line.split('\t')
				queryID_gene = tok[0].split('|')[ len(tok[0].split('|')) - 1 ].strip()
				subjectID_gene = tok[1].split('|')[ len(tok[1].split('|')) - 1 ].strip()
				if queryID_gene != queryID_gene_prev:
					# process if there is no better 'BestHit' for the previous queryID_gene
					# in this case, the first BestHit will be repeated at the end, with the 'BHx' tag 
					if BestHitPair_written == False and FirstEntry == False:
						fout_BestHitPairs.write(queryID_gene_prev + '\t' + subjectID_gene_marked + '\t' + 'BHx' + '\n')
					if FirstEntry == True:
						FirstEntry = False
						
					subjectID_gene_marked = subjectID_gene
					BestHit_rank = 1
					BestHitPair_processed = False
					BestHitPair_written = False
					# decision making at a new query entry
					if bestHit_CLtype_dict.get(queryID_gene + '___' + subjectID_gene, "na") != 'tr':
						BestHitPair_processed = True
					elif bestHit_CLtype_dict.get(subjectID_gene + '___' + queryID_gene, "na") != 'na':
						BestHitPair_processed = True
				elif subjectID_gene != subjectID_gene_prev and BestHitPair_processed == False:
					BestHit_rank = BestHit_rank + 1
					if bestHit_CLtype_dict.get(subjectID_gene + '___' + queryID_gene, "na") != 'na':
						subjectID_gene_marked = subjectID_gene
						BestHitPair_processed = True

				# mark the updated BestHitPair
				if BestHitPair_processed == True and BestHitPair_written == False:
					fout_BestHitPairs.write(queryID_gene + '\t' + subjectID_gene_marked + '\tBH' + str(BestHit_rank) + '\n')
					BestHitPair_written = True
				elif subjectID_gene != subjectID_gene_prev or queryID_gene != queryID_gene_prev:
					fout_BestHitPairs.write(queryID_gene + '\t' + subjectID_gene + '\n')
					
				queryID_gene_prev = queryID_gene
				subjectID_gene_prev = subjectID_gene
				
			fin_BestHitPairs.close()
			fout_BestHitPairs.close()
			
print "\ndone\n"