#!/usr/bin/env python
import os, sys, subprocess, math, numpy, argparse
from argparse import RawTextHelpFormatter

		
###################################################
### 0. script description and parsing arguments ###
###################################################

synopsis1 = "\
  combine CLfinder results for multiple genomes into a single file.\n"

synopsis2 = "detailed description:\n\
 0. Pre-requisite: a 'CL_finder_multi.py' or 'update_OrthNet_afterMCL.py' run,\n\
 1. Input files and parameters:\n\
  - './<Project>.list' includes all species IDs (spcsIDs), one per each line.\n\
  - '-p Path2CLfm_output': path to the 'CL_finder_multi.py' output folder\n\
     that contains 'spcsID.CL_compared2*.<format>.txt' files for all species,\n\
     created by a 'CL_finder_multi.py -r' run; defaults = './'\n\
 2. Output and options\n\
  - for all genome pairs, CLfinder results, including any additional information\n\
     added by 'update_OrthNet_afterMCL.py' etc., are combined to tab-delimited\n\
     <output> file,\n\
  - when the query and target genomes are the same, cells are left empty\n\
  - '-F CLfm_nameFmt': the CLfinder output file name format string, after the\n\
     '<Query spcsID>.CL_compared2<Target spcsID1/2/3/.../n>' part;\n\
     default="".20.4.20.txt""\n\
  - '-O'|'--ORFsize': add three more columns, median and stdev of ortholog CDS length\n\
     (mdCDS_l and sdCDS_l) and proportion of CDS_l compared to the median (%mdCDS_l),\n\
	 for each query gene; all geneIDs should be unique; default=False\n\
 by ohdongha@gmail.com ver0.1 20180307\n"
 
#version_history
#20180307 ver 0.1 '-O' option added to compare ORF sizes among BestHits,
#20180306 ver 0.0

parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)

## positional arguments
parser.add_argument('Project', type=str, help="'./<Project>.list' includes spcsIDs being compared")
parser.add_argument('Output', type=str)

## options to receive PATH to input files
parser.add_argument('-p', dest="Path2CLfm_output", type=str, default=".", help="see below")

## other options
parser.add_argument('-F', dest="CLfm_nameFmt", type=str, default=".20.4.20.txt", help="see below")
parser.add_argument('-O', '--ORFsize', action="store_true", default=False, help="see below")
args = parser.parse_args()

# defining PATH and the temporary file name
path_CLfm_output = args.Path2CLfm_output
if path_CLfm_output[-1] != "/": path_CLfm_output = path_CLfm_output + "/"


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


###########################################################
### 2. creating input file names and output file header ###
###########################################################
# assembling input file names
print "\nCombining following files to %s" % args.Output
CLfm_outputPath_list = []
n = 0

for query_spcsID in spcsID_list:
	CLfm_outputPath_list.append( path_CLfm_output + query_spcsID + ".CL_compared2" )
	for target_spcsID in spcsID_list:
		if query_spcsID != target_spcsID:
			CLfm_outputPath_list[n] = CLfm_outputPath_list[n] + target_spcsID
	CLfm_outputPath_list[n] = CLfm_outputPath_list[n] + args.CLfm_nameFmt
	if os.path.isfile( CLfm_outputPath_list[n] ):
		print CLfm_outputPath_list[n]
	else:
		print "error: %s does not exist" % CLfm_outputPath_list[n]
	n += 1
	
# creating the header string by modifying the header of the first CLfinder output file
with open( CLfm_outputPath_list[0] ) as f:
	header = f.readline().strip()
	colName_list = header.split('\t')

CLfm_param_string = colName_list[14].split('_')[-1] # e.g. CLfm_param_string = "CL.20.3.20"
colName_list.insert(14, spcsID_list[0] + "TDid")
colName_list.insert(14, spcsID_list[0] + "locusID")
colName_list.insert(14, spcsID_list[0] + "_PG")
colName_list.insert(14, spcsID_list[0] + "_geneID")
colName_list.insert(14, "vs" + spcsID_list[0] + "_" + CLfm_param_string)
colName_list.insert(0, "spcs")
colName_list[12] = "PG"
colName_list[13] = "locusID"
colName_list[14] = "TDid"

for i in range(15, len(colName_list)):
	colName_list[i] = colName_list[i].replace(spcsID_list[0] + '-', "vs", 1)

print "\n%s will contain the following columns:" % args.Output
for i in range(0, len(colName_list)):
	print"%d. %s" % ( i+1, colName_list[i])
	
new_header = '\t'.join(colName_list)


#####################################
## 3. read input and print Output ###
#####################################
# open output file (or temp file in case of -O option)
if args.ORFsize:
	fout = open(args.Output + '.' + projectID + '.temp', "w")
else:
	fout = open(args.Output, "w")

# print header
fout.write(new_header + '\n')

# start reading input files and print to the output
colNum2insert = 14

# with -O option, read in CDS_l for all geneIDs
if args.ORFsize:
	CDSlen_dict = dict() # key = geneID, value = cds_l

for i in range(0, len(CLfm_outputPath_list)):
	fin = open( CLfm_outputPath_list[i], 'r')
	header = True
	for line in fin:
		if header:
			header = False
		else:
			tok = line.strip().split('\t')
			for j in range(0,5):
				tok.insert(colNum2insert, "")
			fout.write(spcsID_list[i] + '\t' + '\t'.join(tok) + '\n')
			
			# with -O option, read in CDS_l for all geneIDs 
			if args.ORFsize:
				if tok[0] not in CDSlen_dict:
					CDSlen_dict[ tok[0] ] = int( tok[10] )
				else:
					print "warning: geneID %s is not unique!" % tok[0]
					
	fin.close()
	colNum2insert += 5

fout.close()	


#################################
## 4. process the '-O' option ###
#################################
if args.ORFsize:
	print "\nAdding ORF size comparison results among orthologs:"
	fin2 = open(args.Output + '.' + projectID + '.temp', "r")
	fout2 = open(args.Output, "w")
	header = True
	colNum4orthologGeneID = 16 
	query_geneID = ""
	query_CDSlen = 0
	ortholog_geneID = ""
	ortholog_CDSlen_list = []
	
	for line in fin2:
		if header:
			fout2.write(line.strip() + "\tmdCDS_l\tsdCDS_l\t%mdCDS_l\n")
			header = False
		else:
			# initializing
			colNum4orthologGeneID = 16 
			ortholog_CDSlen_list = []
			tok = line.split('\t')
			query_geneID = tok[1]
			query_CDSlen = CDSlen_dict[ query_geneID ]
			# scanning 
			for i in range(0, len(spcsID_list) ):
				ortholog_geneID = tok[ colNum4orthologGeneID ]
				if ortholog_geneID != "" and ortholog_geneID != "na":
					ortholog_CDSlen_list.append( CDSlen_dict[ ortholog_geneID ] )
				colNum4orthologGeneID += 5
			# printing
			if len(ortholog_CDSlen_list) > 0:
				fout2.write(line.strip() + "\t%d\t%d\t%d\n" % ( \
						int( math.ceil(numpy.median(ortholog_CDSlen_list)) ), \
						int( math.ceil(numpy.std(ortholog_CDSlen_list)) ), \
						int( math.ceil( query_CDSlen / numpy.median(ortholog_CDSlen_list) * 100.0 ) ) ) )
			else:
				fout2.write(line.strip() + "\t\t\t\n")
	fin2.close()
	fout2.close()
	
	#cleaning the temporary file
	try:
		os.remove(args.Output + '.' + projectID + '.temp')
	except OSError:
		pass
	
print "\ndone\n"
