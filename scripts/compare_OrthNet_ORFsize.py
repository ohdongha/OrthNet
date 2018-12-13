#!/usr/bin/env python
import os, sys, numpy, math, argparse
from argparse import RawTextHelpFormatter

		
###################################################
### 0. script description and parsing arguments ###
###################################################

synopsis1 = "\
  report number of nodes in OrthNets that contain complete ORFs;\n\
  report the median, mean, and stdev of ORF sizes in each OrthNet."

synopsis2 = "detailed description:\n\
 0. Pre-requisite: an output file from 'combine_CLfm_result.py -O'\n\
 1. Input files and parameters:\n\
  - './<Project>.list' includes all species IDs (spcsIDs), one per each line.\n\
  - <Input>: the output file of a 'combine_CLfm_result.py -O' run or any file\n\
     with the 'spcs', 'geneID' (if '-O' is on), 'CDS_l', 'OrthNetID' and \n\
     '%mdCDS_l' columns, tab-delimited,\n\
  - '-m min_p_mdCDSlen': minimum '%mdCDS_l' value to decide a 'complete ORF';\n\
     expect a zero or a positive integer; default=50;\n\
 2. Output and options\n\
  - '-o output': output files are '<output>.nodeCounts' and '<output>.mclOutput'\n\
     ;default='<ProjectID>.OrthNet.cORFs_<min_p_mdCDSlen>'\n\
  - for all OrthNets, the number of nodes with complete ORFs (defined by '-m')\n\
     as well as the median, mean, and stdev of ORF lengths for all nodes with\n\
     complete ORFs, are printed to the .nodeCounts output file,\n\
  - '-O | --mclOutput' : create also .mclOutput file, which can be annotated by\n\
     'parse_mclOutput.py -p _internal_ -a ...'; see 'parse_mclOutput.py -h' for\n\
     details,\n\
  - Use '-m 0' to print out values for all nodes.  In this case, the default\n\
     <output> is '<ProjectID>.OrthNet.all'\n\
 by ohdongha@gmail.com ver0.0.2 20180726\n"
 
#version_history
#20180726 ver 0.0.2 another minor fix; OrthNetID "na" is not counted as a valid OrthNetID
#20180423 ver 0.0.1 minor fix for backward compatibility (look for "CDS_len" columns if "CDS_l" doesn't exist
#20180329 ver 0.0

parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)

## positional arguments
parser.add_argument('Project', type=str, help="'./<Project>.list' includes spcsIDs being compared")
parser.add_argument('Input', type=str)

## options and parameters
parser.add_argument('-m', dest="min_p_mdCDSlen", type=int, default=50, help="see below")
parser.add_argument('-o', dest="output", type=str, default="_ND_", help="see below")
parser.add_argument('-O', '--mclOutput', action="store_true", default=False, help="see below")

args = parser.parse_args()


#######################################################
### 1. reading the list file and process parameters ###
#######################################################
# reading the list file
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

# processing the parameters
if args.output == "_ND_":
	if args.min_p_mdCDSlen == 0:
		output_filename = projectID + ".OrthNet.all"
	else:
		output_filename = projectID + ".OrthNet.cORF_" + str( args.min_p_mdCDSlen )
else:
	output_filename = args.output

	
#################################
### 2. reading the input file ###
#################################
fin = open(args.Input, 'r')
print "\nreading input file %s" % fin.name		
header = True

colIndex_spcs = 0
colIndex_CDSl = 0
colIndex_ONid = 0
colIndex_pmCl = 0

spcsID = ""
CDS_l = 0
OrthNetID = ""
pmCl = 0

cORFcounts_OrthNet_dict = dict() # key = OrthNetID, value = a dictionary with spcsID as key and num_cORF as value
cORFlen_OrthNet_dict = dict() # key = OrthNetID, value = a list of the length of all cORFs in the OrthNet
OrthNetID_set = set()

if args.mclOutput:
	colIndex_gene = 0
	geneID = ""
	geneID_OrthNet_dict = dict() # key = OrthNetID, value = a list of geneIDs with cORFs in the OrthNet

for line in fin:
	tok = line.split('\t')
	# parsing the input file header
	if header:
		colIndex_spcs = tok.index("spcs")
		try:
			colIndex_CDSl = tok.index("CDS_l")
		except ValueError:
			colIndex_CDSl = tok.index("CDS_len")
		colIndex_ONid = tok.index("OrthNetID")
		colIndex_pmCl = tok.index("%mdCDS_l")
		if args.mclOutput:
			colIndex_gene = tok.index("geneID")
		header = False
	else:
		try:
			spcsID = tok[ colIndex_spcs ].strip()
			CDS_l = int(tok[ colIndex_CDSl ].strip())
			OrthNetID =  tok[ colIndex_ONid ].strip()
			if tok[ colIndex_pmCl ].strip() != "":
				pmCl =  int(tok[ colIndex_pmCl ].strip())
		except ValueError:
			print "CDS_l %s or " % tok[ colIndex_CDSl ] + "%" + \
					"mdCDS_l %s is not an integer?" % tok[ colIndex_pmCl ].strip() 
		except IndexError:
			print "missing columns in: %s" % line.strip()
		
		if OrthNetID != "" and OrthNetID != "na":
			if args.mclOutput:
				geneID = spcsID + '|' + tok[ colIndex_gene ].strip()
			
			if pmCl >= args.min_p_mdCDSlen:
				if OrthNetID not in cORFcounts_OrthNet_dict:
					cORFcounts_OrthNet_dict[OrthNetID] = dict()
					cORFlen_OrthNet_dict[OrthNetID] = []
					if args.mclOutput:
						geneID_OrthNet_dict[OrthNetID] = []
				cORFcounts_OrthNet_dict[OrthNetID][spcsID] = cORFcounts_OrthNet_dict[OrthNetID].get(spcsID, 0) + 1
				cORFlen_OrthNet_dict[OrthNetID].append(CDS_l)
				if args.mclOutput:
					geneID_OrthNet_dict[OrthNetID].append(geneID)
					
			OrthNetID_set.add(OrthNetID)

print "finished reading input file."
fin.close()


#########################
### 3. writing output ###
#########################
fout = open(output_filename + ".nodeCounts", 'w')
print "\ncollecting nodes with %" + "mdCDS_l >= %d," % args.min_p_mdCDSlen
print "writing "".nodeCounts"" output to %s," % fout.name

if args.mclOutput:
	fout_mclOutput = open(output_filename + ".mclOutput", 'w')
	print "\nwriting "".mclOutput"" output to %s" % fout_mclOutput.name

#constructing and writing the header
if args.min_p_mdCDSlen == 0:
	header_item = ["ONid", "#nd_all" ]
else:
	header_item = ["ONid", "#nd_cORF_" + str( args.min_p_mdCDSlen ) ]
header_item.append('.'.join(spcsID_list))

for sID in spcsID_list:
	header_item.append(sID)
	
header_item = header_item + ["cORFl_av", "cORFl_md", "cORFl_sd"]
fout.write('\t'.join(header_item) + '\n')

# writing the content for both .nodeCounts and .mclOutput files
for ONid in sorted( OrthNetID_set ):
	#initializing 
	line_item = []
	nodeCounts_list = []
	
	if ONid in cORFcounts_OrthNet_dict:
		line_item = [ONid, str( len( cORFlen_OrthNet_dict[ONid] ) ) ] 
		for sID in spcsID_list:
			nodeCounts_list.append( str(cORFcounts_OrthNet_dict[ONid].get(sID, 0) ) )
		line_item.append( '.'.join(nodeCounts_list) )
		line_item.append( '\t'.join(nodeCounts_list) )
		line_item.append( str( int( math.ceil(numpy.mean( cORFlen_OrthNet_dict[ONid] ) ) ) ) )
		line_item.append( str( int( math.ceil(numpy.median( cORFlen_OrthNet_dict[ONid] ) ) ) ) )
		line_item.append( str( int( math.ceil(numpy.std( cORFlen_OrthNet_dict[ONid] ) ) ) ) )
		#line_item.append( str( cORFlen_OrthNet_dict[ONid] ) ) # for debugging purpose

	else:
		line_item = [ONid, '0' ] 
		for sID in spcsID_list:
			nodeCounts_list.append( '0' )
		line_item.append( '.'.join(nodeCounts_list) )
		line_item.append( '\t'.join(nodeCounts_list) )
		line_item = line_item + [ "0", "0", "0" ]
	
	fout.write('\t'.join( line_item ) + '\n')
	if args.mclOutput:
		if ONid in cORFcounts_OrthNet_dict:
			fout_mclOutput.write( ONid + ": " + ' '.join( geneID_OrthNet_dict[ONid] ) + '\n')
		else:
			fout_mclOutput.write( ONid + ":\n")

fout.close()
if args.mclOutput:
	fout_mclOutput.close()
print "\ndone\n"
