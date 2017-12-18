#!/usr/bin/env python
import sys, os, subprocess, datetime, argparse
from argparse import RawTextHelpFormatter


###################################################
### 0. script description and parsing arguments ###
###################################################
synopsis1 = "\
  parses mRNA and CDS features from a .gtf file and create a .txt summary with\n\
  one gene loci per line.  works best with protein-coding (pc) gene loci"
synopsis2 = "detailed description:\n\
 - from <input.gtf>, find out coordinates of mRNA, CDS, number of exons for\n\
    each transcript model and print to <output.txt>, tab-delimited.\n\
 - for each 'transcript_id' recoreds in the 9th column of <input.gtf>, print\n\
    transcript_id, Chr, Str, mRNA_start(s), mRNA_end(e), #exon_mRNA,mRNA_len(l),\n\
    CDS_s, CDS_e, #exon_CDS, CDS_l, to <output.txt>\n\
 - <output.txt> is sorted based on mRNA_start and then Chr, alphanumerically.\n\
    for the same mRNA_s, entries with longer CDS_l and mRNA_l come first.\n\
 - '-c'|'--collapse': remove gene loci whose coordinates identical with or\n\
    nested in another gene locus. report collapese loci to STDOUT.\n\
 - '-m MARGIN': gene loci coordinates different by less than MARGIN nucleotides\n\
    are considered identical. assume '-c' automatically.\n\
 - '-r'|'--report': report transcript_id overlapping with others\n\
    (e.g. isoforms) to STDOUT.\n\
by ohdongha@gmail.com 20171017 ver 0.3\n\n"

#version_history
#20171017 ver 0.3 # added -c -r options
#20171007 ver 0.2 # added a step to detect and report overlapping entries (e.g. isoforms)
#20170316 ver 0.1 # added a step to sort the output
#20151210 ver 0.0 

#### Example: lines in a <input.gtf> file
##ch1-1   FGenesh++       CDS     24041   24116   4.49    +       0       transcript_id "Tp1g00060";
##ch1-1   transdecoder    exon    24316   25212   .       -       .       transcript_id "Tp1g00070.a.x";
##ch1-1   FGenesh++       CDS     24336   24919   36.22   +       2       transcript_id "Tp1g00060";
##ch1-1   transdecoder    CDS     25153   25212   .       -       0       transcript_id "Tp1g00070.a.x";
##ch1-1   transdecoder    CDS     25285   25500   .       -       0       transcript_id "Tp1g00070.a.x";
##ch1-1   transdecoder    exon    25285   25500   .       -       .       transcript_id "Tp1g00070.a.x";
##ch1-1   transdecoder    CDS     25584   25704   .       -       1       transcript_id "Tp1g00070.a.x";
##ch1-1   transdecoder    exon    25584   25704   .       -       .       transcript_id "Tp1g00070.a.x";
##ch1-1   transdecoder    CDS     25825   26155   .       -       2       transcript_id "Tp1g00070.a.x";
##ch1-1   transdecoder    exon    25825   26155   .       -       .       transcript_id "Tp1g00070.a.x";
##ch1-1   transdecoder    CDS     26541   26922   .       -       0       transcript_id "Tp1g00070.a.x";
##ch1-1   transdecoder    exon    26541   27023   .       -       .       transcript_id "Tp1g00070.a.x";
##ch1-1   FGenesh++       CDS     26686   26880   27.57   +       0       transcript_id "Tp1g00080";
##ch1-1   FGenesh++       CDS     27550   27797   10.96   -       2       transcript_id "Tp1g00090";
##ch1-1   FGenesh++       CDS     28177   28270   13.08   -       0       transcript_id "Tp1g00090";
##ch1-1   FGenesh++       CDS     28633   29832   2009.22 -       0       transcript_id "Tp1g00100";
##ch1-1   FGenesh++       CDS     29907   30257   634.46  -       0       transcript_id "Tp1g00100";
##...

parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)

# positional parameters
parser.add_argument('input_gtf', type=argparse.FileType('r'))
parser.add_argument('outfile', type=argparse.FileType('w'))

# options
parser.add_argument('-c', '--collapse', action="store_true", default=False)
parser.add_argument('-m', dest="margin", type=int, default= -1) 
parser.add_argument('-r', '--report_overlap', action="store_true", default=False)

args = parser.parse_args()
outfile_name = args.outfile.name


#################################
### 1. reading in <input.gtf> ###
#################################
chr_dict = dict()
str_dict = dict()

mRNA_start_dict = dict()
mRNA_end_dict = dict()
mRNA_nExon_dict = dict()
mRNA_len_dict = dict()

CDS_start_dict = dict()
CDS_end_dict = dict()
CDS_nExon_dict = dict()
CDS_len_dict = dict()

chr = ""
type = ""
start = 0
end = 0 
strand = ""
ninthColumn_records = ""
geneID = ""

nGene = 0
nCDS = 0
newline_accepted = 0

print "reading %s as the <input.gtf>:" % sys.argv[1]

for line in args.input_gtf:
	newline_accepted = 0
	tok = line.replace('\"','').split('\t')
	try:
		chr = tok[0]
		type = tok[2]
		start = int(tok[3])
		end = int(tok[4])
		strand = tok[6]
		ninthColumn_records = tok[8].split(';')
		for record in ninthColumn_records:
			if record.strip().split(' ')[0] == 'transcript_id':
				geneID = record.strip().split(' ')[1]
				newline_accepted = 1
	except (ValueError, IndexError) :
		print "There is a non-valid line."
		newline_accepted = 0
	if newline_accepted ==1 and type == "exon" :
		if geneID not in mRNA_start_dict:
			chr_dict[geneID] = chr
			str_dict[geneID] = strand
			mRNA_start_dict[geneID] = start
			mRNA_end_dict[geneID] = end
			mRNA_len_dict[geneID] = end - start + 1
			mRNA_nExon_dict[geneID] = 1
			nGene = nGene + 1
		else:
			mRNA_start_dict[geneID] = min(start, mRNA_start_dict[geneID])
			mRNA_end_dict[geneID] = max(end, mRNA_end_dict[geneID])
			mRNA_len_dict[geneID] = mRNA_len_dict[geneID] + end - start + 1
			mRNA_nExon_dict[geneID] = mRNA_nExon_dict[geneID] + 1
	elif newline_accepted ==1 and type == "CDS" :
		if geneID not in CDS_start_dict:
			chr_dict[geneID] = chr
			str_dict[geneID] = strand
			CDS_start_dict[geneID] = start
			CDS_end_dict[geneID] = end
			CDS_len_dict[geneID] = end - start + 1
			CDS_nExon_dict[geneID] = 1
			nCDS = nCDS + 1			
		else:
			CDS_start_dict[geneID] = min(start, CDS_start_dict[geneID])
			CDS_end_dict[geneID] = max(end, CDS_end_dict[geneID])
			CDS_len_dict[geneID] = CDS_len_dict[geneID] + end - start + 1
			CDS_nExon_dict[geneID] = CDS_nExon_dict[geneID] + 1

print "%d gene models with 'exon' records and %d with 'CDS' records were found in %s.\n" % (nGene, nCDS, sys.argv[1])
args.input_gtf.close()


###############################
### 2. writing <output.txt> ###
###############################
print "writing to %s:" % outfile_name
args.outfile.write("geneID\tChr\tStr\tmRNA_s\tmRNA_e\t#exon_mRNA\tmRNA_l\tCDS_s\tCDS_e\t#exon_CDS\tCDS_l\n")

for key in sorted(chr_dict):
	try:
		if key in CDS_start_dict:		
			args.outfile.write( key + '\t' + \
						chr_dict[key] + '\t' + \
						str_dict[key] + '\t' + \
						## if no records for mRNA, copy records from CDS
						str( mRNA_start_dict.get(key, CDS_start_dict[key]) ) + '\t' + \
						str( mRNA_end_dict.get(key, CDS_end_dict[key]) ) + '\t' + \
						str( mRNA_nExon_dict.get(key, CDS_nExon_dict[key]) ) + '\t' + \
						str( mRNA_len_dict.get(key, CDS_len_dict[key]) ) + '\t' + \
						str( CDS_start_dict[key] ) + '\t' + \
						str( CDS_end_dict[key] ) + '\t' + \
						str( CDS_nExon_dict[key] ) + '\t' + \
						str( CDS_len_dict[key] ) + '\n' )
		else:
			args.outfile.write( key + '\t' + \
						chr_dict[key] + '\t' + \
						str_dict[key] + '\t' + \
						## if no records for mRNA, copy records from CDS
						str( mRNA_start_dict[key] ) + '\t' + \
						str( mRNA_end_dict[key] ) + '\t' + \
						str( mRNA_nExon_dict[key] ) + '\t' + \
						str( mRNA_len_dict[key] ) + '\t' + \
						## if no records for CDS, assume non-coding gene model
						"NA" + '\t' + \
						"NA" + '\t' + \
						"NA" + '\t' + \
						"0" + '\n' )
	except (ValueError, IndexError) :
		print "Something bad just happened while writing, please troubleshoot :p"
#	except KeyError :
#		print key		
args.outfile.close()

## sort the output file
print "sorting %s:" % outfile_name
i = datetime.datetime.now()
temp_filename = "parse_gtf_temp_" + i.strftime('%Y%m%d_%H%M%S')
subprocess.call("awk 'NR == 1; NR > 1 {print $0 | \"sort -k2,2 -k4,4n -k7,7nr -k11,11nr\"}' " + outfile_name + " > " + temp_filename , shell=True)
subprocess.call("mv " + temp_filename + " " + outfile_name , shell=True)


#############################
### 3. optional processes ###
#############################

## with -c option
if args.collapse == True or args.margin != -1:
	print "\ncollapsing gene loci with coordinates identical with or nested in other locus in %s:" % outfile_name
	print "output file before collapsing moved to %s:" % (outfile_name + "_b4collapsing")
	subprocess.call("cp " + outfile_name + " " + outfile_name + "_b4collapsing", shell=True)

	#initializing
	margin = max(0, args.margin) # if only "-c" was used, consider only exact match 
	i = datetime.datetime.now()
	temp_filename = "parse_gtf_collapsing_temp_" + i.strftime('%Y%m%d_%H%M%S')
	fin_output = open(outfile_name, "r")
	fout_collapsed = open(temp_filename, "w")

	prev_line = ""
	prev_chr = ""
	prev_str = ""
	prev_end = 0
	chr = ""
	str = ""
	mRNA_end = 0
	num_removed = 0
	num_line = 0
	num_line_removed = 0

	# reading sorted output file (fin_output) and writing to temp_filename
	for line in fin_output:
		num_line += 1
		try:
			tok = line.split('\t')
			chr = tok[1].strip()
			str = tok[2].strip()
			mRNA_start = int(tok[3].strip())
			mRNA_end = int(tok[4].strip())
	
			if (chr == prev_chr) and (str == prev_str) and (mRNA_start < prev_end) and (mRNA_end <= (prev_end + margin)):
				num_removed += 1
				if (num_line - 1) > num_line_removed: # report the line kept
					print prev_line.strip() + "\tline_%d" % (num_line - 1)
				if mRNA_end > prev_end: # mark lines removed by "-m" option
					print line.strip() + "\tline_%d_removed_-m" % num_line
				else:
					print line.strip() + "\tline_%d_removed" % num_line
				num_line_removed = num_line
			else:
				fout_collapsed.write(line)
	
			prev_line = line
			prev_chr = chr
			prev_str = str
			prev_end = mRNA_end
	
		except (ValueError, IndexError) :
			if num_line > 1 :
				print "line %d appears invalid" % num_line
			else:
				fout_collapsed.write(line)

	print "\nremoved %d identical/nested entries" % num_removed

	fin_output.close()
	fout_collapsed.close()
	subprocess.call("mv " + temp_filename + " " + outfile_name , shell=True)

## with -r option 
if args.report_overlap == True:
	print "\ndetecting overlapping transcripts (e.g. isoforms) in %s:" % outfile_name

	#initializing	
	fin_output = open(outfile_name, "r")
	prev_chr = ""
	prev_str = ""
	prev_end = 0
	chr = ""
	str = ""
	mRNA_start = 0
	mRNA_end = 0
	num_overlap = 0
	num_line = 0
	
	for line in fin_output:
		num_line += 1
		try:
			tok = line.split('\t')
			chr = tok[1].strip()
			str = tok[2].strip()
			mRNA_start = int(tok[3].strip())
			mRNA_end = int(tok[4].strip())
	# if you want to use CDS coordinates to detect overlaps
	#		mRNA_start = int(tok[7].strip())
	#		mRNA_end = int(tok[8].strip())
			
			if (chr == prev_chr) and (str == prev_str) and ( mRNA_start < prev_end):
				num_overlap += 1
				print line.strip() + "\toverlap_size: %d" % (prev_end - mRNA_start + 1)
				
			prev_chr = chr
			prev_str = str
			prev_end = mRNA_end
	
		except (ValueError, IndexError) :
			if num_line > 1 :
				print "line %d appears invalid" % num_line
	
	print "\nfound %d overlapping entries" % num_overlap
	
	fin_output.close()

print "all done\n"