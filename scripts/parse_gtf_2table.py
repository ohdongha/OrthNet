#!/usr/bin/env python
import sys, os, math, subprocess, datetime, argparse
from argparse import RawTextHelpFormatter


###################################################
### 0. script description and parsing arguments ###
###################################################
synopsis1 = "\
  parses mRNA and CDS features from a .gtf file and create a .txt summary with\n\
  one gene loci per line.  works best with protein-coding (pc) gene loci"
synopsis2 = "detailed description:\n\
 1. Input files and options:\n\
  - './<Project>.list' including all species IDs (spcsIDs), one per line,\n\
     that appear in OrthNets.\n\
 1. Parsing .gtf file:\n\
  - from <input.gtf>, find out coordinates of mRNA, CDS, number of exons for\n\
     each transcript model and print to <output.txt>, tab-delimited.\n\
  - for each 'transcript_id' record in the 9th column of <input.gtf>, print\n\
     transcript_id, Chr, Str, mRNA_start(s), mRNA_end(e), #exon_mRNA,mRNA_len(l),\n\
     CDS_s, CDS_e, #exon_CDS, CDS_l, to <output.txt>\n\
  - <output.txt> is sorted based on mRNA_start and then Chr, alphanumerically.\n\
     for the same mRNA_s, entries with longer CDS_l and mRNA_l come first.\n\
  - '-g'|'--gene_id': use 'gene_id' record instead of 'transcript_id'\n\
 2. Options for multiple gene models in a locus (ex. isoforms):\n\
  - '-c'|'--collapse': remove gene loci whose coordinates identical with or\n\
     nested in another gene locus. report collapese loci to STDOUT.\n\
  - '-m MARGIN': gene loci coordinates different by less than MARGIN nucleotides\n\
     are considered identical. assume '-c' automatically.\n\
  - '-r'|'--report': report transcript_id overlapping with others\n\
     (e.g. isoforms) to STDOUT.\n\
  - '-l'|'--cluster': instead of collapsing overlapping transcript models\n\
     cluster them and report numerical cluster IDs as the last column; all\n\
     transcript models who overlap and in the same direction\strain, are\n\
     grouped in the same cluster; -c, -r, and -l options are mutually\n\
     exclusive.\n\
  - '-L'|'--LongestORF': after clustering, leave the one with the longest\n\
     ORFs for each cluster; '-l' is assumed and ignores '-c' or '-r'\n\
  - these options may work similar as 'gffread -M'; check 'gffread -h',\n\
 3. Option to filter .gtf file:\n\
  - '-p'|'--protein_coding': print only protein-coding gene models (i.e. with\n\
     CDS records);\n\
  - '-e <transcriptID.list>': given a list of transcriptIDs, one per line,\n\
     print .gtf file containing only those in <transcriptID.list>; <output.txt>\n\
     is the filtered .gtf file, instead of a .gtfParsed.txt file.\n\
by ohdongha@gmail.com 20190409 ver 0.5.2\n\n"

#version_history
#20190409 ver 0.5.2 # -e option takes only the first field as the transcriptID from <transcriptID.list> 
#20190220 ver 0.5.1 # added -p option to print only protein-coding gene models; added -L option to select the one with longest ORF for each cluster detected by -l option 
#20190125 ver 0.5 # added -e option to print .gtf file containing only transcripts listed in an input
#20180816 ver 0.4.1 # fixed a bug with the -l option
#20180622 ver 0.4 # added -l option
#20171223 ver 0.3.1 # non-critical bug fix on STDOUT messages; added -g option
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
parser.add_argument('-g', '--gene_id', action="store_true", default=False)
parser.add_argument('-c', '--collapse', action="store_true", default=False)
parser.add_argument('-m', dest="margin", type=int, default= -1) 
parser.add_argument('-r', '--report_overlap', action="store_true", default=False)
parser.add_argument('-l', '--cluster', action="store_true", default=False)
parser.add_argument('-L', '--Longest_ORF', action="store_true", default=False)
parser.add_argument('-p', '--protein_coding', action="store_true", default=False)
parser.add_argument('-e', dest="tID_list", type=str, default= "")

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

print "reading %s as the <input.gtf>:" % args.input_gtf.name


#############################################################
### 1.1 if a list of transcriptIDs given with '-e' option ###
#############################################################
if args.tID_list != "":
	# 1.5.1 read in the tID.list file
	tID_set = set()
	try:
		fin_tID_list = open(args.tID_list, 'r')
		for line in fin_tID_list:
			tID_set.add(line.split()[0].strip())
	except:
		print "Failed to read in %s, exiting..." % fin_tID_list.name
		sys.exit()
	
	# 1.5.2 print line to outfile if transcriptID is in the set
	for line in args.input_gtf:
		tok = line.replace('\"','').split('\t')
		# read transcript_id (or gene_id if '-g' is on)
		if len(tok) >= 9:
			ninthColumn_records = tok[8].split(';')
			for record in ninthColumn_records:
				if args.gene_id == True:
					if record.strip().split(' ')[0] == 'gene_id':
						geneID = record.strip().split(' ')[1]
						newline_accepted = 1			
				else:
					if record.strip().split(' ')[0] == 'transcript_id':
						geneID = record.strip().split(' ')[1]
						newline_accepted = 1
			# filter the line that is present in the tID.list
			if geneID in tID_set:
				args.outfile.write(line) # that's it!!
		elif line[0] != '#': # ignore header line
			print "an invalid line detected in %s " % args.input_gtf.name # well, ignore anyway ...
else:			
	###########################################
	### 1.2 if '-e' option is off, continue ###
	###########################################
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
				if args.gene_id == True:
					if record.strip().split(' ')[0] == 'gene_id':
						geneID = record.strip().split(' ')[1]
						newline_accepted = 1			
				else:
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
	
	print "## %d gene models with 'exon' records and %d with 'CDS' records were found in %s.\n" % (nGene, nCDS, args.input_gtf.name)
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
			elif not args.protein_coding: # if '-p' option is on, skip those without CDS records
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
		prev_strand = ""
		prev_end = 0
		chr = ""
		strand = ""
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
				strand = tok[2].strip()
				mRNA_start = int(tok[3].strip())
				mRNA_end = int(tok[4].strip())
		
				if (chr == prev_chr) and (strand == prev_strand) and (mRNA_start < prev_end) and (mRNA_end <= (prev_end + margin)):
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
				prev_strand = strand
				prev_end = mRNA_end
		
			except (ValueError, IndexError) :
				if num_line > 1 :
					print "line %d appears invalid" % num_line
				else:
					fout_collapsed.write(line)
	
		print "\n## removed %d identical/nested entries" % num_removed
	
		fin_output.close()
		fout_collapsed.close()
		subprocess.call("mv " + temp_filename + " " + outfile_name , shell=True)
	
	## with -r option
	elif args.report_overlap:
		print "\ndetecting overlapping transcripts (e.g. isoforms) in %s:" % outfile_name
	
		#initializing	
		fin_output = open(outfile_name, "r")
		prev_chr = ""
		prev_strand = ""
		prev_end = 0
		chr = ""
		strand = ""
		mRNA_start = 0
		mRNA_end = 0
		num_overlap = 0
		num_line = 0
		
		for line in fin_output:
			num_line += 1
			try:
				tok = line.split('\t')
				chr = tok[1].strip()
				strand = tok[2].strip()
				mRNA_start = int(tok[3].strip())
				mRNA_end = int(tok[4].strip())
		# if you want to use CDS coordinates to detect overlaps
		#		mRNA_start = int(tok[7].strip())
		#		mRNA_end = int(tok[8].strip())
				
				if (chr == prev_chr) and (strand == prev_strand) and ( mRNA_start < prev_end):
					num_overlap += 1
					print line.strip() + "\toverlap: %d in %d" % ( (prev_end - mRNA_start + 1), (mRNA_end - mRNA_start) )
					
				prev_chr = chr
				prev_strand = strand
				prev_end = mRNA_end
		
			except (ValueError, IndexError) :
				if num_line > 1 :
					print "line %d appears invalid" % num_line
		
		print "\n## found %d overlapping entries" % num_overlap
		
		fin_output.close()
	
	## with -l or -L option 
	elif args.cluster or args.Longest_ORF:
		print "\ndetecting clusters of transcripts in %s, based on genomic locations:" % outfile_name
		
		#initializing	
		fin_output = open(outfile_name, "r")
		temp_filename = "parse_gtf_c_temp_" + i.strftime('%Y%m%d_%H%M%S')
		fout_clustered = open(temp_filename, "w")
		
		digit4cIDs = int(math.log(nCDS/2,10)) + 1
		#print "digit4cIDs = %d" % digit4cIDs
		
		prev_chr = ""
		prev_end_plus = 0
		prev_end_minus = 0
		chr = ""
		strand = ""
		mRNA_start = 0
		mRNA_end = 0
		num_c_plus = 0
		num_c_minus = 0
		num_line = 0
	
		for line in fin_output:
			num_line += 1
			if num_line == 1:
				fout_clustered.write(line.strip() + "\tcID\n")
			else:
				try:
					tok = line.split('\t')
					chr = tok[1].strip()
					strand = tok[2].strip()
					mRNA_start = int(tok[3].strip())
					mRNA_end = int(tok[4].strip())
				except (ValueError, IndexError) :
					print "line %d appears invalid" % num_line
				if chr != prev_chr: # fixed in 0.4.1: initialize if new chr
					prev_end_plus = 0
					prev_end_minus = 0
				if strand == "+":
					if mRNA_start < prev_end_plus:
						prev_end_plus = max( prev_end_plus, mRNA_end ) # extend the end of the cluster
					else:
						num_c_plus += 1
						prev_end_plus = mRNA_end
						prev_chr = chr
					fout_clustered.write(line.strip() + "\tp%s\n" % str(num_c_plus).rjust(digit4cIDs, '0'))
				elif strand == "-":
					if mRNA_start < prev_end_minus:
						prev_end_minus = max( prev_end_minus, mRNA_end ) # extend the end of the cluster
					else:
						num_c_minus += 1
						prev_end_minus = mRNA_end
						prev_chr = chr
					fout_clustered.write(line.strip() + "\tm%s\n" % str(num_c_minus).rjust(digit4cIDs, '0'))
				else: # if neither + or - strand, just print the line
					fout_clustered.write(line)
		print "## identified %d and %d clusters in the plus and minus strand, respectively," % (num_c_plus, num_c_minus)	
		fin_output.close()
		fout_clustered.close()
		
		if args.Longest_ORF:
			# this will leave only the line with the longest CDS per each cluster;
			subprocess.call("sort -k11,11nr " + temp_filename + " | sort -k12,12 -u - | awk 'NR == 1; NR > 1 {print $0 | \"sort -k2,2 -k4,4n\"}' - | cut -f1-11 - > " + outfile_name, shell = True )
			subprocess.call("rm " + temp_filename , shell=True)
		else:
			subprocess.call("mv " + temp_filename + " " + outfile_name , shell=True)
		
print "all done\n"