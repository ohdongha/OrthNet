#!/usr/bin/env python
import os, sys
import argparse
from argparse import RawTextHelpFormatter

###################################################
### 0. script description and parsing arguments ###
###################################################

synopsis1 = "\
  convert a BLAST tabular output into a table of one query-subject pair per line\n\
  by consolidating HSPs and calculate proportion of query and subject covered by\n\
  HSPs (COV) and approximate proportion of identical sequences within HSPs (IDN).\n"

synopsis2 = "detailed description:\n\
 1. Input:\n\
  - <input> is a tabulated blast+ output using -outfmt '6 std qlen slen'\n\
 2. Output:\n\
  - <output> contains the following for each query-subject pair, tab-delimited:\n\
     query(q), subject(s), num_HSPs, qHSP_nt, qIDN_nt, qHSP_cov, qHSP_idn,\n\
     sHSP_nt, sIDN_nt, sIDN_cov, and sHSP_idn,\n\
  - total_sc == sum of blast hit scores for all HSPs,\n\
  - HSP_nt == total nucleotides (nt) in HSPs,\n\
  - IDN_nt == approximate identical nt, i.e., HSP_nt * average%IDN / 100,\n\
  - HSP_ovl == length of overlap among HSPs,\n\
  - qHSP_cov == qHSP_nt / qlen; sHSP_cov == sHSP_nt / slen,\n\
  - HSP_idn == IDN_nt / HSP_nt,\n\
 3. Options to filter results:\n\
  - '--min_qHSP_cov': minimum qHSP_cov to keep (0~1), default==0,\n\
  - '--min_sHSP_cov': minimum sHSP_cov to keep (0~1), default==0,\n\
  - '--min_qHSP_idn': minimum qHSP_idn to keep (0~1), default==0,\n\
  - '--min_sHSP_idn': minimum sHSP_idn to keep (0~1), default==0,\n\
 4. Misc:\n\
  - for BLASTP results, '_nt' becomes '_aa (amino acids)',\n\
  - ignores strands of HSPs and calculates just the total coverages.\n\
 by ohdongha@gmail.com ver0.0.2 20180105\n"
 
#version_history
#20180103 ver 0.0.1 added total_sc, qHSP_ovl, sHSP_ovl; deal with HSPs in reverse strand of subject,
#20180102 ver 0.0

parser = argparse.ArgumentParser(description = synopsis1, epilog = synopsis2, formatter_class = RawTextHelpFormatter)

## positional arguments
parser.add_argument('input', type=argparse.FileType('r'))
parser.add_argument('output', type=argparse.FileType('w'))

## options to filter results
parser.add_argument('--min_qHSP_cov', dest="min_qHSP_cov", type=float, default=0.0)
parser.add_argument('--min_sHSP_cov', dest="min_sHSP_cov", type=float, default=0.0)
parser.add_argument('--min_qHSP_idn', dest="min_qHSP_idn", type=float, default=0.0)
parser.add_argument('--min_sHSP_idn', dest="min_sHSP_idn", type=float, default=0.0)

args = parser.parse_args()

min_qHSP_cov = args.min_qHSP_cov
min_sHSP_cov = args.min_sHSP_cov
min_qHSP_idn = args.min_qHSP_idn
min_sHSP_idn = args.min_sHSP_idn


##############################################
### 1. reading, consolidating, and writing ###
##############################################
first_line = True
output_line = list()
output_header = "query,subject,num_HSPs,total_sc,qHSP_nt,qHSP_ovl,qIDN_nt,qHSP_cov,qHSP_idn,sHSP_nt,sHSP_ovl,sIDN_nt,sHSP_cov,sHSP_idn"
print "tab-delimited output: %s" % output_header

query = ""
subject = ""
percent_idn = 0.0
qS = 0
qE = 0
sS = 0 
sE = 0
sc = 0.0
qlen = 0
slen = 0

num_HSPs = 0
total_sc = 0.0

qHSP_coords = list()
qHSP_nt = 0
qHSP_nt_2bAdded = 0
qHSP_ovl = 0
qIDN_nt = 0
qHSP_cov = 0.0
qHSP_idn = 0.0

sHSP_coords = list()
sHSP_nt = 0
sHSP_nt_2bAdded = 0
sHSP_ovl = 0
sIDN_nt = 0
sHSP_cov = 0.0
sHSP_idn = 0.0

for line in args.input:
	tok = line.split('\t')
	if tok[0] != query or tok[1] != subject:
		# print the previous query-species pair, if it is not the first line (don't forget to also print the last line later)
		if first_line == False:
			qHSP_cov = float(qHSP_nt) / qlen
			qHSP_idn = float(qIDN_nt) / qHSP_nt
			sHSP_cov = float(sHSP_nt) / slen
			sHSP_idn = float(sIDN_nt) / sHSP_nt
			if qHSP_cov >= min_qHSP_cov and qHSP_idn >= min_qHSP_idn and \
					sHSP_cov >= min_sHSP_cov and sHSP_idn >= min_sHSP_idn:
				output_line = [query, subject, '%d'%num_HSPs, '%.1f'%total_sc,\
						'%d'%qHSP_nt, '%d'%qHSP_ovl, '%d'%qIDN_nt, '%.2f'%qHSP_cov, '%.2f'%qHSP_idn,\
						'%d'%sHSP_nt, '%d'%sHSP_ovl, '%d'%sIDN_nt, '%.2f'%sHSP_cov, '%.2f'%sHSP_idn ]
				args.output.write('\t'.join(output_line) + '\n')
		else:
			first_line = False
		# initializing
		query = tok[0]
		subject = tok[1]
		qlen = int(tok[12])
		slen = int(tok[13])
		qHSP_coords = [0] * qlen
		sHSP_coords = [0] * slen
		
		num_HSPs = 0
		total_sc = 0.0
		qHSP_nt = 0
		qHSP_ovl = 0
		qIDN_nt = 0
		qHSP_cov = 0.0
		qHSP_idn = 0.0
		sHSP_nt = 0
		sHSP_ovl = 0
		sIDN_nt = 0
		sHSP_cov = 0.0
		sHSP_idn = 0.0
		
	# now process each HSP ...			
	percent_idn = float(tok[2])
	qs = int(tok[6])
	qe = int(tok[7]) 
	ss = int(tok[8])
	se = int(tok[9])
	sc = float(tok[11])
	
	num_HSPs += 1
	total_sc += sc
	
	qHSP_nt_2bAdded = 0
	for i in range (qs - 1, qe):
		if qHSP_coords[i] == 0:
			qHSP_coords[i] = 1
			qHSP_nt_2bAdded += 1
		else:
			qHSP_ovl += 1
	qHSP_nt += qHSP_nt_2bAdded
	qIDN_nt += int( qHSP_nt_2bAdded * percent_idn / 100.0 )

	sHSP_nt_2bAdded = 0
	for i in range ( min(ss, se) - 1, max(ss, se)): # subject HSP coords can be in reverse direction ...
		if sHSP_coords[i] == 0:
			sHSP_coords[i] = 1
			sHSP_nt_2bAdded += 1
		else:
			sHSP_ovl += 1
	sHSP_nt += sHSP_nt_2bAdded
	sIDN_nt += int( sHSP_nt_2bAdded * percent_idn / 100.0 )
	
# process the last HSP:
qHSP_cov = float(qHSP_nt) / qlen
qHSP_idn = float(qIDN_nt) / qHSP_nt
sHSP_cov = float(sHSP_nt) / slen
sHSP_idn = float(sIDN_nt) / sHSP_nt
if qHSP_cov >= min_qHSP_cov and qHSP_idn >= min_qHSP_idn and \
		sHSP_cov >= min_sHSP_cov and sHSP_idn >= min_sHSP_idn:
	output_line = [query, subject, '%d'%num_HSPs, '%.1f'%total_sc,\
			'%d'%qHSP_nt, '%d'%qHSP_ovl, '%d'%qIDN_nt, '%.2f'%qHSP_cov, '%.2f'%qHSP_idn,\
			'%d'%sHSP_nt, '%d'%sHSP_ovl, '%d'%sIDN_nt, '%.2f'%sHSP_cov, '%.2f'%sHSP_idn ]
	args.output.write('\t'.join(output_line) + '\n')
	
args.input.close()
args.output.close()	
print "\ndone\n"