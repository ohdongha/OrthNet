#!/usr/bin/env python

synopsis = "\n\
TD_finder.py <input.PG.txt> <species_code> <N> <output.txt>\n\n\
 - to a tab-delimited <input.PG.txt>, add locusID, find Tandem Duplication (TD)\n\
     events and assign TDid, and print to <output.txt>.\n\
 - requires a sorted protein-coding (pc) gene coordinate file, parsed from\n\
     a genome annotation file (.gtf) and an Ortholog Group (PG) file, parsed\n\
     from e.g. an OrthoMCL output run within the species. \n\n\
 1. Input file(s)\n\
 - <input.PG.txt> is the final output of the following processes:\n\
     parse_gtf_2table.py <genome.gtf> <genome_parsed.txt>\n\
     parse_OrthoMCL.py <geneID.list> <mclOutput> <mclOutput_parsed.PG>\n\
     join_files_by_NthCol.py <genome_parsed.txt> 1 1 <mclOutput_parsed.PG> <input.PG.txt>\n\
 - for detail, refer to the synopsis of each script.\n\
 - make sure <input.PG.txt> accords with the following:\n\
     sorted based on coordinates\n\
     includes only one representative pc gene model per locus\n\n\
 2. Finding tandem duplication (TD)\n\
 - all pc genes are given numerical locusID based on their positions in the genome.\n\
 - if a locus has gene(s) with the same PGid (defined by 'parse_OrthoMCL.py'), within the adjacent <N> loci,\n\
     then the locus is tandem duplicated (TD), and given a TDid. \n\
 - TDid contains the <species_code>, PGid, 'TD', and a unique number (e.g. 'Spc|PG_xxxxx__TDyy')\n\
 by ohdongha@gmail.com 20170317 ver 0.2\n"

#version_history
#20171225 ver 0.2.2 nomenclature changed: OG (Ortholog Group) -> PG (Paralog Group)
#20170510 ver 0.2.1 bug fixed (header)
#20170317 ver 0.2 synopsis modified
#20160422 ver 0.1 TDid format changed
#20151230 ver 0.0
 
#### Example: lines in a <input.PG.txt> file
##geneID  Chr     Str     mRNA_start      mRNA_end        #exon_mRNA      mRNA_len        CDS_start       CDS_end #exon_CDS       CDS_len PG
##Sp1g00010       ch1-1   -       1159    1863    1       705     1159    1863    1       705     PG_02995
##Sp1g00020       ch1-1   +       2166    5101    6       1062    2166    5101    6       1062    na
##Sp1g00030       ch1-1   -       5839    7378    6       795     5839    7378    6       795     PG_01696
##Sp1g00040       ch1-1   -       8513    11989   7       2706    8513    11989   7       2706    PG_00120
##Sp1g00050.b     ch1-1   -       16149   20251   10      2795    17002   20250   10      1941    PG_00598
##Sp1g00060       ch1-1   +       24041   24919   2       660     24041   24919   2       660     PG_03939
##Sp1g00070.a.x   ch1-1   -       24316   27023   5       2048    25153   26922   5       1110    PG_07541
##Sp1g00080       ch1-1   +       26686   26880   1       195     26686   26880   1       195     na
##Sp1g00090       ch1-1   -       27550   28270   2       342     27550   28270   2       342     PG_18626
##Sp1g00100       ch1-1   -       28633   30257   2       1551    28633   30257   2       1551    PG_03940
##Sp1g00110       ch1-1   -       30562   34373   20      1464    30562   34373   20      1464    PG_03941
##Sp1g00120       ch1-1   +       34636   35286   2       423     34636   35286   2       423     PG_03942
##Sp1g00130       ch1-1   +       36939   37603   4       420     36939   37603   4       420     na
##Sp1g00140       ch1-1   -       38602   39220   3       447     38602   39220   3       447     PG_03943

locusID_gap_between_chromosome = 200

import sys
	
try: 
	fin_PG = open(sys.argv[1], "rU")
	species_code = sys.argv[2]
	max_TD_loci_distance = int(sys.argv[3])
	fout = open(sys.argv[4], "w")
#	fout_test = open(sys.argv[2] + "_test.151230_1117.txt","w")
except (ValueError, IndexError) :
	print synopsis
	sys.exit(0)

## reading in <input.PG.txt> and assign locusID
locusID = 0
number_loci = 0
previous_chr = ""
current_chr = ""
current_PG = ""
locus_line_dict = dict()
locus_PG_dict = dict()
PG_loci_dict = dict()
first_line = True
header = ""

print "reading %s as the <input.PG.txt>:" % sys.argv[1]

for line in fin_PG:
	newline_accepted = 0
	tok = line.split('\t')
	try:
		if first_line: # 170510 bug fix
			header = line.strip()
			first_line = False # 170510 bug fix
		elif tok[11][0:2] == 'na' or tok[11][0:3] == 'PG_' :
			previous_chr = current_chr
			current_chr = tok[1].strip()
			current_PG = tok[11].strip()
			if current_chr == previous_chr:
				locusID = locusID + 1
			elif first_line == 1:
				locusID = 1
				first_line = 0
			else:
				locusID = locusID + locusID_gap_between_chromosome
			locus_line_dict[locusID] = line
			locus_PG_dict[locusID] = current_PG
			if current_PG != 'na':
				if current_PG in PG_loci_dict:
					PG_loci_dict[current_PG].append(locusID)
				else:
					PG_loci_dict[current_PG] = [locusID]	
			number_loci = number_loci + 1
	except (ValueError, IndexError) :
		print "There is a non-valid line."
fin_PG.close()
		
## assigning TDid, iterating over PGids,
#TDid = 0
#number_TD = 0
#current_TD_PG = ""
#previous_TD_PG = ""
#TDid_dict = dict()
#nearby_PGs = set()

previous_locus = 0
current_locus = 0
TD_continue = 0
number_TD_in_PG = 0
number_TD_total = 0
number_gene_in_TD = 0
TDid = ""
TDid_dict = dict()

for key in sorted(PG_loci_dict):
#	fout_test.write(key +' : ' + ','.join(map(str,PG_loci_dict[key])) + '\n')
	previous_locus = PG_loci_dict[key][0]
	number_TD_in_PG = 0
	for n in range(1, len(PG_loci_dict[key])):
		current_locus = PG_loci_dict[key][n]
		if (current_locus - previous_locus) <= max_TD_loci_distance : 
			if n==1 or TD_continue == 0 :
				number_TD_in_PG = number_TD_in_PG + 1
#				TDid = key + '_' + species_code + '_TD' + '%02d' % number_TD_in_PG
				TDid = species_code + '|' + key + '__TD' + '%02d' % number_TD_in_PG
				number_TD_total = number_TD_total + 1
				TD_continue = 1
			TDid_dict[current_locus] = TDid
			number_gene_in_TD = number_gene_in_TD + 1
			if previous_locus not in TDid_dict:
				TDid_dict[previous_locus] = TDid
				number_gene_in_TD = number_gene_in_TD + 1
		else:
			TD_continue = 0
			TDid = ""
		previous_locus = current_locus

print "Total %d genes were found in %d TD events." % (number_gene_in_TD, number_TD_total)			
#for key in sorted(locus_line_dict):
#	if locus_PG_dict[key] == 'na':
#		TDid_dict[key] = ""
#	else:
#		nearby_PGs.clear()
#		try:
#			for n in range(int(key) - max_TD_loci_distance, int(key) + max_TD_loci_distance):
#				if n != key and n in locus_PG_dict:
#					nearby_PGs.add(locus_PG_dict[n])
#			if locus_PG_dict[key] in nearby_PGs:
#				previous_TD_PG = current_TD_PG
#				current_TD_PG = locus_PG_dict[key]
#				number_TD = number_TD + 1
#				if current_TD_PG != previous_TD_PG:
#					TDid = TDid + 1
#				TDid_dict[key] = species_code + '_TD' + '%04d' % TDid
#		except (ValueError, IndexError, KeyError) :
#			print "Something not nice happened.  Try to troubleshoot..."
#
## output

print "printing to %s \n" % sys.argv[4]
PG_header = header.split('\t')[-1].strip() # 170510 bug fix
fout.write(header.replace(PG_header, species_code + "_PG") + '\t' + species_code + "_locusID" + '\t' + species_code + "_TDid" + '\n')

for key in sorted(locus_line_dict):
	fout.write(locus_line_dict[key].strip() + '\t' + str(key) + '\t' + TDid_dict.get(key, "-") + '\n')
fout.close()
