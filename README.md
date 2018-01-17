# CLfinder-OrthNet
## Synopsis
A pipeline to detect co-linearity in gene orders among multiple closely-related genomes (CLfinder) and build networks (OrthNets) connecting orthologs with co-linearity information (co-linear, transposed, or unable to determine) among them as edge properties.
- Ideally and in most cases, each OrthNet consists of orthologs derived from a single ancestral locus, including those duplicated specifically in a subset of genomes.
- In addition to visualize the evolutionary context/history of each ortholog group, users can search ortholog groups (as OrthNets) based on an evolutionary history pattern/context they represent.
- For example, users can retrieve all OrthNets (i.e. ortholog groups) that have underwent:
	* tandem duplication events unique to a genome or a subset of genomes, either mono-, para-, or polyphyletic
	* orthologs deleted, duplicated, or transposed uniquely in a genome or a subset of genomes,
	* orthologs transposed and duplicated uniquely in a genome or a subset of genomes, etc.
- CLfinder and OrthNet are separable modules.  Users can use only the CLfinder module to quickly obtain co-linearity information and a summary matrix of the pairwise comparisons for multiple genomes.  The OrthNet module is optimized for working with CLfinder, but also can accept co-linearity information from other programs.

---
## Before starting
### Prerequisites - programs
- pyhton 2.x - all scripts were written in python 2.7, in linux.
- mcl (https://github.com/JohannesBuchner/mcl) - other clustering options will be added in future updates.
- I assume users are familiar with basic linux commands.
- Optional
	* gffread (https://github.com/gpertea/gffread) - if you start with a .gff file to generate input #1
	* orthoMCL (http://orthomcl.org/common/downloads/software/) - to generate input #2 (see below for an alternative)
	* blast+ (https://www.ncbi.nlm.nih.gov/books/NBK1762/) - to generate input #3
	* Cytoscape (http://cytoscape.org/) - to visualize and print OrthNets

### Prerequisites -genome data
- genome annotation of "representative" gene models (_.gff_, .gff3, or .gtf)
- genome or gene model sequences (.fasta)

### Installing
Copy all python scripts to a folder, add the folder to $PATH, and made them executable:

```
export PATH=<folder>:$PATH # in bash, add this to ~/.bashrc
chmod 755 <folder>/*.py
```

---
## Preparing input files
CLfinder-OrthNet accept three inputs  1. gene model coordinates (genome annotation), 2. within-species paralog groups, and 3. between species "best-hit" pairs for all pair of genomes
### 0. List of genomes
*ProjectID.list* includes all *GenomeIDs* that you want to compare, one per line.  I recommend *GenomeIDs* to be simple (2~5 alphanumeric) and *ProjectID* to be unique by adding date or time-stamp. For example, below, *180101_crucifers* will be the *ProjectID* to compare six crucifer genomes included in the tutorial:
```
echo 'Aly Ath Cru Esa Sir Spa' | tr ' ' '\n' > 180101_Crucifers.list
```

### 1. Input #1: gene model coordinates (genome annotation)
For each genome, coordinates of representative gene models were parsed from genome annotations in *.gtf* format.  The parsed file will have strand, coordinates, number of exons, and length of the mRNA and CDS (ORF), one gene model per line.

1. If genome annotation is in _.gff_ or _.gff3_ format, convert it to _.gtf_:
	```
	gffread input.gff -T -o output.gtf
	```
	If converting multiple files:

	```
	for file in *.gff; do gffread $file -T -o ${file%%.gff}.gtf; done
	```
2. Parse the _.gtf_ file into a _.gtfParsed.txt_ file.  Name the output file as "GenomeID.gtfParsed.txt".  Repeat for all *GenomeIDs*:

	```
	parse_gtf_2table.py -r input.gtf GenomeID.gtfParsed.txt > GenomeID.gtfParsed.log
	```

	#### Important! Genomes should contain one representative gene/transcript model per each locus. See [Note 1](### 1. Obtaining one representative gene model per locus)


### 2. Input #2: within-species paralog groups
A tab-delimited text file with *GeneID* and paralog group ID (*PGID*), one gene per line, for each genome.  Input #2 can be prepared by various methods.  Below are two example options:
1. **method #1** If orthoMCL is available, you can run it for each genome and get "in-paralog" groups. Convert the orthoMCL output (_mclOutput_GenomeID.txt_) to input #2:
	```
	parse_mclOutput.py -rH mclOutput_GenomeID.txt PG > GenomeID.PG
	```
	if converting multiple orthoMCL output files for genomes in *ProjectID.list*:
	```
	while read g; do parse_mclOutput.py mclOutput_${g}.txt PG -rH > ${g}.PG; done < ProjectID.list
	```
	The resulting *GenomeID.PG* can be used as input #2.

	**method #2** Input #2 can be also generate from all-by-all blast results for either representative CDS or deduced amino acid sequences of each genome. For example, if the representative CDS sequences for _GenomeID_ is _GenomeID.cds.rep.fa, first run an all-by-all blastn as follows:

	```
	makeblastdb -in GenomeID.cds.rep.fa -dbtype nucl
	blastn -evalue 1e-5 -max_target_seqs 100000 -outfmt '6 std qlen slen stitle' -db GenomeID.cds.rep.fa -query GenomeID.cds.rep.fa -out out__GenomeID.cds__vs__self.txt
	```
	The blast output can be filtered by `consolidate_blast_HSPs.py` with user-defined thresholds for High-scoring Segment Pair (HSP) coverage on query, subject, or both.  Afterwards, use `create_hard_clusters.py` to cluster paralogs and assign paralog group.  See help screens for both scripts with the _-h_ option for details.

	As long as formatted correctly and generated for each genome, input #2 can be created by any method deemed appropriated by the user to detect and cluster similar sequences.

2. For each genome, add input #2 to input #1 for each genome.  For example, if input #1 and #2 for *GenomeID* is named *GenomeID.gtfParsed.txt* and *GenomeID.PG*, respectively:
	```
	join_files_by_NthCol.py GenomeID.gtfParsed.txt 1 1 GenomeID.PG GenomeID.gtfParsed.PG.txt
	```
	To process multiple genomes listed in *ProjectID.list*:
	```
	while read g; do join_files_by_NthCol.py ${g}.gtfParsed.rep.txt 1 1 ${g}.PG ${g}.gtfParsed.PG.txt; done < ProjectID.list
	```

### 3. Input #3: between-species best-hit pairs
A tab-delimited text file with the GeneID of the query gene and its 'best-hit' or best-hit candidate GeneID in the target genome, one pair per line, for all possible pairs of genomes in *ProjectID.list*.  

1. For all *GenomeIDs* in *ProjectID.list*, create a blast database for the representative CDS sequences in *GenomeID.cds.rep.fa*:
```
makeblastdb -in GenomeID.cds.rep.fa -dbtype nucl
```
To process multiple genomes listed in *ProjectID.list*:
```
while read g; do makeblastdb -in ${g}.cds.rep.fa -dbtype nucl; done < ProjectID.list
```
2. Create blast commands for all possible pair of genomes in *ProjectID.list*.  For blastn:
```
create_pairwiseBLAST_commands.py ProjectID.list -n "-task blastn -evalue 1e-5 -max_target_seqs 10 -outfmt '6 std qlen slen'" > ProjectID_pairwiseBLASTN.sh
```

	Check `create_pairwiseBLAST_commands.py -h` for detailed options to designate folders for CDS sequences or blastn output files, as well as options to use blastp on deduced peptide sequences instead.

	Once blast commands were created, users will want to run it in the background (e.g. using the linux _screen_ command) and multiplex if possible, depending on the computational resource.  Users can add *-num_threads* option to the string given with *-n* option in the example above.

	After running all pairwise blastn, you will have output files named as *out\__GenomeID1\__vs\__GenomeID2.bln.txt*, for all pairs with GenomeID1 != GenomeID2.

	If the user choose to add filters for HSP_cov and/or add HSP_cov and HSP_idn (see `consolidate_blast_HSPs -h`) in the co-linearity information, see [Note 2](### 2. Filtering 'best-hit' pairs based on blast HSP_cov), instead of proceeding to the item 3 below.

3. Convert the blastn output to input #3:
```
for f in out__*__vs__*.bln.txt; do f2=${f##*out__}; cut -f1,2 $f | uniq > BestHits__${f2%%.bln.txt}.list; done
```
	This will generate input #3 for *GenomeID1* and *GenomeID2* as *BestHits\__GenomeID1\__vs\__GenomeID2.list* for all genome pairs.  As long as the file names and formats are correct, input #3 can be created by other methods to detect similar sequences, such as blastp.

---
## Running CLfinder
At this point, the following should be ready:
- List of _GenomeIDs_ in _ProjectID.list_
- Input #1 and #2 for each _GenomeID_, as *GenomeID_gtfParsed.PG.txt*
- Input #3 for all pairs of genomes, as *BestHits\__GenomeID1\__vs\__GenomeID2.list*

Now CLfinder module is ready to run.

1. Add *LocusID* and tandem duplication information to *GenomeID_gtfParsed.PG.txt*:
```
working ... 




---
## Notes
### 1. Obtaining one representative gene model per locus
- To detect co-linearity correctly, CLfinder needs genome coordinates of one gene model per each locus. If possible, select the gene model annotatoin file that inculdes "primary transcript" or "representative gene/isoform".
- `parse_gtf_2table.py -r` reports all gene models in the *.gtf* files whose genomic coordinates are overlapping.  If number of such gene models are small (less than <1%), probably the genome annotation has only representative gene models.
- If genome annotation include isoforms, select only the primary isoforms.  Often isoforms were named as *geneID.1*, *geneID.2*, etc., in which case you can choose only those ending with *.1* in the *.gtfParsed.txt* file.  Check `select_primary_fromGTFparsed.py -h` if you have a list of primary/representative transcript/isoform/gene model IDs and select them in the *.gtfParsed.txt* file.
- If there is no better way to select representative gene/transcript/isoform, `parse_gtf_2table.py -c` will collapse all gene models that have overlapping or identical coordinates, keeping only the longest one.  See the script help for the detail.
- After obtaining the _.gtfParsed.txt_ file, make sure that only gene models included in this file are used for generating Input #2 and #3.  *GeneIDs* should be consistent over all three inputs.

### 2. Filtering 'best-hit' pairs based on blast HSP_cov
- coming soon
