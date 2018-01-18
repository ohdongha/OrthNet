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

Jump to:
 [**Before starting**](https://github.com/ohdongha/CL_finder#before-starting);
 [**Preparing input files**](https://github.com/ohdongha/CL_finder#preparing-input-files);
 [**Running CLfinder**](https://github.com/ohdongha/CL_finder#running-clfinder);
 [**Running OrthNet**](https://github.com/ohdongha/CL_finder#running-orthnet);
 [**Searching OrthNets**](https://github.com/ohdongha/CL_finder#searching-orthnets);
 [**Notes**](https://github.com/ohdongha/CL_finder#notes);
 [**Tutorial**](https://github.com/ohdongha/CL_finder#tutorial)

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
### ProjectID and list of genomes
*ProjectID.list* includes all *GenomeIDs* that you want to compare, one per line.  I recommend *GenomeIDs* to be simple (2~5 alphanumeric) and *ProjectID* to be unique by adding date or time-stamp. For example, below, *180101_crucifers* will be the *ProjectID* to compare six crucifer genomes included in the tutorial:
```
echo 'Aly Ath Cru Esa Sir Spa' | tr ' ' '\n' > 180101_Crucifers.list
```

### Input #1: gene model coordinates (genome annotation)
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

	#### Important! Genomes should contain one representative gene/transcript model per each locus. See [Note 1](https://github.com/ohdongha/CL_finder#1-obtaining-one-representative-gene-model-per-locus)


### Input #2: within-species paralog groups
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

### Input #3: between-species best-hit pairs
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

	If the user choose to add filters for HSP_cov and/or add HSP_cov and HSP_idn (see `consolidate_blast_HSPs -h`) in the co-linearity information, see [Note 2](https://github.com/ohdongha/CL_finder#2-filtering-best-hit-pairs-based-on-blast-hsp_cov), instead of proceeding to the item 3 below.
3. Convert the blastn output to input #3:
	```
	for f in out__*__vs__*.bln.txt; do f2=${f##*out__}; cut -f1,2 $f | uniq > BestHits__${f2%%.bln.txt}.list; done
	mkdir ./BHPairs.bln; mv BestHits__*.list ./BHPairs.bln
	```
	This will generate input #3 for *GenomeID1* and *GenomeID2* as *BestHits\__GenomeID1\__vs\__GenomeID2.list* for all genome pairs in the folder _./BHPairs.bln_ (the scripts work for any folder name; .bln stands for blastn).  As long as the file names and formats are correct, input #3 can be created by other methods to detect similar sequences, such as blastp.

---
## Running CLfinder
At this point, the following should be ready:
- List of _GenomeIDs_ in _ProjectID.list_
- Combined input #1 and #2 for each _GenomeID_, as *GenomeID_gtfParsed.PG.txt*
- Input #3 for all pairs of genomes, as *BestHits\__GenomeID1\__vs\__GenomeID2.list* all located in the folder _./BHPairs.bln_.

Now CLfinder module is ready to run:

1. Add *LocusID* and tandem duplication information to *GenomeID_gtfParsed.PG.txt*:
	```
	while read g; do TD_finder.py ${g}.gtfParsed.PG.txt $g 4 ${g}.gtfParsed.TD.txt; done < ProjectID.list > ProjectID_TD.log 2>&1
	```
	This will identify all paralogs in the same paralog group and within 4 loci (min_TD_loci) as tandem duplicated (_td_), and report the number of _td_ events and genes in those events as _ProjectID_TD.log_. Users can modify the *min_TD_loci* parameter as needed. _GenomeID.gtfParsed.TD.txt_ files now include columns for _td_ events as well as numerical _LocusIDs_ for each genome.

2. Run CLfinder:
	```
	CL_finder_multi.py ProjectID -nr -T .gtfParsed.TD.txt -b ./BHPairs.bln -W 20 -N 3 -G 20 -o ./ProjectID_bln
	```
	This will run CLfinder for each _GenomeID_, comparing it to all other genomes, and report (with _-r_ option) numbers of co-linear, lineage-specific, and transposed best-hit pairs, shared _td_ events, etc. Users can modify the _window_size_ (_-W_), _num_loci_trshld_ (_-N_), and _gap_loci_trshld_ (_-G_) as needed.  With _-n_ option, it also create an output _ProjectID.4OrthNet.input_ which can be used for updating best-hit pairs (next section) or building OrhNets.  See `CL_finder_multi.py -h` for details on options and parameters.

3. Update best-hit pairs:
	```
	update_BestHitPairs.py ProjectID ./ProjectID_bln/ProjectID.4OrthNet.input -b BHPairs.bln -o BHPairs.bln.1
	```
	CLfinder can search for an alternative best-hit among many best-hit candidates with less similarity scores, if such alternative best-hit can achieve a reciprocal best-hit pairs.  This process will update best-hit pairs to prefer reciprocal relationship.
	```
	CL_finder_multi.py ProjectID -unrp -T .gtfParsed.TD.txt -b BHPairs.bln.1 -o ProjectID_bln.1 -W 20 -N 3 -G 20
	```
	Re-run CLfinder based on the updated best-hit pairs (_-u_ option).  With _-p_ option, CLfinder also print reciprocal best-hit pairs and co-linearity between them for all pairs of genomes, which can be useful for determining synonymous (_Ks_) or four degenerated site (_4d_) substitution rates between pairs ([Note 3](https://github.com/ohdongha/CL_finder#3-calculating-substitution-rates-between-best-hit-pairs-with-codeml)).
4. Creating a summary report for all pairwise CLfinder analyses:
 	```
	create_CLfm_summary.py ProjectID CLfinder_summary.txt -p ProjectID_bln.1
	```
	This script looks into a CLfinder output folder and create a summary matrix for all query-target genome pairs, reporting number of co-linear (_cl_), lineage-specific (_ls_), transposed (_tr_), and not determined (_nd_) due to too fragmented genome scaffold assembly.  See `create_CLfm_summary.py -h` for details.
---
## Running OrthNet
The OrthNet module accept a tab-delimited text file with two genes, i.e., best-hit pairs from different genomes or tandem duplicated paralogs from the same genome, and their co-linearity relationship. The CLfinder model generate such a file (_ProjectID.4OrthNet.input_), as described above.  As long as formatted correctly, the OrthNet modeule can accept co-linearity information from other pipeline or sources.  See `create_OrthNet.py -h` for the required input file format.

1. Create initial hard clusters:
	```
	create_OrthNet.py ProjectID -sd -o ./
	```
	This step connect all best-hit pairs or tandem duplicated paralogs to create initial clusters.

	With _-d_ option, tandem duplicated paralogs will be included in OrthNets.  Note that the script expects _GenomeID.gtfParsed.TD.txt_ files in the ./ folder when run with _-d_ option. With _-s_ option, it also identifies clusters that need to be further separated into sub-clusters. See `create_CLfm_summary.py -h` for details on options and parameters.

2. Markov clustering (mcl) of hard clusters:
	```
	mkdir ./mcl
	mcl_OrthNet.py ProjectID -o mcl -sc -w weights4mcl.list -I 1.2
	```
	Since this process may create a large number of temporary files, I suggest to create a separate working folder (_./mcl_).  The _weights4mcl.list_ includes edge weights users can assign to each type of edges. Users can also select the inflation rate for mcl with _-I_ option. See `mcl_OrthNet.py -h` for details.

	The output filename includes edge weight information.  If default values were used, the output file will be *ProjectID_TD1.5_rC1.2_rNC0.5_uC0.6_uNC0.25_I2.0_mclOut.PC.txt*.

3. Update best-hit pairs and OrthNets after mcl:
	```
	update_OrthNet_after_mcl.py ProjectID mcl/ProjectID_TD1.5_rC1.2_rNC1.0_uC0.3_uNC0.25_I1.2_mclOut.PC.txt -b BHPairs.bln.1 -o1 BHPairs.bln.2 -o2 170316_C.2  -u
	```


---
## Searching OrthNets

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

### 3. Calculating substitution rates between best-hit pairs with codeml
- coming soon

---
## Tutorials
