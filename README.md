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

## News
- 2018-06-05 Added options to use MMseqs2 (https://github.com/soedinglab/mmseqs2) to generate inputs #2 and #3.  MMseqs2 is much faster than orthoMCL or BLASTP and also generate one HSP (High-scoring Segment Pair) per each sequence comparison. 

Jump to:
 [**Before starting**](https://github.com/ohdongha/CL_finder#before-starting);
 [**Preparing input files**](https://github.com/ohdongha/CL_finder#preparing-input-files);
 [**Running CLfinder**](https://github.com/ohdongha/CL_finder#running-clfinder);
 [**Running OrthNet**](https://github.com/ohdongha/CL_finder#running-orthnet);
 [**Searching OrthNets**](https://github.com/ohdongha/CL_finder#searching-orthnets);
 [**Annotating OrthNets**](https://github.com/ohdongha/CL_finder#annotating-orthnets-optional);
 [**Notes**](https://github.com/ohdongha/CL_finder#notes);
 [**Tutorial**](https://github.com/ohdongha/CL_finder#tutorial)

---
## Before starting
### Prerequisites - programs
- python 2.7
- mcl (https://github.com/JohannesBuchner/mcl) - other clustering options will be added in future updates.
- MMseqs2 (https://github.com/soedinglab/mmseqs2) - to generate input #2 and #3 (see below for alternatives)	
- I assume users are familiar with basic linux commands.
- Optional
	* gffread (https://github.com/gpertea/gffread) - if you start with a .gff file to generate input #1
	* orthoMCL (http://orthomcl.org/common/downloads/software/) - an alternative to generate input #2
	* blast+ (https://www.ncbi.nlm.nih.gov/books/NBK1762/) - an alternative to generate input #3 and annotation
	* Cytoscape (http://cytoscape.org/) - to visualize and print OrthNets

### Prerequisites -genome data
- genome annotation of "representative" gene models (.gff, .gff3, or .gtf)
- representative gene model sequences for all loci (.fasta)

### Installing
Copy all python scripts to a folder, add the folder to $PATH, and made them executable:

```
export PATH=<folder>:$PATH # in bash, add this to ~/.bashrc
chmod 755 <folder>/*.py
```

---
## Preparing input files
CLfinder-OrthNet accept three inputs: 1. gene model coordinates (genome annotation), 2. within-species paralog groups, and 3. between species "best-hit" pairs for all pair of genomes
### ProjectID and list of genomes
*ProjectID.list* includes all *GenomeIDs* that you want to compare, one per line.  I recommend *GenomeIDs* to be simple (2~5 alphanumeric) and *ProjectID* to be unique by adding date or time-stamp. For example, below, *180101_crucifers* will be the *ProjectID* to compare six crucifer genomes included in the first CLfinder-OrthNet article (https://doi.org/10.1101/236299):
```
echo 'Aly Ath Cru Esa Sir Spa' | tr ' ' '\n' > 180101_Crucifers.list
```

### Input #1: gene model coordinates (genome annotation)
For each genome, coordinates of representative gene models are parsed from genome annotations in *.gtf* format.  The parsed file will have strand, coordinates, number of exons, and length of the mRNA and CDS (ORF), one gene model per line.

1. If genome annotation is in _.gff_ or _.gff3_ format, convert it to _.gtf_:
	```
	gffread input.gff -T -o output.gtf
	```
2. Parse the _.gtf_ file into a _.gtfParsed.txt_ file.  Name the output file as "GenomeID.gtfParsed.txt".  Repeat for all *GenomeIDs*:

	```
	parse_gtf_2table.py -r input.gtf GenomeID.gtfParsed.txt > GenomeID.gtfParsed.log
	```

	#### Important! Genomes should contain one representative gene/transcript model per each locus. See [Note 1](https://github.com/ohdongha/CL_finder#1-obtaining-one-representative-gene-model-per-locus)


### Input #2: within-species paralog groups (PGs)
A tab-delimited text file with *GeneID* and paralog group ID (*PG*), one gene per line for each genome:
```
GeneID	PG
Gene1	PGxxxx
Gene2	PGxxxy
Gene3	PGxxxz
...
```
Input #2 can be prepared by various methods.  Below are two example options:

**method #1** Cluster all representative protein sequences using MMseqs2 for each genome and get "in-paralog" groups. Convert the .tsv output to input #2. Assuming all representative peptide sequences for each genome is _GenomeID.pep.rep.fa_:
```
mkdir tmp_mms # temporary folder for MMseqs2 runs
mmseqs createdb GenomeID.pep.rep GenomeID_DB
mmseqs createindex GenomeID_DB tmp_mms
mmseqs cluster GenomeID_DB GenomeID_c tmp_mms --max-seqs 50000 -c 0.5
mmseqs createtsv GenomeID_DB GenomeID_DB GenomeID_c GenomeID_c.tsv
parse_mmseqs_clusters.py -H GenomeID_c.tsv GenomeID.PG
```
To run MMseqs2 clustering for all genomes in _ProjectID.list_ :
```
mkdir tmp_mms # temporary folder for MMseqs2 runs
while read g; do 
mmseqs createdb ${g}.pep.rep ${g}_DB
mmseqs createindex ${g}_DB tmp_mms
mmseqs cluster ${g}_DB ${g}_c tmp_mms --max-seqs 50000 -c 0.5
mmseqs createtsv ${g}_DB ${g}_DB ${g}_c ${g}_c.tsv
parse_mmseqs_clusters.py -H ${g}_c.tsv ${g}.PG
done < ProjectID.list
```
Proceed to the next step with _GenomeID.PG_ (input #2).

**method #2** (originally the default method) If orthoMCL is available, you can run it for each genome and get "in-paralog" groups. Convert the orthoMCL output (_mclOutput_GenomeID.txt_) to input #2:
```
parse_mclOutput.py -rH mclOutput_GenomeID.txt PG > GenomeID.PG
```
if converting multiple orthoMCL output files for genomes in *ProjectID.list*:
```
while read g; do parse_mclOutput.py mclOutput_${g}.txt PG -rH > ${g}.PG; done < ProjectID.list
```
Proceed to the next step with *GenomeID.PG* (input #2).


### Input #3: between-species best-hit pairs (BHPairs)
A tab-delimited text file with the GeneID of the query gene and its 'best-hit' or best-hit candidate GeneID in the target genome, one pair per line, for all possible pairs of genomes in *ProjectID.list*.  Below, I describe two methods using BLASTN (default) and MMseqs2 (alternative): 

**method #1** For all *GenomeIDs* in *ProjectID.list*, create a blast database for the representative CDS sequences in *GenomeID.cds.rep.fa*:

```
makeblastdb -in GenomeID.cds.rep.fa -dbtype nucl
```
To process multiple genomes listed in *ProjectID.list*:
```
while read g; do makeblastdb -in ${g}.cds.rep.fa -dbtype nucl; done < ProjectID.list
```
Create blast commands for all possible pair of genomes in *ProjectID.list*.  For blastn:
```
create_pairwiseBLAST_commands.py ProjectID.list -n "-task blastn -evalue 1e-5 -max_target_seqs 10 -outfmt '6 std qlen slen'" > ProjectID_pairwiseBLASTN.sh
```
Check `create_pairwiseBLAST_commands.py -h` for detailed options to designate folders for CDS sequences or blastn output files, as well as options to use blastp or MMseqs2 on deduced peptide sequences instead.

Once blast commands (_ProjectID_pairwiseBLASTN.sh_) were created, users will want to run it in the background (e.g. using the linux _screen_ command) and multiplex if possible, depending on the computational resource.  Users can add *-num_threads* option to the string given with *-n* option in the example above. 
	
**method #2** Users can use peptide sequences deduced from representative gene models (_GenomeID.pep.rep.fa_) to generate input #3 using MMseqs2.  First create and index DB for all genomes listed in _ProjectID.list_ (if MMseqs2 was used for the input #2, this step must have been already done):

```
mkdir tmp_mms # temporary folder for MMseqs2 runs
while read g; do 
mmseqs createdb ${g}.pep.rep ${g}_DB
mmseqs createindex ${g}_DB tmp_mms
done < ProjectID.list
```

And create MMseqs2 commands for all pairwise genomes using `create_pairwiseBLAST_commands.py -M`:

```
create_pairwiseBLAST_commands.py ProjectID -M -n "--max-seqs 10" > ProjectID_pairwiseMMseqs2.bash
```
Run the MMseqs2 command (_ProjectID_pairwiseMMseqs2.sh_) as a background process (e.g. using the linux _screen_ command), etc.

**After running all pairwise comparison** you will have output files named as *out\__GenomeID1\__vs\__GenomeID2.txt*, for all pairs with GenomeID1 != GenomeID2.  Convert these blastn (or MMseqs2) results to input #3:
```
for f in out__*.txt; do 
f2=${f##*out__}; cut -f1,2 $f | uniq > BestHits__${f2%%.txt}.list
done
mkdir ./BHPairs; mv BestHits__*.list ./BHPairs
```
This will generate input #3 for all genome pairs as *BestHits\__GenomeID1\_vs\_GenomeID2.list* in the folder _./BHPairs_.  As long as the file names and formats are correct, input #3 can be created by other methods to detect similar sequences, such as blastp.

* Users can add filters to BHPairs at this step, to remove BHPairs with too low HSP_cov, HSP_len, and/or HSP_idn (see `consolidate_blast_HSPs -h`). See [Note 2](https://github.com/ohdongha/CL_finder#2-filtering-best-hit-pairs-based-on-hsp_cov)

---
## Running CLfinder
At this point, the following should be ready:
- List of _GenomeIDs_ in _ProjectID.list_
- Input #1 as _GenomeID.gtfParsed.txt_ for all _GenomeIDs_ in _ProjectID.list_
- Input #2 as _GenomeID.PG_ for all _GenomeIDs_ in _ProjectID.list_
- Input #3 for all pairs of genomes, as *BestHits\__GenomeID1\_vs\_GenomeID2.list* in the folder _./BHPairs_

Now CLfinder module is ready to run:

1. Combine Input #1 (*GenomeID.gtfParsed.txt*) and Input #2 (*GenomeID.PG*) and add tandem duplication information:
	```
	join_files_by_NthCol.py GenomeID.gtfParsed.txt 1 1 GenomeID.PG GenomeID.gtfParsed.PG.txt
	TD_finder.py GenomeID.gtfParsed.PG.txt GenomeID 4 GenomeID.gtfParsed.TD.txt
	```
	To process all genomes listed in *ProjectID.list*:
	```
	while read g; do 
	join_files_by_NthCol.py ${g}.gtfParsed.rep.txt 1 1 ${g}.PG ${g}.gtfParsed.PG.txt
	TD_finder.py ${g}.gtfParsed.PG.txt $g 4 ${g}.gtfParsed.TD.txt
	done < ProjectID.list
	```
	This will identify all paralogs in the same paralog group and within 4 loci (min_TD_loci) as tandem duplicated (_td_), and report the number of _td_ events and genes in those events as _ProjectID_TD.log_. Users can modify the *min_TD_loci* parameter as needed. _GenomeID.gtfParsed.TD.txt_ files now include columns for _td_ events as well as numerical _LocusIDs_ for each genome.

2. Run CLfinder:
	```
	CL_finder_multi.py ProjectID -nr -T .gtfParsed.TD.txt -b ./BHPairs -W 20 -N 3 -G 20 -o ./ProjectID_out
	```
	This will run CLfinder for each _GenomeID_, comparing it to all other genomes, and report (with _-r_ option) numbers of co-linear, lineage-specific, and transposed best-hit pairs, shared _td_ events, etc. Users can modify the _window_size_ (_-W_), _num_loci_trshld_ (_-N_), and _gap_loci_trshld_ (_-G_) as needed.  With _-n_ option, it also create an output _ProjectID.4OrthNet.input_ which can be used for updating best-hit pairs (next section) or building OrhNets.  See `CL_finder_multi.py -h` for details on options and parameters.

3. Update best-hit pairs:
	```
	update_BestHitPairs.py ProjectID ./ProjectID_out/ProjectID.4OrthNet.input -b BHPairs -o BHPairs.1
	```
	CLfinder can search for an alternative best-hit among many best-hit candidates with less similarity scores, if such alternative best-hit can achieve a reciprocal best-hit pairs.  This process will update best-hit pairs to prefer reciprocal relationship.
	```
	CL_finder_multi.py ProjectID -unrp -T .gtfParsed.TD.txt -b BHPairs.1 -o ProjectID_out.1 -W 20 -N 3 -G 20
	```
	Re-run CLfinder based on the updated best-hit pairs (_-u_ option).  With _-p_ option, CLfinder also print reciprocal best-hit pairs and co-linearity between them for all pairs of genomes, which can be useful for determining synonymous (_Ks_) or four degenerated site (_4d_) substitution rates between pairs ([Note 3](https://github.com/ohdongha/CL_finder#3-calculating-substitution-rates-between-best-hit-pairs-with-codeml)).
4. Creating a summary report for all pairwise CLfinder analyses:
 	```
	create_CLfm_summary.py ProjectID CLfinder_summary.txt -p ProjectID_out.1
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
	update_OrthNet_after_mcl.py ProjectID mcl/ProjectID_TD1.5_rC1.2_rNC1.0_uC0.3_uNC0.25_I1.2_mclOut.PC.txt -b BHPairs.1 -o1 BHPairs.2 -o2 ProjectID_out.2 -u
	format_OrthNetEdges_4SIF.py ./ProjectID_out.2/ProjectID.clstrd.afterMCL.edges
	```
	This process again searches for alternative best-hit pairs to maximize pairing within each OrthNet after the _mcl_ clustering. Then, reformat the _.edges_ file to a _.sif_ file to finalize OrthNets:

	The resulting _ProjectID.clstrd.afterMCL.edges.sif_ file includes all OrthNets identified by the OrthNet module in _.sif_ format for Cytoscape. However, I advise not trying to open the entire OrthNets at once in Cytoscape, since the file is expected to be quite large. See the next section and find how to search and extract subsets of OrthNets using search by GeneID, OrthNetID, or evolutionary context search.

4. Run CLfinder one last time to reflect updated best-hit pairs and OrthNet information:
	```
	CL_finder_multi.py ProjectID -unr -b BHPairs.2 -o ProjectID_out.2 -W 20 -N 3 -G 20
	```
	Users may want to re-create the CLfinder summary report at this point:
 	```
	create_CLfm_summary.py ProjectID CLfinder_summary.afterOrthNet.txt -p ProjectID_out.2
	```
---
## Searching OrthNets
OrthNets are stored as .sif (simple interaction file) format, i.e., tab-delimited text with NodeID1, CLtype, NodeID2, and OrthNetID. NodeIDs are formatted as 'GenomeID|GeneID'. Users can use simple linux `grep` commands to retrieve OrthNets using OrthNetIDs. 

For example, to find out OrthNetID of the OrthNet to which a gene of interest belongs:
```
grep 'GeneID' ProjectID.clstrd.afterMCL.edges.sif
```
To retrieve an OrthNet using OrthNetID, or group of OrthNets with OrthNetIDs listed in _OrthNetID.list_, one per line:
```
grep -P "\tOrthNetID$" ProjectID.clstrd.afterMCL.edges.sif > ProjectID.OrthNetID.sif
while read oid ; do 
grep -P "\t${oid}$" ProjectID.clstrd.afterMCL.edges.sif > ProjectID.${oid}.sif
done < OrthNetID.list
```
Finally, users can search OrthNets for nodes representing a specific evolutionary context (see examples in the Synopsis).  For this, `search_OrthNet.py` uses CLfinder results files (_CL files_) created by the item 4 of the previous section and regular expression patterns created by users as the query.  See `search_OrthNet -h` for details. 
	
---
## Annotating OrthNets (optional)

This step combines the tabulated results of a CLfinder-OrthNet run, to a single table.  This table includes, in addition to CLfinder results, the following information for each gene locus:
- a functional annotation (i.e. the best blast-hit from a model species, as well as the blast HSP coverage and identity),
- the OrthNet ID and a one-line summary of the properties of all edges connecting the gene locus,   
- the median (md) and standard deviation (sd) of CDS (ORF) lengths for the gene and all of its best-hits.  This is to determine whether the gene locus contains a truncated or complete ORF.

See the [Supplementry Dataset S1](https://figshare.com/articles/DatasetS1_xlsx/5825937) of the first CLfinder-OrthNet manuscript for an example.

Separately, this step generate a summary including the following information for each OrthNet:
- number of nodes in each species, counting either all or nodes with complete ORFs only,
- top three most prevalent functional annotations among OrthNet nodes,
- the average (av), md, and sd of ORF sizes of either all nodes or nodes with complete ORFs.


1. Prepare the anotation file:

	This step requires a BLASTN or BLASTP output generated by local blast+ with ```-outfmt '6 std qlen slen stitle'``` option, comparing either CDS or protein sequences for all loci as the query, to sequences from a functionally-annotated reference model species (e.g., for flowering plants, users can obtain <i>Arabidopsis thaliana</i> reference sequences from these Cyberse links: [CDS](http://datacommons.cyverse.org/browse/iplant/home/araport/public_data/Araport11_Release_201606/annotation/Araport11_genes.201606.cds.repr.fasta.gz) / [protein](http://datacommons.cyverse.org/browse/iplant/home/araport/public_data/Araport11_Release_201606/annotation/Araport11_genes.201606.pep.repr.fasta.gz)).  
	
	If GeneIDs are not unique in target genomes, we recommend query sequence names to be formatted as 'SpeciesID|GeneID' for the BLAST.

	Assuming <i>Annotation_blast_output.txt</i> is the blast output:

	```
	consolidate_blast_HSPs.py -Hs Annotation_blast_output.txt temp
	cut -f 1,8,9,15 temp > Annotation_consolidated.txt
	cut -f 1,15 temp > Annotation_consolidated.short.txt
	```
	
	These files are to be used later for annotating both gene loci and OrthNets.
	
2. Combine CLfinder results and add annotation:
	First, combine all CLfinder results in the output folder of the last ```CL_finder_multi.py``` run. The '-O' option adds the information for ORF lengths for each locus and its best-hits.  See ```combine_CLfm_results.py -h``` for more details. 
	```
	combine_CLfm_results.py ProjectID temp_output_combined -p ./ProjectID_out.2 -F .20.3.20.afterMCL.txt -O
	join_files_by_NthCol.py temp_output_combined 2 1 Annotation_consolidated.txt ProjectID.combined.annotated.txt e
	```
	The <i>ProjectID.combined.annotated.txt</i> is the final output. 
 
3. Generate a summary of OrthNet node copy numbers and annotation:

	The ```update_OrthNet_after_mcl.py``` script (see [Running OrthNet](https://github.com/ohdongha/CL_finder#running-orthnet) item #3) generates an <i>.mclOutput</i> output file, which can be used for creating an annotated summary of OrthNets. See ```parse_mclOutput.py -h``` for more details.
	```
	parse_mclOutput.py -Hsrx -a Annotation_consolidated.short.txt -p _internal_ ./ProjectID_out.2/ProjectID.clstrd.afterMCL.nodes.mclOutput ON
	```
	This will create <i>ProjectID.all.mclOutput.summary.txt</i>, which contains counts of all nodes for each species, as well as top three representative annotations (i.e. annotations appearing most frequently among the members of the OrthNet), for each OrthNet.
	
4. Count the number of nodes with complete ORFs:

	The ```compare_OrthNet_ORFsize.py``` script report the number of nodes with complete ORFs, as well as the av, md, and sd of ORF sizes of all nodes in each OrthNet.  It uses the combined CLfinder output file generated in item #2 in this section as the input.

	```
	compare_OrthNet_ORFsize.py ProjectID ProjectID.combined.annotated.txt -m 0
	```
	The '-m' option defines the minimum % of the ORF length of a gene compared to the median of its best-hits.  '-m 0' will count all nodes and report to <i>ProjectID.all.nodeCounts</i>, while '-m 70' only counts nodes with ORF size at least 70% of the median of its best-hits and report to <i>ProjectID.cORF70.nodeCounts</i>.  See ```compare_OrthNet_ORFsize.py -h``` for more details.
	
---
## Notes
### 1. Obtaining one representative gene model per locus
- To detect co-linearity correctly, CLfinder needs genome coordinates of one gene model per each locus. If possible, select the gene model annotatoin file that inculdes "primary transcript" or "representative gene/isoform".
- `parse_gtf_2table.py -r` reports all gene models in the *.gtf* files whose genomic coordinates are overlapping.  If number of such gene models are small (less than <1%), probably the genome annotation has only representative gene models.
- If genome annotation include isoforms, select only the primary isoforms.  Often isoforms were named as *geneID.1*, *geneID.2*, etc., in which case you can choose only those ending with *.1* in the *.gtfParsed.txt* file.  Check `select_primary_fromGTFparsed.py -h` if you have a list of primary/representative transcript/isoform/gene model IDs and select them in the *.gtfParsed.txt* file.
- If there is no better way to select representative gene/transcript/isoform, `parse_gtf_2table.py -c` will collapse all gene models that have overlapping or identical coordinates, keeping only the longest one.  See the script help for the detail.  One can also modify the *.gtf* file, using [gtf2gtf.py in CGAT](https://www.cgat.org/downloads/public/cgat/documentation/scripts/gtf2gtf.html), to condense isoforms in each locus (thanks Keeley Adams for finding out this tool).
- After obtaining the _.gtfParsed.txt_ file, make sure that only gene models included in this file are used for generating Input #2 and #3.  *GeneIDs* should be consistent over all three inputs.

### 2. Filtering 'best-hit' pairs based on HSP_cov
- coming soon

### 3. Calculating substitution rates between best-hit pairs with codeml
- coming soon

---
## Tutorials
- coming soon
