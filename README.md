# CLfinder-OrthNet
A pipeline to detect co-linearity in gene orders among multiple closely-related genomes (CLfinder) and build networks (OrthNets) connecting orthologs with co-linearity information (co-linear, transposed, or unable to determine) among them as edge properties.  
- Ideally and in most cases, each OrthNet consists of orthologs derived from a single ancestral locus, including those duplicated specifically in a subset of genomes.
- In addition to visualize the evolutionary context/history of each ortholog group, users can search ortholog groups (as OrthNets) based on an evolutionary history pattern/context they represent.
- For example, users can retrieve all OrthNets (i.e. ortholog groups) that have underwent: 
-- tandem duplication events unique to a genome or a subset of genomes, either mono-, para-, or polyphyletic
-- orthologs transposed uniquely in a genome or a subset of genomes, 
-- orthologs transposed and duplicated uniquely in a genome or a subset of genomes, etc.

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
- genome annotation of "representative" gene models (.gtf) - see below on how to get "representative" gene models
- genome or gene model sequences (.fasta)

### Installing
Copy all python scripts to a folder, add the folder to $PATH, and made them executable:

`export PATH=<folder>:$PATH # in bash, add this to ~/.bashrc`

`chmod 755 <folder>/*.py`


## Preparing input files
CLfinder-OrthNet accept three inputs  1. gene model coordinates (genome annotation), 2. within-species paralog groups, and 3. between species "best-hit" pairs for all pair of genomes
### 0. List of genomes
*ProjectID.list* includes all *GenomeIDs* that you want to compare, one *GenomeID* per line.  I recommend *GenomeIDs* to be simple (2~5 alphanumeric) and *ProjectID* to be unique by adding the date. For example, *180101_crucifers* will be the *ProjectID* to compare six crucifer genomes included in the tutorial:

`echo 'Aly Ath Cru Esa Sir Spa' | tr ' ' '\n' > 180101_Crucifers.list`


### 1. Input #1 Gene model coordinates (genome annotation)
1. If genome annotation was given as a _.gff_ or _.gff3_ file, convert it to a _.gtf_ file:

	`gffread input.gff -T -o output.gtf`

	If converting multiple files:

	`for file in *.gff; do gffread $file -T -o ${file%%.gff}.gtf; done`

2. Convert the _.gtf_ file into a _.gtfParsed.txt_ file.  Name the output file as "GenomeID.gtfParsed.txt".  Repeat for all *GenomeIDs*:

	`parse_gtf_2table.py -r input.gtf GenomeID.gtfParsed.txt > GenomeID.gtfParsed.log`

#### Important: one representative gene model per locus
- To detect co-linearity correctly, CLfinder needs genome coordinates of one gene model per each locus. If possible, select the gene model annotatoin file that inculdes "primary transcript" or "representative gene/isoform".
- `parse_gtf_2table.py -r` reports all gene models in the *.gtf* files whose genomic coordinates are overlapping.  If number of such gene models are small (less than <1%), probably the genome annotation has only representative gene models.
- If genome annotation include isoforms, select only the primary isoforms.  Often isoforms were named as *geneID.1*, *geneID.2*, etc., in which case you can choose only those ending with *.1* in the *.gtfParsed.txt* file.  Check `select_primary_fromGTFparsed.py -h` if you have a list of primary/representative transcript/isoform/gene model IDs and select them in the *.gtfParsed.txt* file.
- If there is no better way to select representative gene/transcript/isoform, `parse_gtf_2table.py -c` will collapse all gene models that have overlapping or identical coordinates, keeping only the longest one.  See the script help for the detail. 
- After obtaining the _.gtfParsed.txt_ file, make sure that only gene models included in this file are used for generating Input #2 and #3.  *GeneIDs* should be consistent over all three inputs.

### 2. Input #2 Within-species paralog groups
1. If OrthoMCL is available, you can run it for each genome and get "in-paralog" groups.  

... to be continued soon, 