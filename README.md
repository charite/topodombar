topodombar
==========

Phenotypic analysis of microdeletions and topological chromosome domain boundaries. These scripts are meant to document the analysis performed in 

> Ibn-Salem J et al.<br>
> Deletions of Chromosomal Regulatory Boundaries are Associated with Congenital Disease.<br>
> [Genome Biology 2014 15:423](http://genomebiology.com/2014/15/9/423) 


## Requirements

### Python and Java
- Python 2.7 
- Python package 'numpy' 
- Java '1.7.0 55'

### Input data
The following input data is needed for the analysis (see below for information on file formats):
* CNV file
* Domain file
* Boundary file
* Gene file
* HPO term to gene file
* Enhancer file
* HPO file
* Target term file

#### CNV file
Tab-separated file with CNVs. One CNV per row with the following columns:

1. Chromosome
2. CNV start coordinate
3. CNV end coordinate
4. Unique identifier for CNV/patient
5. Type of CNV (deletion/insertion)
6. Phenotype annotation of the patient -- a list of HPO terms. Term IDs separated by ';'
7. Genes within the CNV as list of Entrez Gene IDs separated by ';'
8. Genes upstream of the CNV between CNV breakpoint and end of underlying topological domain as list of Entrez Gene ID separated by ';'
9. Genes downstream of the CNV between CNV breakpoint and end of underlying topological domain as list of Entrez Gene ID separated by ';'
10. Genes upstream of the CNV within a distance window of 400kb as list of Entrez Gene ID separated by ';'
11. Genes downstream of the CNV within a distance window of 400kb as list of Entrez Gene ID separated by ';'
12. Phenotype category (target term) as single HPO term ID.

#### Domain file
BED file with non-overlapping topological domains with the following columns:
1. Chromosome
2. Start
3. End
4. Unique identifier

#### Boundary file
BED file with topological domain boundaries with the following columns:
1. Chromosome
2. Start
3. End
4. Unique identifier

#### Gene file
Tab-separated file with one gene per row and the following columns:
1. Chromosome
2. Start
3. End
4. EntrezGene ID
5. Associated HPO terms separated by ';'

#### HPO term to gene file
Mapping for each HPO term to its associated genes.
Tab separated file with HPO term ID in first and EntrezGene ID in second column.
Only one term and gene per line.

#### Enhancer files
For each tissue a BED file with tissue specific enhancers. 
The files should be named <tissueName>.tab and should have a unique ID in the fourth column.

#### HPO file
The human phenotype ontology as .obo file


#### Target term file
Tab-separated file with all target terms and corresponding tissue names analysed.
Each line corresponds to one target term and should have the three columns:
1. HPO term ID
2. HPO term name
3. Tissue name


## Usage
The analysis consists of three main steps.

### 1. Computation of PhenoGram score.

```
java -Xmx2G -jar bin/CnvStatistics.jar \
	-d <CNV file> -c 6,7,8,9,10 -k 1 -l 0 \
	-u <HPO file> \
	-a ALL_SOURCES_TYPICAL_FEATURES_genes_to_phenotype.txt \
	-o <CNV file>.hpo_phenoScore
```

The file ALL_SOURCES_TYPICAL_FEATURES_genes_to_phenotype.txt and the HPO file can be downloaded from the Human Phenotype Ontology project page: http://www.human-phenotype-ontology.org/

### 2. Get the maximal PhenoGram score per gene in each region

```
python phenogram_score.py -i <CNV file>.hpo_phenoScore -f max \
    -o <CNV file>.hpo_phenoScore.max_scores
```

### 3. Analyse CNVs for topological boundary disruption (TDBD)

```
python barrier_analysis.py \
	-c <CNV file>.hpo_phenoScore.max_scores \
	-d <Domain file> \
	-b <Boundary file> \
	-g <Gene fil> \
	-hg <HPO term to gene file> \
	-e <enhancer directory>  \
	-hpo <HPO file> \
	-p <Target term file> \
	-o <CNV file>.hpo_phenoScore.max_scores.barrier
	
```


Usage information for the Python scripts can be seen by executing the script with '-h' option.

```
python barrier_analysis.py -h
usage: barrier_analysis.py [-h] -c CNV_FILE -d DOMAIN_FILE -b BOUNDARY_FILE -g
                           GENES_FILE -hg TERM_TO_GENE_FILE [-e ENHANCER_DIR]
                           [-ef ENHANCER_FILE] -hpo HPO_FILE -p
                           TARGET_PHENOTYPE_FILE
                           [-of {complete,any,percent50}] [-w WINDOW_SIZE]
                           [-bs BIN_SIZE] -o OUTPUT_FILE [-sf]

optional arguments:
  -h, --help            show this help message and exit
  -c CNV_FILE, --cnv_file CNV_FILE
                        input CNV file in tab separated format. With columns
                        chr, start, end, id
  -d DOMAIN_FILE, --domain_file DOMAIN_FILE
                        Domain file in .bed format
  -b BOUNDARY_FILE, --boundary_file BOUNDARY_FILE
                        Domainboundary file in .bed format
  -g GENES_FILE, --genes_file GENES_FILE
                        Genes file in .tab format
  -hg TERM_TO_GENE_FILE, --term_to_gene_file TERM_TO_GENE_FILE
                        Tab separated file, that maps each phenotype term
                        (including decendants) to its associated genes
  -e ENHANCER_DIR, --enhancer_dir ENHANCER_DIR
                        path to directory with enhancer data matching the
                        target tissues. Assume files with <tissue>.bed.id
  -ef ENHANCER_FILE, --enhancer_file ENHANCER_FILE
                        path to a single file with enhancers.
  -hpo HPO_FILE, --hpo_file HPO_FILE
                        Human Phenotype Ontology file in .obo format
  -p TARGET_PHENOTYPE_FILE, --target_phenotype_file TARGET_PHENOTYPE_FILE
                        Tab separated file with target phenotypes as HPO ID in
                        first column
  -of {complete,any,percent50}, --overlap_function {complete,any,percent50}
                        Function to compute the overlap of boundaries.
                        'complete' requires the CNV to completely overlap a
                        boundary, 'any' requires only a partial overlap and
                        'percent50' requires at least 50 percent of the
                        boundary affect.
  -w WINDOW_SIZE, --window_size WINDOW_SIZE
                        Window size for testing enhancer adaption mechanism
                        without the boundary disruption effect. That is search
                        for matching enahncer and gene in a fixed distance
                        window in the flanking regions of the deletion.
  -bs BIN_SIZE, --bin_size BIN_SIZE
                        bin size for faster acces to regions while computiong
                        region overlaps. Default is 10^6
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        output file
  -sf, --sparse_output_format
                        Write sparse output format, that is only CNVs that
                        match the target_phenotype and only some (see soruce
                        code) columns.
```
	
