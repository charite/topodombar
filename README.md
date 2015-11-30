topodombar
==========

Phenotypic analysis of microdeletions and topological chromosome domain boundaries. 
This tool implements the  analysis performed in 
Ibn-Salem J et al. (2014) Deletions of Chromosomal Regulatory Boundaries are Associated with Congenital Disease, *Genome Biology*.
http://www.ncbi.nlm.nih.gov/pubmed/25315429



## Requirements

### Java
 - TODO

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

## Development setup

This section describes the setup of java, NetBeans, maven and git on an Ubuntu 14.04 system.

### Java
 - Install Oracle Java 8 JDK by downloadding the latest file from here http://www.oracle.com/technetwork/java/javase/downloads/
The file name should be '''jdk-8uVERSION-linux-x64.tar.gz'''
 - Unpack file and copy it with root permissions to  **/opt/Oracle_Java/**
 - Close all running webbrowsers and ensure that no other java browser plugins are active.
 - Use the following commands to setup the alternative system and configure it:

'''
sudo update-alternatives --install "/usr/bin/java" "java" "/opt/Oracle_Java/jdk1.8.0_VERSION/bin/java" 1
sudo update-alternatives --install "/usr/bin/javac" "javac" "/opt/Oracle_Java/jdk1.8.0_VERSION/bin/javac" 1
sudo update-alternatives --install "/usr/bin/javaws" "javaws" "/opt/Oracle_Java/jdk1.8.0_VERSION/bin/javaws" 1
sudo update-alternatives --install "/usr/bin/jar" "jar" "/opt/Oracle_Java/jdk1.8.0_VERSION/bin/jar" 1 

sudo update-alternatives --set "java" "/opt/Oracle_Java/jdk1.8.0_VERSION/bin/java"
sudo update-alternatives --set "javac" "/opt/Oracle_Java/jdk1.8.0_VERSION/bin/javac"
sudo update-alternatives --set "javaws" "/opt/Oracle_Java/jdk1.8.0_VERSION/bin/javaws"
sudo update-alternatives --set "jar" "/opt/Oracle_Java/jdk1.8.0_VERSION/bin/jar" 
'''

## Netbeans
 - Download the installation file from the NetBeans website: https://netbeans.org/downloads/
 - Run it with root permissions and put as installation path **/opt/netbeans-VERSION**

## Maven
Install maven with the following command:
''' 
sudo apt-get install maven
'''

## Git
Install git with the following command:
''' 
sudo apt-get install git
'''

