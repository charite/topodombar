topodombar
==========

Phenotypic analysis of structural variations (SVs) affecting topologically associating domains (TADs) and enhancers.
This tool implements and extends the  analysis performed in 
Ibn-Salem J et al. (2014) Deletions of Chromosomal Regulatory Boundaries are Associated with Congenital Disease, *Genome Biology*.
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0423-1

## Requirements
 * Java 8
 * Mavon >= 3.3.9

## Module structure
This project is structured into three modules:

 * topodombar.core
 * topodombar.commandline
 * topodombar.gui
 
 Most of the functionallity and classes are implemented in the `topodombar.core` module.

## Build commandline tool from source
```
mvn package
```
This should create the following .jar file:
`topodombar.commandline/target/topodombar.commandline-0.0.1-SNAPSHOT-jar-with-dependencies.jar`

##  Example usage
```
java -jar topodombar.commandline/target/topodombar.commandline-0.0.1-SNAPSHOT-jar-with-dependencies.jar  \
  -i topodombar.core/src/test/resources/example_CNV.tab \
  -d topodombar.core/src/test/resources/example_domains.tab \
  -g topodombar.core/src/test/resources/example_genes.tab \
  -e topodombar.core/src/test/resources/example_enhancer.tab \
  -O topodombar.core/src/test/resources/example_ontology.obo \
  -a topodombar.core/src/test/resources/example_genes_to_penotype.txt \
  -o example_out.tab
```
This creates several output file with the prefix `example_out.tab`. For more information on input arguments use `-h` option.


```
java -jar topodombar.commandline/target/topodombar.commandline-0.0.1-SNAPSHOT-jar-with-dependencies.jar -h
```

## Input data
The following input data is needed for the analysis (see below for information on file formats):

 * CNV file
 * Domain file
 * Gene file
 * Enhancer file
 * HPO file
 * HPO term to gene file

### CNV file
Tab-separated file with CNVs. One CNV per row with the following columns:

1. Chromosome
2. CNV start coordinate
3. CNV end coordinate
4. Unique identifier for CNV/patient
5. Type of CNV (deletion/insertion)
6. Phenotype annotation of the patient as list of HPO terms. Term IDs separated by ';'
7. Phenotype category (target term) as single HPO term ID.

### Domain file
BED file with non-overlapping topological domains with the following columns:

1. Chromosome
2. Start
3. End
4. Unique identifier

### Boundary file
BED file with topological domain boundaries with the following columns:

1. Chromosome
2. Start
3. End
4. Unique identifier

### Gene file
Tab-separated file with one gene per row and the following columns:

1. Chromosome
2. Start
3. End
4. Strand
5. EntrezGene ID
6. Associated HPO term IDs separated by ';'

### Enhancer files
BED file with enhancers and following columns:

1. Chromosome
2. Start
3. End
4. Unique identifier

### HPO file
The human phenotype ontology as .obo file

### HPO term to gene file
Mapping for each HPO term to its associated genes.
Tab separated file with HPO term ID in first and EntrezGene ID in second column.
Only one term and gene per line.


## Development setup

This section describes the setup of java, NetBeans, maven and git on an Ubuntu 14.04 system.

### Java
 - Install Oracle Java 8 JDK by downloadding the latest file from here http://www.oracle.com/technetwork/java/javase/downloads/
The file name should be '''jdk-8uVERSION-linux-x64.tar.gz'''
 - Unpack file and copy it with root permissions to  **/opt/Oracle_Java/**
 - Close all running webbrowsers and ensure that no other java browser plugins are active.
 - Use the following commands to setup the alternative system and configure it:

```
sudo update-alternatives --install "/usr/bin/java" "java" "/opt/Oracle_Java/jdk1.8.0_VERSION/bin/java" 1
sudo update-alternatives --install "/usr/bin/javac" "javac" "/opt/Oracle_Java/jdk1.8.0_VERSION/bin/javac" 1
sudo update-alternatives --install "/usr/bin/javaws" "javaws" "/opt/Oracle_Java/jdk1.8.0_VERSION/bin/javaws" 1
sudo update-alternatives --install "/usr/bin/jar" "jar" "/opt/Oracle_Java/jdk1.8.0_VERSION/bin/jar" 1 

sudo update-alternatives --set "java" "/opt/Oracle_Java/jdk1.8.0_VERSION/bin/java"
sudo update-alternatives --set "javac" "/opt/Oracle_Java/jdk1.8.0_VERSION/bin/javac"
sudo update-alternatives --set "javaws" "/opt/Oracle_Java/jdk1.8.0_VERSION/bin/javaws"
sudo update-alternatives --set "jar" "/opt/Oracle_Java/jdk1.8.0_VERSION/bin/jar" 
```

### Netbeans
 - Download the installation file from the NetBeans website: https://netbeans.org/downloads/
 - Run it with root permissions and put as installation path **/opt/netbeans-VERSION**

### Maven
Install maven with the following command:
```
sudo apt-get install maven
```

### Git
Install git with the following command:
```
sudo apt-get install git
```

