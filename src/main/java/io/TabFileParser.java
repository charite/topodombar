/*
 * Copyright (c) 2014, Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

package io;

import genomicregions.CNV;
import genomicregions.Gene;
import genomicregions.GenomicElement;
import genomicregions.GenomicElementSet;
import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Reads a tab-separated file with genomic elements like CNVs or enhancers 
 * into an GenomicElementSet object.
 * 
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class TabFileParser {
    
    private File file;
    private Path path;
    
    /**
     * Construct a {@link TabFileParser} object form an input {@link File} object
     * 
     * @param inputFile Input {@link File} object 
     */
    public TabFileParser(File inputFile){
        this.file = inputFile;
        this.path = inputFile.toPath();
    }
    
    /**
     * Construct a {@code TabFileParser} object from an input path.
     * @param path  path to input file.
     */
    public TabFileParser(String path){
        this.path = Paths.get(path);
        this.file = new File(path);
    }
    
    /**
     * Reads a TAB separated file with genomic elements like enhancers or topological domains.
     * Assumes each line in the file to represent a genomic element. The first
     * four columns should contain the following: chromosome, start, end, and name.
     * Note, the genomic coordinates are assumed in 0-based half-open format 
     * like in the BED format specifications 
     * (See {@link http://genome.ucsc.edu/FAQ/FAQformat.html#format1}).
     * 
     * @return {@link GenomicRegionSet} with all elements from the input file.
     * 
     * @throws IOException if file can not be read. 
     */
    public GenomicElementSet parse() throws IOException{
        
        // construct new set of genomic reigons:
        GenomicElementSet ges = new GenomicElementSet();
        
        for ( String line : Files.readAllLines( path, StandardCharsets.UTF_8 ) ){
            
            // split line by TABs
            String [] cols = line.split("\t");
            
            
            String chr = cols[0];
            int start = Integer.parseInt(cols[1]);
            int end = Integer.parseInt(cols[2]);
            String name = cols[3];

            // create new element and add it to the set
            ges.put(name, new GenomicElement(chr, start, end, name));
        }
        
        return ges;
    }
    
    /**
     * Reads a TAB separated file with CNVs.
     * Assumes each line in the file to represent a CNV. The first
     * four columns should contain the following: chromosome, start, end, and name.
     * Three additional columns are are assumed. the fifth column should contain the CNV type.
     * In the sixth column we assume the 
     * phenotypes annotation of the patient in HPO terms separated by semicolon ';'.
     * The seventh column should contain the so called target term, a single more general 
     * HPO term that represent the more abstract phenotypical class ore tissue 
     * type to which the patient belongs to.
     * Note, the genomic coordinates are assumed in 0-based half-open format 
     * like in the BED format specifications 
     * (See {@link http://genome.ucsc.edu/FAQ/FAQformat.html#format1}).
     * 
     * @return {@link GenomicRegionSet} with {@link CNV} objects from the input file.
     * 
     * @throws IOException if file can not be read. 
     */
    public GenomicElementSet<CNV> parseCNV() throws IOException{
        
        // construct new set of genomic reigons:
        GenomicElementSet<CNV> cnvs = new GenomicElementSet();
        
        for ( String line : Files.readAllLines( path, StandardCharsets.UTF_8 ) ){
            
            // split line by TABs
            String [] cols = line.split("\t");
            
            // parse columns
            String chr = cols[0];
            int start = Integer.parseInt(cols[1]);
            int end = Integer.parseInt(cols[2]);
            String name = cols[3];
            
            // parse CNV specific columns
            String type = cols[4];
            ArrayList<String> phenotypes = new ArrayList(Arrays.asList( cols[5].split(";")));
            String targetTerm = cols[6];
            
            // create new {@link CNV} object
            CNV cnv = new CNV(chr, start, end, name, type, phenotypes, targetTerm);
            
            // add it to the set
            cnvs.put(name, cnv);
        }
        
        return cnvs;
    }
    
    /**
     * Reads a TAB separated file with genes.
     * Assumes the .tab format form the barrier project.
     * Note, that the file should contain genes with unique Entrez Gene IDs in the
     * fourth column.
     * 
     * @return  a set of genes as {@link GenomicElmentSet} of {@link Gene} objects.
     * @throws IOException if the file can not be read
     */
    public GenomicElementSet<Gene> parseGene() throws IOException{
        
        // construct new set of genomic reigons:
        GenomicElementSet<Gene> genes = new GenomicElementSet();
        
        for ( String line : Files.readAllLines( path, StandardCharsets.UTF_8 ) ){
            
            // split line by TABs
            String [] cols = line.split("\t");
            
            // parse columns
            String chr = cols[0];
            int start = Integer.parseInt(cols[1]);
            int end = Integer.parseInt(cols[2]);
            String name = cols[3];

            // create new {@link Gene} object
            Gene g = new Gene(chr, start, end, name);

            // parse Gene specific columns
            g.strand = cols[4];
            g.phenotpyes = Arrays.asList( cols[5].split(";"));
            //Gene Symbol is not contained in the .tab format of the barrier project 
            g.symbol = "."; 
            
            // add it to the set
            genes.put(name, g);
        }
        
        return genes;
    }
}
