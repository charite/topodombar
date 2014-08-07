/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package io;

import genomicregions.CNV;
import genomicregions.GenomicElement;
import genomicregions.GenomicSet;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.apache.commons.lang3.StringUtils;
import phenotypeontology.PhenotypeData;

/**
 * Writes to tab separated output file.
 * 
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 * @param <T> Type parameter can be {@link GenomicElement} or a subclass of it like {@link CNV}.
 */
public class TabFileWriter<T extends GenomicElement> {
    
    /** Path to the output file */
    private final Path path;
    
    // define charset for nio.Files.write function
    private final Charset charset = Charset.forName("utf-8");

    
    /**
     * Construct a {@link TabFileWriter} object.
     * 
     * @param outFile output file
     */
    public TabFileWriter(File outFile){
        this.path = outFile.toPath();

    }
    

    /**
     * Construct a {@link TabFileWriter} object.
     * 
     * @param path path to output file
     */
    public TabFileWriter(String path){
        this.path = Paths.get(path);
    }
    
    /**
     * 
     * @param elements
     * @throws IOException 
     */
    public void write(GenomicSet<T> elements) throws IOException{
        
        // convert elements to list of output lines
        List lines = new ArrayList<String>();

        
        String headerLine = T.getOutputHeaderLine();
        /*
        Construct header line. Note, that the function getOutputHeaderLine is a
        static memberfunction that is however specidfic for the type T (GenomicElement or one
        of its subclasses)       
        */
//        if (elements.isEmpty()){
//            // if element set is empty, take generic header of a sample GenomicElement
//            headerLine = new GenomicElement("chrTmp", 0, 1, "tmp").getOutputHeaderLine();
//        }else{
//            T anyElement = elements.get(elements.keySet().toArray()[0]);
//            headerLine = anyElement.getOutputHeaderLine();
//        }
        
        // add header to output lines:
        lines.add(headerLine);
        
        // iterate over each element and add a line for it to the output lines
        for ( T e : elements.values() ){
            
            // call the memberfunction toOutputLine to convert each element to 
            // one output line in the appropriate format. 
            lines.add(e.toOutputLine());
        
        }
                
        // wirte output lines to output file
        java.nio.file.Files.write(path, lines, charset);
    
    }
    
    /**
     * Writes CNVs to a tab-separated output file including a header line.
     * 
     * @param cnvs CNVs to be written to the output file
     * @param phenotypeData
     * 
     * @throws IOException 
     */
    public void writeCNVs(GenomicSet<CNV> cnvs, PhenotypeData phenotypeData) throws IOException{
        
        // Create a header line as tab-separated string with all column identifiers
        final String headerLineCNV = StringUtils.join(new String[]{
                    "#chr",
                    "start",
                    "end",
                    "name",
                    "type", 
                    "phenotypes", 
                    "targetTerm",
                    "boundaryOverlap",
                    "overlapGenes",
                    "overlapScore",
                    "leftAdjacentGenes",
                    "leftAdjacentScore",
                    "rightAdjacentGenes",
                    "rightAdjacentScore",
                    "leftAdjacentEnhancers",
                    "rightAdjacentEnhancers",
                    "EffectMechanismTDBD",
                    "EffectMechanismEA",
                    "EffectMechanismEAlowG"
                }, '\t');
        
        // get all output lines from the GenomicSet object:
//        ArrayList<String> outLines = cnvs.getOutputLines();
        ArrayList<String> outLines = new ArrayList<String>();

        // sort CNVs by there effect mechanism class
        ArrayList<CNV> sortedCNVs = new ArrayList<CNV>(cnvs.values());
        Collections.sort( sortedCNVs, CNV.EFFECTMECHANISM_TDBD_ORDER);

        for (CNV c : sortedCNVs){
            outLines.add(c.getOutputLineWithRelevantGenes(phenotypeData));
        }
        // add the header line to the first position
        outLines.add(0, headerLineCNV);
        
        // write all lines to the output file.
        java.nio.file.Files.write(path, outLines, charset);
    }

}
