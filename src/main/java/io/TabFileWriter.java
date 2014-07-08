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
import java.util.List;

/**
 *
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class TabFileWriter<T extends GenomicElement> {
    
    private File file;
    private Path path;
    
    // define charset for nio.Files.write function
    private Charset charset = Charset.forName("utf-8");

    /**
     * 
     * @param outFile 
     */
    public TabFileWriter(File outFile){
        this.file = outFile;
        this.path = outFile.toPath();

    }
    

    /**
     * 
     * @param path 
     */
    public TabFileWriter(String path){
        this.path = Paths.get(path);
        this.file = new File(path);
    }
    
    /**
     * 
     * @param elements
     * @throws IOException 
     */
    public void write(GenomicSet<T> elements) throws IOException{
        
        // convert elements to list of output lines
        List lines = new ArrayList<String>();

        
        //T tmp = (T) new GenomicElement("chrTmp", 0, 1, "tmp");
        
        /*
        Construct header line. Note, that the function getOutputHeaderLine is a
        non-static memberfunction and therfore any of the proper type (subcalss
        of GenomicElements) is needed to get the header string.
        */
        String headerLine;
        if (elements.isEmpty()){
            // if element set is empty, take generic header of a sample GenomicElement
            headerLine = new GenomicElement("chrTmp", 0, 1, "tmp").getOutputHeaderLine();
        }else{
            T anyElement = elements.get(elements.keySet().toArray()[0]);
            headerLine = anyElement.getOutputHeaderLine();
        }
        
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
    
}
