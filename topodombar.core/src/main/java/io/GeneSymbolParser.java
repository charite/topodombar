/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;

/**
 * This parser reads the annotation files provided by the HPO project to create
 * a mapping of EntrezGeneIDs to GeneSymbols.
 * The annotaion file for the Human Penotype Ontology (HPO) is named: 
 * "ALL_SOURCES_TYPICAL_FEATURES_genes_to_phenotype.txt"
 * For the cross species Uberpheno ontology it is:
 * "HSgenes_crossSpeciesPhenoAnnotation.txt"
 * Both files can be optained form the HPO download page: 
 * <a href="http://www.human-phenotype-ontology.org/contao/index.php/downloads.html">http://www.human-phenotype-ontology.org/contao/index.php/downloads.html</a>
 * 
 * @author jonas
 */
public class GeneSymbolParser {
    
    /** {@link Path} to the file that is read by this parser */
    private final Path path;

    public GeneSymbolParser(String annotationFile){
       this.path = Paths.get(annotationFile);        
    }
    
    public HashMap<String, String> parseEntrezToSymbol() throws IOException{
        
        HashMap<String, String> entrezToSymbol = new HashMap<String, String>();
        String entrezID = null;
        String symbol = null;
        
        try {

            // read the first line of the file to descide wether it it is the HPO
            // annotation fiel or the Uberpheno annotaion file
            BufferedReader in = new BufferedReader(new FileReader(this.path.toString()));

            String firstLine = in.readLine();
            boolean isHpo = false;
            boolean isUpheno = false;
            if (firstLine.startsWith("#Entrez Gene ID of human gene ; Gene symbol ; Annotated Uberpheno")) {
                    isUpheno = true;
            } else if (firstLine.startsWith("#Format: entrez-gene-id<tab>entrez-gene-symbol<tab>HPO-Term-Name<tab>HPO-Ter")) {
                    isHpo = true;
            } else {
                    throw new IOException("Can't handle annotation-file format!");
            }

            // read all lines of the file
            for ( String line : Files.readAllLines( this.path, StandardCharsets.UTF_8 ) ){
                
                // ignore comment lines
                if (line.startsWith("#"))
                    continue;
                
                // in case of HPO annotation format split line by TAB chars 
                // and fine entrez ID in first column and gene symbol in second.
                if (isHpo){
                    // split line by TABs
                    String [] cols = line.split("\t");
                    entrezID = cols[0];
                    symbol = cols[1];
                }
                
                // in case of Uberpheno annotaion format
                if (isUpheno){
                    String [] cols = line.split(";");
                    entrezID = cols[0];
                    symbol = cols[1];                    
                }
                
                // add IDs to mapping
                entrezToSymbol.put(entrezID, symbol);
            }
            
        }catch (IOException e){
            System.err.println("[Error] while parsing file:" + this.path);
            throw e;
        }
        
        return entrezToSymbol;
    }
}
