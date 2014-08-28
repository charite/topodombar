/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package annotation;

import genomicregions.Gene;
import genomicregions.GenomicSet;
import java.util.HashMap;

/**
 * Provides functionallity to annotate {@GenomicSet} of {@Gene}s.
 * 
 * @author jonas
 */
public class AnnotateGenes {
    
    /**
     * Annotates genes with GeneSymbol identifier.
     * 
     * @param genes Genes to be annotated
     * @param idToSymbol mapping of EntrezGenes IDs to GeneSymbols
     */
    public static void addGeneSymbol(GenomicSet<Gene> genes, HashMap<String, String> idToSymbol){
       
        for (Gene g: genes.values()){
            
            // check if GeneSymbol is available for this gene
            if (idToSymbol.containsKey(g.getName())){
  
                // set gene symbol
                g.setSymbol(idToSymbol.get(g.getName()));
            
            }else{
                // if no symbol was found set a dot as default.
                g.setSymbol(".");
            }
        }
    }
}
