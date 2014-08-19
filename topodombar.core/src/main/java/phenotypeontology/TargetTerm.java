/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package phenotypeontology;

import genomicregions.GenomicElement;
import genomicregions.GenomicSet;
import io.TabFileParser;
import java.io.IOException;
import ontologizer.go.Term;

/**
 * Target terms represent general phenotype groups or tissues for which specific
 * genomic data, such as enhancers, are available.
 * 
 * @author jonas
 */
public class TargetTerm{
    
    /**Term object representing this target term. */
    private Term term;
    
    /**Enhancer (or other regulatory elements) associated with this target term. */
    private GenomicSet<GenomicElement> enhancers;
    
    /**A name or representing this target term. */
    private String name;
    
    /**
     * Consturctor for a {@link TargetTerm} instance by reading enhancer data
     * form the given path. 
     * 
     * @param term the term object form the corresponding phenotype ontology
     * @param name name of the target term or tissue
     * @param enhancerFilePath  path to enhancer data specific to the target term
     * @throws IOException if the enhancer file cannot be read.
     */
    public TargetTerm(Term term, String name, String enhancerFilePath) throws IOException{
        
        this.term = term;
        this.name = name;
        
        // parse the enhancer data
        TabFileParser enhancerParser = new TabFileParser(enhancerFilePath);
        this.enhancers = enhancerParser.parse();

    }

    /**
     * Consturctor for a {@link TargetTerm} instance.
     * @param term the term object form the corresponding phenotype ontology
     * @param name name of the target term or tissue
     * @param enhancers set of enhancers associated with this target term.
     */
    public TargetTerm(Term term, String name, GenomicSet<GenomicElement> enhancers) {
        this.term = term;
        this.name = name;
        this.enhancers = enhancers;
    }
    
    /**
     * Returns the target term ID as a string.
     * @return 
     */
    public String getIDAsString(){
        return this.term.getIDAsString();
    }
    
    /**
     * Term object representing this target term.
     * @return the term
     */
    public Term getTerm() {
        return term;
    }

    /**
     * Enhancer (or other regulatory elements) associated with this target term.
     * @return the enhancers
     */
    public GenomicSet<GenomicElement> getEnhancers() {
        return enhancers;
    }

    /**
     * A name or representing this target term.
     * @return the name
     */
    public String getName() {
        return name;
    }
    
}
