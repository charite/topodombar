/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package annotation;

import genomicregions.Gene;
import genomicregions.GenomicElement;
import genomicregions.GenomicSet;
import java.util.ArrayList;

/**
 * This class implements an effect mechanism of a CNV/SV that involves the loss 
 * or gain of regulatory interaction between enhancer elements and genes that are
 * associated the phenotype of the patient. Furthermore, deletion or duplication
 * of topological domain boundaries might be involved in this effect mechanism.
 * 
 * @author jonas
 */
public class InteractionChange extends EffectAnnotation {
    
    /**
     * List of enhancers that are involved in the CNV effect mechanism
     */
    private GenomicSet<GenomicElement> enhancers;
    
    /**
     * List of topological domain boundary elements involved in the effect mechanism
     */
    private GenomicSet<GenomicElement> boundaries;
    
    /**
     * Interaction change between {@link Gene} and enhancers. 
     * This String sould be "gain" or "loss".
     */
    private String interactionChange;
    
    /**
     * Constructor with all needed arguments.
     * 
     * @param method
     * @param gene
     * @param score
     * @param enhancers
     * @param boundaries
     * @param interactionChange 
     */
    public InteractionChange(String method, Gene gene, Double score, 
            GenomicSet<GenomicElement> enhancers,
            GenomicSet<GenomicElement> boundaries, 
            String interactionChange
            ){

        super(method, gene, score);        
        this.enhancers = enhancers;
        this.boundaries = boundaries;
        this.interactionChange = interactionChange;
        
    }

    /**
     * List of enhancers that are involved in the CNV effect mechanism
     * @return the enhancers
     */
    public GenomicSet<GenomicElement> getEnhancers() {
        return enhancers;
    }

    /**
     * List of topological domain boundary elements involved in the effect mechanism
     * @return the boundaries
     */
    public GenomicSet<GenomicElement> getBoundaries() {
        return boundaries;
    }

    /**
     * Interaction change between {@link Gene} and enhancers.
     * This String sould be "gain" or "loss".
     * @return the interactionChange
     */
    public String getInteractionChange() {
        return interactionChange;
    }
    
    
}
