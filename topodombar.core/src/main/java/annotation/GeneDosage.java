/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package annotation;

import genomicregions.Gene;

/**
 * Implements a gene dosage effect mechanism by eigther "loss" or "gain" of the 
 * affected gene.
 * 
 * @author jonas
 */
public class GeneDosage extends EffectAnnotation {
    
    /**
     * String indicating if there is a "loss" or "gain" in the dosage of the affected gene.
     */
    private String dosageChange;
    
    /**
     * Cunstructor for a {@link GeneDosage} object.
     * 
     * @param method
     * @param gene
     * @param score
     * @param dosageChange 
     */
    public GeneDosage(String method, Gene gene, Double score,
            String dosageChange){
        
        super(method, gene, score);
        this.dosageChange = dosageChange;
                
    }

    /**
     * String indicating if there is a "loss" or "gain" in the dosage of the affected gene.
     * @return the dosageChange
     */
    public String getDosageChange() {
        return dosageChange;
    }

}
