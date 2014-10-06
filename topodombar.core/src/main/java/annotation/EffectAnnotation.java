/*
 * Copyright (C) 2014 ibnsalem
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package annotation;

import genomicregions.Gene;

/**
 * This class implements an effect mechanism annotation of a structural variant, 
 * such as deletions, duplication or inversion with respect to affected genes.
 * A structural variation (SV) or CNV might hold a list of several annotation 
 * objects. Each annotation object is related to a gene that is affected by the SV.
 * More specific {@link EffectAnnotation} such as {@link InteractionChange} 
 * or {@link GeneDosage} can be implemented as subclasses to this class.
 * 
 * @author ibnsalem
 */
public class EffectAnnotation {
    
    /**
     * The {@link Gene} that is affected by the CNV
     */
    private Gene gene;
    
    /** 
     * A score to quantify this effect mechanism. This score should be the 
     * Phenomatch score that quantifies the similarity of phenotypes associated
     * with the gene and phenotypes observed in the patient.
     */
    private Double score;
    
    /**
     * An indicator {@link String} for the method and or settings used to create 
     * this annotation.
     */
    private String method;
    
    /**
     * Constructor with all needed input parameters
     * 
     * @param method
     * @param gene
     * @param score 
     */
    public EffectAnnotation(String method, Gene gene, Double score){
        
        this.method = method;
        this.gene = gene;
        this.score = score;
        
    }

    /**
     * The {@link Gene} that is affected by the CNV
     * @return the gene
     */
    public Gene getGene() {
        return gene;
    }

    /**
     * A score to quantify this effect mechanism. This score should be the
     * Phenomatch score that quantifies the similarity of phenotypes associated
     * with the gene and phenotypes observed in the patient.
     * @return the score
     */
    public Double getScore() {
        return score;
    }

    /**
     * An indicator {@link String} for the method and or settings used to create
     * this annotation.
     * @return the method
     */
    public String getMethod() {
        return method;
    }
}
    

