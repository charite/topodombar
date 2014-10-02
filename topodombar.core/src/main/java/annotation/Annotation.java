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
import genomicregions.GenomicElement;
import java.util.ArrayList;

/**
 * This class implements an effect mechanism annotation of a structural variant, 
 * such as deletions, duplication or inversion with respect to affected genes, 
 * enhancers and regulatory interaction between them. 
 * A structural variation (SV) or CNV might hold a list of serveral annotation 
 * object. Each annotation object is related to a gene that is affected by the SV.
 * 
 * @author ibnsalem
 */
public class Annotation {
    
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
     * List of enhancers that are involved in the CNV effect mechanism
     */
    private ArrayList<GenomicElement> enhancers;
    
    /**
     * List of topological domain boundary elements involved in the effect mechanism
     */
    private ArrayList<GenomicElement> boundaries;
    
    /**
     * Interaction change between {@link Gene} and enhancers. 
     * This String sould be "gain" or "loss".
     */
    private String interactionChange;
    
    /**
     * An indicator {@link String} for the method and or settings used to create 
     * this annotation.
     */
    private String method;
    
}
    

