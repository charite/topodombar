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

package genomicregions;

import java.util.ArrayList;
import java.util.List;

/**
 * Implements a copy number variation (CNV) object.
 * Extends the {@link GenomicElement} class and has additional annotations like
 * type, HPO terms, and target term.
 *
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class CNV extends GenomicElement {
    
    /**
     * Type of CNV (loss or gain).
     */
    public String type;
    /**
     * List of phenotypes as HPO term IDs.
     */
    public List<String> phenotpyes;
    /**
     * Target term or phenotype category  as single general HPO term ID.
     */
    public String targetTerm;
    
    // annotations:
    
    /**
     * List of overlapping boundaries.
     */
    public List<GenomicElement> boundaryOverlap;
    
    /**
     * true if CNV overlaps a boundary element.
     */
    public boolean hasBoundaryOverlap;
    
    /**
     * List of overlapping genes (any overlap)
     */
    public List<Gene> geneOverlap;
    
    /**
     * Constructor for CNV object.
     * Construct a {@link GenomicElement} and sets all {@link CNV} specific 
     * annotations to default values.
     * 
     * @param chr   Chromosome identifier
     * @param start Start coordinate (0-based)
     * @param end   End coordinate (0-based, half-open)
     * @param name  Name or ID of the CNV
     * 
     * @throws IllegalArgumentException 
     */
    public CNV(String chr, int start, int end, String name) throws IllegalArgumentException {
        super(chr, start, end, name);
        
        // set default values for annotations
        type = ".";
        phenotpyes = new ArrayList();
        targetTerm = ".";
    }
    
    /**
     * Constructor for {@link CNV} object.
     * Construct an CNV object by taking all annotations as arguments.
     * 
     * @param chr   Chromosome identifier
     * @param start Start coordinate (0-based)
     * @param end   End coordinate (0-based, half-open)
     * @param name  Name or ID of the CNV
     * @param type  CNV type (loss or gain)
     * @param phenotypes    List of HPO term IDs that represent the phenotypes used to annotate the patient carrying the CNV.
     * @param targetTerm A unspecific target term as a HPO term ID that is used to group patients in cohort or associate patients to tissues for which enhancer data his available.
     */
    public CNV(String chr, int start, int end, String name, String type, List<String> phenotypes, String targetTerm){

        // consturct an CVN object using the constructor of the {@link GenomicElement} super calss
        super(chr, start, end, name);
        
        // add annotations
        this.type = type;
        this.phenotpyes = phenotypes;
        this.targetTerm = targetTerm;        
    }
    
    
}
