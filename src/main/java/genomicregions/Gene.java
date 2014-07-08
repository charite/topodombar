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
 * Implements a gene object.
 * Extends the {@link GenomicElement} class and has additional annotations like
 * type, HPO terms, and target term.
 *
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class Gene extends GenomicElement {
    
    /**
     * List of phenotypes as HPO term IDs that are associated with this gene.
     */
    public List<String> phenotpyes;
    /**
     * Target term or phenotype category  as single general HPO term ID.
     */
    
    /**
     * Coding strand ("+" or "-").
     */
    public String strand;
    
    /**
     * The GeneSymbol identifier for the gene.
     * Note, the {@code name} field should hold the Entrez Gene ID.
     */
    public String symbol;
    
    
    /**
     * Constructor for {@link Gene} object.
     * Construct a {@link GenomicElement} and sets all {@link Gene} specific 
     * annotations to default values.
     * 
     * @param chr   Chromosome identifier
     * @param start Start coordinate (0-based)
     * @param end   End coordinate (0-based, half-open)
     * @param name  Name or ID of the CNV
     * 
     * @throws IllegalArgumentException 
     */
    public Gene(String chr, int start, int end, String name) throws IllegalArgumentException {
        super(chr, start, end, name);
        
        // set default values for annotations
        phenotpyes = new ArrayList<String>();
        strand = ".";
        symbol=".";
    }
    
    /**
     * Constructor for {@link Gene} object.
     * Construct an CNV object by taking all annotations as arguments.
     * 
     * @param chr   Chromosome identifier
     * @param start Start coordinate (0-based)
     * @param end   End coordinate (0-based, half-open)
     * @param name  Name or ID of the CNV
     * @param phenotypes    List of HPO term IDs that are associated with the gene
     */
    public Gene(String chr, int start, int end, String name, List<String> phenotypes){

        // consturct an CVN object using the constructor of the {@link GenomicElement} super calss
        super(chr, start, end, name);
        
        // add annotations
        this.phenotpyes = phenotypes;
        strand = ".";
        symbol=".";
    }
    
}
