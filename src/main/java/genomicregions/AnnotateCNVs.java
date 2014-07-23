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

import phenotypeontology.OntologyWrapper;

/**
 * This class implements functionality to annotate CNVs.
 *
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class AnnotateCNVs {
    
    /**
     * Annotates CNVs with overlapping boundaries.
     * 
     * @param cnvs  copy number variations (CVNs)
     * @param boundaries Topological domain boundaries
     */
    public static void boundaryOverlap(GenomicSet<CNV> cnvs, GenomicSet<GenomicElement> boundaries){
        
        // iterate over all CNVs:
        for (CNV cnv : cnvs.values()){
            
            GenomicSet<GenomicElement> overlap = boundaries.completeOverlap(cnv);
            cnv.setBoundaryOverlap( overlap );
            cnv.setHasBoundaryOverlap( cnv.getBoundaryOverlap().isEmpty() );            
        }
    }
    
    /**
     * Annotates all input CNVs with all genes that have any overlap with the CNV.
     * For each {@link CNV} object the variable {@link CNV.geneOverlap} is filled 
     * with a {@link GenomicSet} of {@link Gene} objects.
     * 
     * @param cnvs
     * @param genes 
     */
    public static void geneOverlap(GenomicSet<CNV> cnvs, GenomicSet<Gene> genes){
        
        // iterate over all CNVs:
        for (CNV cnv : cnvs.values()){
            
            GenomicSet<Gene> overlap = genes.anyOverlap(cnv);
            cnv.setGeneOverlap( overlap );
            
        }
    }
    
    /**
     * Compute the phenogram score for genes overlapped by the input {@link CNV}s.
     * It writes the memeber variables {@link CNV.overlapPhenogramScore} in each {@link CNV} object
     * 
     * @param cnvs CNVs for which the phenogram score should be calculated.
     * @param ontolgyWrapper the phenotype ontology
     */
    public static void phenogramScore(GenomicSet<CNV> cnvs, OntologyWrapper ontolgyWrapper){
        
        for (CNV cnv: cnvs.values()){
            cnv.setOverlapPhenogramScore( ontolgyWrapper.phenoGramScore( cnv.getPhenotypeTerms(), cnv.getGeneOverlap()) );
        }
    }
}
