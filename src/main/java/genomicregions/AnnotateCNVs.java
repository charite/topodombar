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
     * Annotates CNVs with adjacent regions, defined as the region between the 
     * CNV breakpoint and the end of the underling domain. If a break-point do not
     * lie inside a domain region (but in a domain border or unorganized chromatin),
     * an zero-length region at the breakpoint will be annotated. Note, this function 
     * assumes non-overlapping topological domain regions as input. 
     * 
     * @param cnvs CNVs to annotate
     * @param domains non-overlapping domain regions
     */
    public static void defineAdjacentRegionsByDomains(GenomicSet<CNV> cnvs, GenomicSet<GenomicElement> domains){
        
        // iterate over all CNVs:
        for (CNV cnv : cnvs.values()){
            
            String chr = cnv.getChr();
            int start = cnv.getStart();
            int end = cnv.getEnd();
            
            // get domain regions underling the left (start) and right (end) borders of the CNV
            GenomicSet<GenomicElement> leftDomains = domains.anyOverlap(new GenomicElement(chr, start, start, "cnvStart"));
            GenomicSet<GenomicElement> rightDomains = domains.anyOverlap(new GenomicElement(chr, end, end, "cnvEnd"));
            
            // if left CNV border lies not in a domain region (but in boundary or unorganized chromatin),
            // an zero length region will be defined.
            if (leftDomains.isEmpty()){
                cnv.setLeftAdjacentRegion(new GenomicElement(chr, start, start, "leftAdjacent"));
            }else{
                
                // get start of adjacent region as start of left adjacten domain
                GenomicElement leftDomain = leftDomains.values().iterator().next();
                int leftAdjecentRegionStart = leftDomain.getStart();

                // construct element for the left adjacent regions
                cnv.setLeftAdjacentRegion(new GenomicElement(chr, leftAdjecentRegionStart, cnv.getStart(), "leftAdjacentRegion"));
            }
            
            // same for the right site:
            if(rightDomains.isEmpty()){
                // zero-lenght region
                cnv.setRightAdjacentRegion(new GenomicElement(chr, end, end, "rightAdjacent"));
            }else{
                
                // get end of the query region as end of right adjacten domain
                GenomicElement rightDomain = rightDomains.values().iterator().next();
                int rightAdjecentRegionEnd = rightDomain.getEnd();

                // construct query element for the left adjacent regions
                cnv.setRightAdjacentRegion( new GenomicElement(chr, cnv.getEnd(), rightAdjecentRegionEnd, "rightAdjacentRegion"));
                
            }
        }
    }

    /**
     * Annotates all input CNVs with adjacent genes. Thereby a gene is adjacent,
     * if it is located in the region between the CNV boundary and the end of the
     * underling domain.
     * For each {@link CNV} object the variables {@link CNV.adjacentGenesLeft} 
     * and {@link CNV.adjacentGenesLeft} are filled with a {@link GenomicSet} 
     * of {@link Gene} objects.
     * 
     * @param cnvs
     * @param genes 
     */
//    public static void annotateAdjecentGenes(GenomicSet<CNV> cnvs, GenomicSet<Gene> genes, GenomicSet<GenomicElement> domains){
//        
//        // iterate over all CNVs:
//        for (CNV cnv : cnvs.values()){
//            
//            String chr = cnv.getChr();
//            int start = cnv.getStart();
//            int end = cnv.getEnd();
//            
//            // get domain regions underling the left (start) and right (end) borders of the CNV
//            GenomicSet<GenomicElement> leftDomains = domains.anyOverlap(new GenomicElement(chr, start, start, "cnvStart"));
//            GenomicSet<GenomicElement> rightDomains = domains.anyOverlap(new GenomicElement(chr, end, end, "cnvEnd"));
//            
//            // if left CNV border lies not in a domain region (but in boundary or unorganized chromatin),
//            // no adjacent genes will be assigned.
//            if (leftDomains.isEmpty()){
//                GenomicSet leftAdjacentElements = new GenomicSet();
//            }else{
//                
//                // get start of the query region as start of left adjecten domain
//                GenomicElement leftDomain = leftDomains.values().iterator().next();
//                int leftAdjecentRegionStart = leftDomain.getStart();
//
//                // construct query element for the left adjacent regions
//                GenomicElement leftAdjacentRegion = new GenomicElement(chr, leftAdjecentRegionStart, cnv.getStart(), "leftAdjacentRegion");
//                GenomicSet leftAdjacentElements = genes.completeOverlap(leftAdjacentRegion);
//            }
//            
//            // same for the right site:
//            if(rightDomains.isEmpty()){
//                GenomicSet rightAdjacentElements = new GenomicSet();
//            }else{
//                
//                // get end of the query region as end of right adjecten domain
//                GenomicElement rightDomain = rightDomains.values().iterator().next();
//                int rightAdjecentRegionEnd = rightDomain.getEnd();
//
//                // construct query element for the left adjacent regions
//                GenomicElement rightAdjacentRegion = new GenomicElement(chr, cnv.getEnd(), rightAdjecentRegionEnd, "reightAdjacentRegion");
//                GenomicSet leftAdjacentElements = genes.completeOverlap(rightAdjacentRegion);
//                
//            }
//        }
//    }
//    
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

    public static void phenogramScoreAdjacentGenes(GenomicSet<CNV> cnvs, GenomicSet<Gene> genes, OntologyWrapper ontologyWrapper) {
        for (CNV cnv: cnvs.values()){
            
            // get genes in the left adjacent region
            GenomicSet<Gene> leftAdjacentGenes = genes.completeOverlap(cnv.getLeftAdjacentRegion());

            // calculate phenogram score of these genes
            cnv.setLeftAdjacentPhenogramScore( ontologyWrapper.phenoGramScore( cnv.getPhenotypeTerms(), leftAdjacentGenes) );

            // get genes in the right adjacent region
            GenomicSet<Gene> rightAdjacentGenes = genes.completeOverlap(cnv.getRightAdjacentRegion());

            // calculate phenogram score of these genes
            cnv.setRightAdjacentPhenogramScore( ontologyWrapper.phenoGramScore( cnv.getPhenotypeTerms(), rightAdjacentGenes) );
        
        }
    }
    
    /**
     * This function annotates the CNVs as toplological domain boundary disruption (TDBD).
     * Note, it assumes that the CNVs are are annotated with boundaryOverlap, overlapPhenogramscore, 
     * adjacentPhenogram scores and adjacentRegions.
     * 
     * @param cnvs
     * @param enhancer 
     */
    public static void annotateTDBD(GenomicSet<CNV> cnvs, GenomicSet<GenomicElement> enhancer){
        
        for (CNV cnv: cnvs.values()){
            
            // initialize with False
            cnv.setIsTDBD(false);
            
            // get maximum of left and right phenogram score.
            // TODO: take the left one for enhancer on the right and the other way rund after
            // reproducing the python script results.
            Double maxAdjacentScore = Math.max(cnv.getLeftAdjacentPhenogramScore(), cnv.getRightAdjacentPhenogramScore());

            if ( cnv.isHasBoundaryOverlap() & maxAdjacentScore > 0.0){
                
                // check if enhancers are available in left and/or right adjacent regions
                boolean hasLeftEnhancer = ! enhancer.completeOverlap(cnv.getLeftAdjacentRegion()).isEmpty();
                boolean hasRightEnhancer = ! enhancer.completeOverlap(cnv.getRightAdjacentRegion()).isEmpty();
            
                // enhancer in left adjacent region and right phenoGram > 0
                if(hasLeftEnhancer & cnv.getRightAdjacentPhenogramScore() > 0){
                    cnv.setIsTDBD(true);
                }

                // enhancer in right adjacent region and left phenoScore > 0
                if(hasRightEnhancer & cnv.getLeftAdjacentPhenogramScore() > 0){
                    cnv.setIsTDBD(true);
                }
                
            }
        }
        
    }
}
