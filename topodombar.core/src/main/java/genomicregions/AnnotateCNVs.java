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

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import ontologizer.go.Term;
import phenotypeontology.PhenotypeData;

/**
 * This class implements functionality to annotate CNVs and predict their most
 * likely effect mechanisms.
 *
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class AnnotateCNVs {
    
    
    /**
     * This function annotates the input CNVs with overlapped topological domain
     * boundaries and genes as well as genes and enhancers in the adjacent regions
     * of the CNVs and calculates the phenogram score of genes within and adjacent of the
     * CNV. Note, this function assumes that the input CNVs have already defined
     * adjacent regions.
     * 
     * @param cnvs
     * @param boundaries
     * @param genes
     * @param enhancers
     * @param phenotypeData 
     */
    public static void annoateCNVsForTDBD(GenomicSet<CNV> cnvs, GenomicSet<GenomicElement> domains, 
            GenomicSet<GenomicElement> boundaries, GenomicSet<Gene> genes, 
            GenomicSet<GenomicElement> enhancers, PhenotypeData phenotypeData){
        
        // define adjacent regions as the interval form CNV breakpoint to the end
        // of the underlying toplological domain.
        defineAdjacentRegionsByDomains(cnvs, domains);
        
        // annotate boundaries:
        boundaryOverlap(cnvs, boundaries);
        
        // annotate CNVs with genes that are completely overlapped by the CNV
        AnnotateCNVs.annotateOverlappedGenes(cnvs, genes);

        // annotate CNVs with genes that lie in the adjacent regions
        AnnotateCNVs.annotateAdjacentGenes(cnvs, genes);
        
        // compute phenogram score for overlapped and adjacent genes:
        AnnotateCNVs.phenogramScore(cnvs, phenotypeData);
        
        // annotate CNVs with adjacent enhancers
        AnnotateCNVs.annotateAdjacentEnhancers(cnvs, enhancers);
    }
    
    /**
     * This function annotates the input CNVs with overlapped topological domain
     * boundaries and genes together with the corresponding phenogram score of 
     * these genes overlapped by the CNV.
     * 
     * @param cnvs
     * @param boundaries
     * @param genes
     * @param enhancers
     * @param phenotypeData 
     */
    public static void annoateOverlap(GenomicSet<CNV> cnvs, 
            GenomicSet<GenomicElement> boundaries, GenomicSet<Gene> genes, 
            GenomicSet<GenomicElement> enhancers, PhenotypeData phenotypeData){
        
        // annotate boundaries:
        boundaryOverlap(cnvs, boundaries);
        
        // annotate CNVs with genes that are completely overlapped by the CNV
        AnnotateCNVs.annotateOverlappedGenes(cnvs, genes);

        // compute phenogram score for overlapped genes:
        AnnotateCNVs.overlapPhenogramScore(cnvs, phenotypeData);
    }

    /**
     * This function annotates the adjacent regions of CNVs with genes, the corresponding
     * phenogram score and enhancer elements.
     * Note, that this function assumes that the adjacent regions are already
     * defined for all input CNVs.
     * @see defineAdjacentRegionsByDomains
     * @see defineAdjacentRegionsByDistance 
     * 
     * @param cnvs
     * @param genes
     * @param enhancers
     * @param phenotypeData 
     */
    public static void annoateAdjacentRegions(GenomicSet<CNV> cnvs, GenomicSet<Gene> genes, 
            GenomicSet<GenomicElement> enhancers, PhenotypeData phenotypeData){

        // annotate CNVs with genes that lie in the adjacent regions
        AnnotateCNVs.annotateAdjacentGenes(cnvs, genes);
        
        // compute phenogram score for overlapped and adjacent genes:
        AnnotateCNVs.adjacentPhenogramScore(cnvs, phenotypeData);
        
        // annotate CNVs with adjacent enhancers
        AnnotateCNVs.annotateAdjacentEnhancers(cnvs, enhancers);
        
    }
    
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
        }
    }
    
    /**
     * Annotates all input CNVs with all genes that have any overlap with the CNV.
     * For each {@link CNV} object the variable {@link CNV.geneOverlap} is filled 
     * with a {@link GenomicSet} of {@link Gene} objects.
     * 
     * @param cnvs CNVs that should be annotated
     * @param genes Set of genes
     */
    public static void annotateOverlappedGenes(GenomicSet<CNV> cnvs, GenomicSet<Gene> genes){
        
        // iterate over all CNVs:
        for (CNV cnv : cnvs.values()){
            
            GenomicSet<Gene> overlap = genes.anyOverlap(cnv);
            cnv.setGenesInOverlap( overlap );
            
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
            GenomicSet<GenomicElement> leftDomains = domains.anyOverlap(new GenomicElement(chr, start-1, start, "cnvStart"));
            GenomicSet<GenomicElement> rightDomains = domains.anyOverlap(new GenomicElement(chr, end, end+1, "cnvEnd"));
                        
            // if left CNV border lies not in a domain region (but in boundary or unorganized chromatin),
            // an zero length region will be defined.
            if (leftDomains.isEmpty()){
                cnv.setLeftAdjacentRegion(new GenomicElement(chr, start, start, "leftAdjacentRegion"));
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
                cnv.setRightAdjacentRegion(new GenomicElement(chr, end, end, "rightAdjacentRegion"));
            }else{
                
                // get end of the query region as end of right adjacten domain
                GenomicElement rightDomain = rightDomains.values().iterator().next();
                int rightAdjecentRegionEnd = rightDomain.getEnd();

                // construct query element for the right adjacent regions
                cnv.setRightAdjacentRegion( new GenomicElement(chr, cnv.getEnd(), rightAdjecentRegionEnd, "rightAdjacentRegion"));
                
            }
        }
    }

    /**
     * Annotates CNVs with adjacent regions, defined as windows with fixed size 
     * on each side of the CNV.
     * 
     * @param cnvs the CNVs to be annotated
     * @param regionSize size of the adjacent regoin in base pairs (bp)
     */
    public static void defineAdjacentRegionsByDistance(GenomicSet<CNV> cnvs, int regionSize){
        
        for (CNV cnv : cnvs.values()){

            int start = cnv.getStart();
            int end = cnv.getEnd();
            String chr = cnv.getChr();
            
            // catch negative coordinate start values if CNV lies to near to the chromosom start.
            int leftRegStart = Math.max(0, start-regionSize);
            
            // construct query element for the left adjacent regions
            cnv.setLeftAdjacentRegion(new GenomicElement(chr, leftRegStart, start, "leftAdjacent"));
            // construct query element for the right adjacent regions
            cnv.setRightAdjacentRegion(new GenomicElement(chr, end, end+regionSize, "rightAdjacent"));
        }
    }
    
    /**
     * Annotates all input CNVs with genes in the left and right adjacent regions of the CNV.
     * For each {@link CNV} object the variable {@link CNV.genesInLeftRegion} 
     * and {@link CNV.genesInRightRegion} is filled with a {@link GenomicSet} of {@link Gene} objects.
     * This function assumes pre-defined the the adjacent regions 
     * (e.g by the function {@link defineAdjacentRegionsByDomains}).  
     * Note, this function requires that a gene lies completely in the adjacent regions,
     * that is the adjacent regions overlaps completely the entire gene.
     * 
     * @param cnvs CNVs that should be annotated
     * @param genes Set of genes
     */    
    public static void annotateAdjacentGenes(GenomicSet<CNV> cnvs, GenomicSet<Gene> genes){

        // iterate over all input CNVs:
        for (CNV cnv : cnvs.values()){

            // get all genes in the left adjacent regions
            GenomicSet<Gene> leftGenes = genes.completeOverlap(cnv.getLeftAdjacentRegion());
            cnv.setGenesInLeftRegion( leftGenes );

            // get all genes in the right adjacent regions
            GenomicSet<Gene> rightGenes = genes.completeOverlap(cnv.getRightAdjacentRegion());
            cnv.setGenesInRightRegion( rightGenes );
            
        }
        
    }

    /**
     * Compute the phenogram score for genes overlapped by the input {@link CNV}s and
     * for genes in the left and right adjacent regions.
     * It writes the memeber variables {@link CNV.overlapPhenogramScore}, 
     * {@link CNVleftAdjacentPhenogramScore} and 
     * {@link CNVrightAdjacentPhenogramScore} in each {@link CNV} object
     * Note, this function assumes that the CNVs are annotated with overlped and
     * adjacent genes by the functions {@link annotateOverlappedGenes} and
     * {@link annotateAdjacentGenes}.
     * 
     * @param cnvs CNVs for which the phenogram score should be calculated.
     * @param phenotypeData the phenotype ontology
     */
    public static void phenogramScore(GenomicSet<CNV> cnvs, PhenotypeData phenotypeData){
        
        for (CNV cnv: cnvs.values()){
            // overlap PhenogramScore
            cnv.setOverlapPhenogramScore( 
                phenotypeData.phenoGramScore( cnv.getPhenotypes(), cnv.getGenesInOverlap()) 
            );
        
            // left adjacent PhenogramScore:
            cnv.setLeftAdjacentPhenogramScore( 
                phenotypeData.phenoGramScore( cnv.getPhenotypes(), cnv.getGenesInLeftRegion()) 
            );

            // right adjacent PhenogramScore:
            cnv.setRightAdjacentPhenogramScore(
               phenotypeData.phenoGramScore( cnv.getPhenotypes(), cnv.getGenesInRightRegion())
            );
        }
    }
    
    /**
     * Compute the phenogram score for genes overlapped by the input {@link CNV}s.
     * It writes the memeber variables {@link CNV.overlapPhenogramScore} in each {@link CNV} object
     * Note, this function assumes that the CNVs are annotated with overlped and
     * adjacent genes by the functions {@link annotateOverlappedGenes}.
     * 
     * @param cnvs CNVs for which the phenogram score should be calculated.
     * @param phenotypeData the phenotype ontology
     */
    public static void overlapPhenogramScore(GenomicSet<CNV> cnvs, PhenotypeData phenotypeData){
        
        for (CNV cnv: cnvs.values()){
            // overlap PhenogramScore
            cnv.setOverlapPhenogramScore( 
                phenotypeData.phenoGramScore( cnv.getPhenotypes(), cnv.getGenesInOverlap()) 
            );
        }
    }
    
    /**
     * Compute the phenogram score for genes in the left and right adjacent regions.
     * It writes the memeber variables {@link CNVleftAdjacentPhenogramScore} and 
     * {@link CNVrightAdjacentPhenogramScore} in each {@link CNV} object
     * Note, this function assumes that the CNVs are annotated adjacent genes by 
     * the function {@link annotateAdjacentGenes}.
     * 
     * @param cnvs CNVs for which the phenogram score should be calculated.
     * @param phenotypeData the phenotype ontology
     */
    public static void adjacentPhenogramScore(GenomicSet<CNV> cnvs, PhenotypeData phenotypeData){
        
        for (CNV cnv: cnvs.values()){
        
            // left adjacent PhenogramScore:
            cnv.setLeftAdjacentPhenogramScore( 
                phenotypeData.phenoGramScore( cnv.getPhenotypes(), cnv.getGenesInLeftRegion()) 
            );

            // right adjacent PhenogramScore:
            cnv.setRightAdjacentPhenogramScore(
               phenotypeData.phenoGramScore( cnv.getPhenotypes(), cnv.getGenesInRightRegion())
            );
        }
    }
    
    /**
     * Annotates all input CNVs with enhancers in the left and right adjacent regions of the CNV.
     * For each {@link CNV} object the variable {@link CNV.enhancersInLeftRegion} 
     * and {@link CNV.enhancersInRightRegion} is filled with a 
     * {@link GenomicSet} of {@link GenomicElemnt} objects representing the enhancers.
     * Note, his function assumes pre-defined the the adjacent regions 
     * (e.g by the function {@link defineAdjacentRegionsByDomains}).  
     * Note, this function requires that an enhancer lies completely in the adjacent regions,
     * that is the adjacent regions overlaps completely the entire enhancer.
     * 
     * @param cnvs CNVs that should be annotated
     * @param enhancers Set of enhancers
     */    
    public static void annotateAdjacentEnhancers(GenomicSet<CNV> cnvs, GenomicSet<GenomicElement> enhancers){

        // iterate over all input CNVs:
        for (CNV cnv : cnvs.values()){

            // get all enhancers in the left adjacent regions
            GenomicSet<GenomicElement> leftEnhancers = enhancers.completeOverlap(cnv.getLeftAdjacentRegion());
            cnv.setEnhancersInLeftRegion( leftEnhancers );

            // get all enhancers in the right adjacent regions
            GenomicSet<GenomicElement> rightEnhancers = enhancers.completeOverlap(cnv.getRightAdjacentRegion());
            cnv.setEnhancersInRightRegion( rightEnhancers );
            
        }
        
    }
    
    /**
     * This function annotates the CNVs as topological domain boundary 
     * disruption (TDBD) by the new definition, that does not require the 
     * concept of target terms.
     * Note, it assumes that the CNVs are are annotated with boundaryOverlap, overlapPhenogramscore, 
     * adjacentPhenogram scores and adjacentRegions. It defines target genes just by
     * any semantic similarity greater zero of all patient terms to the terms associated with genes.
     * 
     * @param cnvs
     */
    public static void annotateTDBDjustByScore(GenomicSet<CNV> cnvs){
        
        for (CNV cnv: cnvs.values()){
            
            // initialize with False
            boolean isTDBD = false;
            
            if ( cnv.hasBoundaryOverlap()){
                
                // check if enhancers are available in left and/or right adjacent regions
                boolean hasLeftEnhancer = ! cnv.getEnhancersInLeftRegion().isEmpty();
                boolean hasRightEnhancer = ! cnv.getEnhancersInRightRegion().isEmpty();
            
                // enhancer in left adjacent region and right phenoGram > 0
                if(hasLeftEnhancer & cnv.getRightAdjacentPhenogramScore() > 0){
                    isTDBD = true;
                }

                // enhancer in right adjacent region and left phenoScore > 0
                if(hasRightEnhancer & cnv.getLeftAdjacentPhenogramScore() > 0){
                    isTDBD = true;
                }
                
            }
            
            // calc maximum of adjacent scores
            Double maxAdjacentScore = Math.max(
                    cnv.getLeftAdjacentPhenogramScore(), 
                    cnv.getRightAdjacentPhenogramScore()
                );
            
            // set newTDBD effect mechanism 
            if (isTDBD){
                // if ther is evidance for TDBD decide if, mixed or TDBD
                if(maxAdjacentScore > cnv.getOverlapPhenogramScore()){
                    cnv.setEffectMechanism("newTDBD", "TDBD");
                }else{
                    cnv.setEffectMechanism("newTDBD", "Mixed");
                }                
            }else{
                // if no evidance for TDBD decide if GDE or no data
                if(cnv.getOverlapPhenogramScore() > 0){
                    cnv.setEffectMechanism("newTDBD", "GDE");
                }else{
                    cnv.setEffectMechanism("newTDBD", "NoData");
                }
            }
            
        }
        
    }
    
    /**
     * This function annotates the CNVs as toplological domain boundary disruption (TDBD).
     * Note, it assumes that the CNVs are are annotated with boundaryOverlap, overlapPhenogramscore, 
     * adjacentPhenogram scores and adjacentRegions.
     * Here, target Genes are defined as genes that are associated to the (more general) target Term of 
     * the patient.
     * 
     * @param cnvs
     * @param term2genes 
     */
    public static void annotateTDBD(GenomicSet<CNV> cnvs,  HashMap<Term,HashSet<String>> term2genes){
        
        for (CNV cnv: cnvs.values()){


            // define target gnes as gens associated with the target term
            HashSet<String> allTargetGenes = term2genes.get(cnv.getTargetTerm());
            
            // is TDBD if cnv overlaps a boundary element and has matching gene
            // in one adjacent region and enhancer in the other adjacent region.
            boolean isTDBD = cnv.hasBoundaryOverlap() && matchingEnhancerAndGene(cnv, allTargetGenes);
                
            
            // Decide to which evvect mechanism class the CNV corresponds to
            if (isTDBD){

                // take the maximum of left and right adjacent phenogram score
                Double maxAdjacentScore = Math.max(cnv.getLeftAdjacentPhenogramScore(), cnv.getRightAdjacentPhenogramScore());

                // if the score in the adjacent regions are higher compaired to the overlapped regions
                if  (maxAdjacentScore > cnv.getOverlapPhenogramScore()){
                    // TDBD only categorie
                    cnv.setEffectMechanism("TDBD", "TDBD");
                }else{

                    // TDBD and GDE evidence
                    cnv.setEffectMechanism("TDBD", "Mixed");
                }

            }else{
                // if no TDBD or TDBD_only can be assigned test, if the CNV can
                // be explained by GDE:
                if (cnv.getOverlapPhenogramScore() > 0){
                    cnv.setEffectMechanism("TDBD", "GDE");
                // ther is no efidence for TDBD or GDE, assign CNV to "No_Data" categorie.
                }else{
                    cnv.setEffectMechanism("TDBD", "NoData");
                }
            }            
        }    
    }
    
    /**
     * This function annotates the CNVs as enhancer adoption mechanism (EA) without
     * the requirement to overlap a topological domain boundary.
     * Note, it assumes that the CNVs are are annotated with overlapPhenogramscore, 
     * adjacentPhenogram, genes in overlap region and adjacent region and
     * enhancers in adjacent regions.
     * Here, target Genes are defined as genes that are associated to the (more general) target Term of 
     * the patient.
     * 
     * @param cnvs
     * @param enhancer 
     * @param genes 
     * @param term2genes 
     */
    public static void annotateEA(GenomicSet<CNV> cnvs,  HashMap<Term,HashSet<String>> term2genes){
        
        for (CNV cnv: cnvs.values()){
            
            // define target gnes as gens associated with the target term
            HashSet<String> allTargetGenes = term2genes.get(cnv.getTargetTerm());
            
            // Decide to which evvect mechanism class the CNV corresponds to
            if (matchingEnhancerAndGene(cnv, allTargetGenes)){

                // take the maximum of left and right adjacent phenogram score
                Double maxAdjacentScore = Math.max(cnv.getLeftAdjacentPhenogramScore(), cnv.getRightAdjacentPhenogramScore());

                // if the score in the adjacent regions are higher compaired to the overlapped regions
                if  (maxAdjacentScore > cnv.getOverlapPhenogramScore()){
                    // TDBD only categorie
                    cnv.setEffectMechanism("EA", "EA");
                }else{

                    // TDBD and GDE evidence
                    cnv.setEffectMechanism("EA", "Mixed");
                }

            }else{
                // if no TDBD or TDBD_only can be assigned test, if the CNV can
                // be explained by GDE:
                if (cnv.getOverlapPhenogramScore() > 0){
                    cnv.setEffectMechanism("EA", "GDE");
                // ther is no efidence for TDBD or GDE, assign CNV to "No_Data" categorie.
                }else{
                    cnv.setEffectMechanism("EA", "NoData");
                }
            }            
        }    
    }
    
    /**
     * This function annotates the CNVs as enhancer adoption mechanism (EA) with 
     * low evidence for GDE.
     * In contrast to EA, here a phenogram score > 0 for the overlap region is 
     * required too.
     * Note, it assumes that the CNVs are are annotated with overlapPhenogramscore, 
     * adjacentPhenogram, genes in overlap region and adjacent region and
     * enhancers in adjacent regions.
     * Here, target Genes are defined as genes that are associated to the (more general) target Term of 
     * the patient.
     * 
     * @param cnvs
     * @param enhancer 
     * @param genes 
     * @param term2genes 
     */
    public static void annotateEAlowG(GenomicSet<CNV> cnvs,  HashMap<Term,HashSet<String>> term2genes){
        
        for (CNV cnv: cnvs.values()){
            
            // define target gnes as gens associated with the target term
            HashSet<String> allTargetGenes = term2genes.get(cnv.getTargetTerm());
            
            // Decide to which evvect mechanism class the CNV corresponds to
            if (matchingEnhancerAndGene(cnv, allTargetGenes)){

                // take the maximum of left and right adjacent phenogram score
                Double maxAdjacentScore = Math.max(cnv.getLeftAdjacentPhenogramScore(), cnv.getRightAdjacentPhenogramScore());

                // if the score in the adjacent regions are higher compaired to the overlapped regions
                // but ther is a score > 0 in the overlap region than it is EAlowG
                if  (maxAdjacentScore > cnv.getOverlapPhenogramScore() && cnv.getOverlapPhenogramScore() > 0){
                    // TDBD only categorie
                    cnv.setEffectMechanism("EAlowG", "EAlowG");
                }else{
                    
                    // TDBD and GDE evidence
                    cnv.setEffectMechanism("EAlowG", "Mixed");
                }

            }else{
                // if no TDBD or TDBD_only can be assigned test, if the CNV can
                // be explained by GDE:
                if (cnv.getOverlapPhenogramScore() > 0){
                    cnv.setEffectMechanism("EAlowG", "GDE");
                // ther is no efidence for TDBD or GDE, assign CNV to "No_Data" categorie.
                }else{
                    cnv.setEffectMechanism("EAlowG", "NoData");
                }
            }            
        }    
    }

    /**
     * This function defines the overlapped inner-domain regions by the CNV.
     * There is a region for each CNV breakpoint that is defined from the break-point
     * to the next boundary that is overlapped by the cnv.
     * <pre>
Domains:               /'''''''''\   /'''\  /''''''''''''\
Boundary:                         ---     --       
CNV:                         ======================
overlappedDomainRegions:     *****          *******
     </pre>
     * These regions are needed to interpret tandem duplications and inversions for 
     * enhancer adoption mechanisms (e.g relevant enhancer in left overlap region
     * might come close and interact with a duplicated gene in the right overlapped region).
     * 
     * @param cnvs
     * @param domains 
     */
    public static void defineOverlappedDomainRegions(GenomicSet<CNV> cnvs, GenomicSet<GenomicElement> domains){
        
        for (CNV cnv: cnvs.values()){
            
            int start = cnv.getStart();
            int end = cnv.getEnd();
            String chr = cnv.getChr();

            // check if CNV ovelaps at leaste one boundary element
            if (cnv.hasBoundaryOverlap()){
                
                // get domain regions underling the left (start) and right (end) borders of the CNV
                GenomicSet<GenomicElement> leftDomains = domains.anyOverlap(new GenomicElement(chr, start-1, start, "cnvStart"));
                GenomicSet<GenomicElement> rightDomains = domains.anyOverlap(new GenomicElement(chr, end, end+1, "cnvEnd"));

                // if left CNV border lies not in a domain region (but in boundary or unorganized chromatin),
                // an zero length region will be defined.
                if (leftDomains.isEmpty()){
                    cnv.setLeftOverlappedDomainRegion(new GenomicElement(chr, start, start, "leftOverlapped"));
                }else{

                    // get end of overlapped region as end of the domain overlapped by the CNV start
                    GenomicElement leftDomain = leftDomains.values().iterator().next();
                    int leftOverlapRegionEnd = leftDomain.getEnd();

                    // construct element for the left adjacent regions
                    cnv.setLeftOverlappedDomainRegion(new GenomicElement(chr, cnv.getStart(), leftOverlapRegionEnd, "leftOverlapped"));
                }

                // same for the right site:
                if(rightDomains.isEmpty()){
                    // zero-lenght region
                    cnv.setRightOverlappedDomainRegion(new GenomicElement(chr, end, end, "rightOverlapped"));
                }else{

                    // get start of the query region as start of right adjacten domain
                    GenomicElement rightDomain = rightDomains.values().iterator().next();
                    int rightOverlapRegionStart = rightDomain.getStart();

                    // construct query element for the right adjacent regions
                    cnv.setRightOverlappedDomainRegion( new GenomicElement(chr, rightOverlapRegionStart, cnv.getEnd(), "rightOverlapped"));

                }
            
            // in case of no boundary overlap set default regions of zero size
            }else{
                cnv.setLeftOverlappedDomainRegion(new GenomicElement(chr, start, start, "leftOverlapped"));
                cnv.setRightOverlappedDomainRegion(new GenomicElement(chr, end, end, "rightOverlapped"));
            }
        }
    }
    
    /**
     * Explains the input tandem duplications by enhancer adoption (TanDupEA)
     * mechanisms. Therefore a domain boundary has to be overlapoed and furthermore a 
     * phenotypic relevant gene has to lie in the left overlapped domain region
     * and an enhancer in the right overlapped domain region (or the other way around).
     * For the definition of the overlapped regions see {@link defineOverlappedDomainRegions}.
     * By assuming the duplication in tandem, the enhancer 'moves' into the vicinity 
     * of the gene and can interact with it, since no boundary element separate them.
     * 
     * @param cnvs CNVs to be annotated, duplications (type = "gain") are assumed to be tandem duplications
     * @param genes genes
     * @param enhancers enhancer elements or regulatory regions of interest
     * @param phenotypeData phenotype data
     */
    public static void tandemDuplicationEnhancerAdoption(GenomicSet<CNV> cnvs,
            GenomicSet<Gene> genes, GenomicSet<GenomicElement> enhancers, PhenotypeData phenotypeData){
        
        for (CNV cnv : cnvs.values()){
            
            // check if type of the CNV is a duplication
            if ("gain".equals(cnv.getType())){
                
                /**
                 * Note, there is no need to test for boundary overlap, since in 
                 * case the CNV does not overlap a boundary the overlapedDomainRegions
                 * are of zero size and can therefore not contain enhancers or genes.
                 * {@see defineOverlappedDomainRegions}.
                */
                
                // get genes in left and right overlapped domain regions
                GenomicSet<Gene> leftGenes = genes.completeOverlap(cnv.getLeftOverlappedDomainRegion());
                GenomicSet<Gene> rightGenes = genes.completeOverlap(cnv.getRightOverlappedDomainRegion());
                
                // calculate phenogram scores separately for genes in left and right regions
                Double leftPhenogramScore = phenotypeData.phenoGramScore(cnv.getPhenotypes(), leftGenes);
                Double rightPhenogramScore = phenotypeData.phenoGramScore(cnv.getPhenotypes(), rightGenes);
                
                // get enhancers in left and right overlapped domain regions:
                GenomicSet<GenomicElement> leftEnhancers = enhancers.completeOverlap(cnv.getLeftOverlappedDomainRegion());
                GenomicSet<GenomicElement> rightEnhancers = enhancers.completeOverlap(cnv.getRightOverlappedDomainRegion());
                
                
                // test for enhancer adoption mechanism
                boolean TanDupEA = 
                        // enhancer left and relevant gene right
                        (! leftEnhancers.isEmpty() && rightPhenogramScore > 0)

                        || // or gene left and enhancer right                        
                        (leftPhenogramScore > 0 && ! rightEnhancers.isEmpty());
                    
                if (TanDupEA){
                    // set effectmechansim to "TanDupEA"
                    cnv.setEffectMechanism("TanDupEA", "TanDupEA");
                    
                // if no evidance for TanDupEA is found   
                }else{
                    
                    // check if gene dosag effect can be assigned 
                    if (cnv.getOverlapPhenogramScore() > 0){
                        
                        cnv.setEffectMechanism("TanDupEA", "onlyGDE");                        
                    
                    }else{
                    
                        // set the mechanism to "NoData" 
                        cnv.setEffectMechanism("TanDupEA", "NoData");
                    }
                }
                
           
            // if CNV type does not fit, set annotation to not available (NA).
            }else{
                cnv.setEffectMechanism("TanDupEA", "NA");
            }
        }
    }
    
    public static void inversionEnhancerAdoption(GenomicSet<CNV> cnvs,
            GenomicSet<Gene> genes, GenomicSet<GenomicElement> enhancers, PhenotypeData phenotypeData){

        for (CNV cnv : cnvs.values()){

            // check if type of the CNV is an inversion
            if ("inversion".equals(cnv.getType())){
                
                // check if CNV overlaps a domain boundary:
                if (cnv.hasBoundaryOverlap()){
                    
                    /*
                        first, consider an inverted enhancer that comes close
                        to an adjacent gene on the other side.
                    */
                    
                    // get enhancers in left and right overlapped domain regions:
                    GenomicSet<GenomicElement> leftOverlappedEnhancers = enhancers.completeOverlap(cnv.getLeftOverlappedDomainRegion());
                    GenomicSet<GenomicElement> rightOverlappedEnhancers = enhancers.completeOverlap(cnv.getRightOverlappedDomainRegion());
                    
                    boolean invEnhancer =  (
                            
                        // enhancer in left overlapped region and right adjacent gene is relevant
                        (!leftOverlappedEnhancers.isEmpty() && cnv.getRightAdjacentPhenogramScore() > 0)
                        
                        // OR enhancer in right overlapped region and left adjacent gene is relevant
                        || (!rightOverlappedEnhancers.isEmpty() && cnv.getLeftAdjacentPhenogramScore() > 0)
                            );

                    
                    /*
                        secondly, consider the case that a relevant gene is inverted
                        to the other side of the boundary (or boundaries).
                    */
                    
                    // check if enhancers are available in left and/or right adjacent regions
                    boolean hasLeftEnhancer = ! cnv.getEnhancersInLeftRegion().isEmpty();
                    boolean hasRightEnhancer = ! cnv.getEnhancersInRightRegion().isEmpty();
                    
                    // get genes in left and right overlapped domain regions
                    GenomicSet<Gene> leftOverlappedGenes = genes.completeOverlap(cnv.getLeftOverlappedDomainRegion());
                    GenomicSet<Gene> rightOverlappedGenes = genes.completeOverlap(cnv.getRightOverlappedDomainRegion());

                    // calculate phenogram scores separately for genes in left and right regions
                    Double leftOverlppedPhenogramScore = phenotypeData.phenoGramScore(cnv.getPhenotypes(), leftOverlappedGenes);
                    Double rightOverlappedPhenogramScore = phenotypeData.phenoGramScore(cnv.getPhenotypes(), rightOverlappedGenes);
                    
                    boolean invGene = ( 
                        
                        // enhancer in left adjacent region and right overlaped gene is relevant
                        (hasLeftEnhancer && rightOverlappedPhenogramScore > 0)
                        
                        // OR enhancer in right adjacent region and left ovleraped gene is relevant
                        || (hasRightEnhancer && leftOverlppedPhenogramScore > 0)
                            );
//                    System.out.println("DEBUG: name : " + cnv.getName());
//                    System.out.println("DEBUG: left E right G:" + hasLeftEnhancer + rightOverlappedPhenogramScore);
//                    System.out.println("DEBUG: right E left G:" + hasRightEnhancer + leftOverlppedPhenogramScore);
//                    System.out.println("DEBUG: invGene : " + invGene);
//                    
                    // set the proper annotation for all combinations:
                    if (!invEnhancer && !invGene){
                        cnv.setEffectMechanism("InvEA", "noInvEA");             
                    }

                    if (!invEnhancer && invGene){
                        cnv.setEffectMechanism("InvEA", "GeneInvEA");             
                    }
                    
                    if (invEnhancer && !invGene){
                        cnv.setEffectMechanism("InvEA", "EnhancerInvEA");   
                    }
                    
                    if (invEnhancer && invGene){
                        cnv.setEffectMechanism("InvEA", "EandGInvEA");   
                    }
                
                // in case of no boundary overlap: 
                }else{
                    cnv.setEffectMechanism("InvEA", "noTDB");                
                }
                
            // if CNV type does not fit, set annotation to not available (NA).
            }else{
                cnv.setEffectMechanism("InvEA", "NA");
            }
        }
        
    }

    
    /**
     * Checks for a single {@link CNV} if there is a signature of matching gene(s)
     * in one adjacent region and enhancer(s) in the other adjacent regions.
     * A gene matches, if it is contained in the input set of target genes.
     * Note this function assumes, that the CNV is annotated with adjacent genes
     * and adjacent enhancers. 
     * 
     * @param cnv the {@link CNV} object to be tested
     * @param targetGenes a set of Entrez Gene IDs that are associated with the patients phenotype
     * @return ture if a matching signature is found
     */
    private static boolean matchingEnhancerAndGene(CNV cnv, HashSet<String> targetGenes){

        // initialize with False
        boolean hasSignature = false;

        // get target genes in the adjacent regions.
        Set<String> leftTargetGenes = cnv.getGenesInLeftRegion().keySet();
        leftTargetGenes.retainAll(targetGenes);

        Set<String> rightTargetGenes = cnv.getGenesInRightRegion().keySet();
        rightTargetGenes.retainAll(targetGenes);

        // check if enhancers are available in left and/or right adjacent regions
        boolean hasLeftEnhancer = ! cnv.getEnhancersInLeftRegion().isEmpty();
        boolean hasRightEnhancer = ! cnv.getEnhancersInRightRegion().isEmpty();

        // enhancer in left adjacent region and target gene in the right adjacent region
        if(hasLeftEnhancer & ! rightTargetGenes.isEmpty()){
            hasSignature = true;
        }

        // enhancer in right adjacent region and target gene in left adjacent regions
        if(hasRightEnhancer & !leftTargetGenes.isEmpty()){
             hasSignature = true;
        }
        
        return hasSignature;
    }
 
    
    
}
