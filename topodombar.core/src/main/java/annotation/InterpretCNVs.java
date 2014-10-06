/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package annotation;

import static annotation.AnnotateCNVs.defineOverlappedDomainRegions;
import genomicregions.CNV;
import genomicregions.Gene;
import genomicregions.GenomicElement;
import genomicregions.GenomicSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import ontologizer.go.Term;
import phenotypeontology.PhenotypeData;

/**
 *
 * @author jonas
 */
public class InterpretCNVs {

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
    private static boolean matchingEnhancerAndGene(CNV cnv, HashSet<String> targetGenes) {
        // initialize with False
        boolean hasSignature = false;
        Set<String> leftTargetGenes = cnv.getGenesInLeftRegion().keySet();
        leftTargetGenes.retainAll(targetGenes);
        Set<String> rightTargetGenes = cnv.getGenesInRightRegion().keySet();
        rightTargetGenes.retainAll(targetGenes);
        boolean hasLeftEnhancer = !cnv.getEnhancersInLeftRegion().isEmpty();
        boolean hasRightEnhancer = !cnv.getEnhancersInRightRegion().isEmpty();
        if (hasLeftEnhancer & !rightTargetGenes.isEmpty()) {
            hasSignature = true;
        }
        if (hasRightEnhancer & !leftTargetGenes.isEmpty()) {
            hasSignature = true;
        }
        return hasSignature;
    }

    public static void inversionEnhancerAdoption(GenomicSet<CNV> cnvs, GenomicSet<Gene> genes, GenomicSet<GenomicElement> enhancers, PhenotypeData phenotypeData) {
        for (CNV cnv : cnvs.values()) {
            if ("inversion".equals(cnv.getType())) {
                if (cnv.hasBoundaryOverlap()) {
                    GenomicSet<GenomicElement> leftOverlappedEnhancers = enhancers.completeOverlap(cnv.getLeftOverlappedDomainRegion());
                    GenomicSet<GenomicElement> rightOverlappedEnhancers = enhancers.completeOverlap(cnv.getRightOverlappedDomainRegion());
                    boolean invEnhancer = (!leftOverlappedEnhancers.isEmpty() && cnv.getRightAdjacentPhenogramScore() > 0) || (!rightOverlappedEnhancers.isEmpty() && cnv.getLeftAdjacentPhenogramScore() > 0);
                    boolean hasLeftEnhancer = !cnv.getEnhancersInLeftRegion().isEmpty();
                    boolean hasRightEnhancer = !cnv.getEnhancersInRightRegion().isEmpty();
                    GenomicSet<Gene> leftOverlappedGenes = genes.completeOverlap(cnv.getLeftOverlappedDomainRegion());
                    GenomicSet<Gene> rightOverlappedGenes = genes.completeOverlap(cnv.getRightOverlappedDomainRegion());
                    Double leftOverlppedPhenogramScore = phenotypeData.phenoGramScore(cnv.getPhenotypes(), leftOverlappedGenes);
                    Double rightOverlappedPhenogramScore = phenotypeData.phenoGramScore(cnv.getPhenotypes(), rightOverlappedGenes);
                    boolean invGene = (hasLeftEnhancer && rightOverlappedPhenogramScore > 0) || (hasRightEnhancer && leftOverlppedPhenogramScore > 0);
                    if (!invEnhancer && !invGene) {
                        cnv.setEffectMechanism("InvEA", "noInvEA");
                    }
                    if (!invEnhancer && invGene) {
                        cnv.setEffectMechanism("InvEA", "GeneInvEA");
                    }
                    if (invEnhancer && !invGene) {
                        cnv.setEffectMechanism("InvEA", "EnhancerInvEA");
                    }
                    if (invEnhancer && invGene) {
                        cnv.setEffectMechanism("InvEA", "EandGInvEA");
                    }
                } else {
                    cnv.setEffectMechanism("InvEA", "noTDB");
                }
            } else {
                cnv.setEffectMechanism("InvEA", "NA");
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
    public static void tandemDuplicationEnhancerAdoption(GenomicSet<CNV> cnvs, GenomicSet<Gene> genes, GenomicSet<GenomicElement> enhancers, PhenotypeData phenotypeData) {
        for (CNV cnv : cnvs.values()) {
            if ("gain".equals(cnv.getType())) {
                GenomicSet<Gene> leftGenes = genes.completeOverlap(cnv.getLeftOverlappedDomainRegion());
                GenomicSet<Gene> rightGenes = genes.completeOverlap(cnv.getRightOverlappedDomainRegion());
                Double leftPhenogramScore = phenotypeData.phenoGramScore(cnv.getPhenotypes(), leftGenes);
                Double rightPhenogramScore = phenotypeData.phenoGramScore(cnv.getPhenotypes(), rightGenes);
                GenomicSet<GenomicElement> leftEnhancers = enhancers.completeOverlap(cnv.getLeftOverlappedDomainRegion());
                GenomicSet<GenomicElement> rightEnhancers = enhancers.completeOverlap(cnv.getRightOverlappedDomainRegion());
                boolean TanDupEA = (!leftEnhancers.isEmpty() && rightPhenogramScore > 0) || (leftPhenogramScore > 0 && !rightEnhancers.isEmpty());
                if (TanDupEA) {
                    cnv.setEffectMechanism("TanDupEA", "TanDupEA");
                } else {
                    if (cnv.getOverlapPhenogramScore() > 0) {
                        cnv.setEffectMechanism("TanDupEA", "onlyGDE");
                    } else {
                        cnv.setEffectMechanism("TanDupEA", "NoData");
                    }
                }
            } else {
                cnv.setEffectMechanism("TanDupEA", "NA");
            }
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
    public static void annotateTDBDjustByScore(GenomicSet<CNV> cnvs, PhenotypeData phenotypeData) {

        for (CNV cnv : cnvs.values()) {
        
            // initialize with False
            boolean isTDBD = false;
            
            if (cnv.hasBoundaryOverlap()) {
                
                // check if CNV has adjacent enhancers
                boolean hasLeftEnhancer = !cnv.getEnhancersInLeftRegion().isEmpty();
                boolean hasRightEnhancer = !cnv.getEnhancersInRightRegion().isEmpty();

                if (hasLeftEnhancer & cnv.getRightAdjacentPhenogramScore() > 0) {
                    isTDBD = true;
                    addInteractionChanges(
                            cnv, 
                            cnv.getEnhancersInLeftRegion(),
                            cnv.getGenesInRightRegion(), 
                            phenotypeData,
                            "gain", "newTDBD"
                    );
                    
                }

                if (hasRightEnhancer & cnv.getLeftAdjacentPhenogramScore() > 0) {
                    isTDBD = true;
                    addInteractionChanges(
                            cnv, 
                            cnv.getEnhancersInRightRegion(),
                            cnv.getGenesInLeftRegion(), 
                            phenotypeData,
                            "gain", "newTDBD"
                    );
                }
            }
            
            // annotate potential gene dosage effects
            addGeneDosageAnnotation(cnv, cnv.getGenesInOverlap(), phenotypeData, "newTDBD");
            
            // assigne effect mechanism class
            Double maxAdjacentScore = Math.max(cnv.getLeftAdjacentPhenogramScore(), cnv.getRightAdjacentPhenogramScore());
            if (isTDBD) {
                if (maxAdjacentScore > cnv.getOverlapPhenogramScore()) {
                    cnv.setEffectMechanism("newTDBD", "TDBD");
                } else {
                    cnv.setEffectMechanism("newTDBD", "Mixed");
                }
            } else {
                if (cnv.getOverlapPhenogramScore() > 0) {
                    cnv.setEffectMechanism("newTDBD", "GDE");
                } else {
                    cnv.setEffectMechanism("newTDBD", "NoData");
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
     * @param term2genes
     */
    public static void annotateEAlowG(GenomicSet<CNV> cnvs, HashMap<Term, HashSet<String>> term2genes) {
        for (CNV cnv : cnvs.values()) {
            HashSet<String> allTargetGenes = term2genes.get(cnv.getTargetTerm());
            if (matchingEnhancerAndGene(cnv, allTargetGenes)) {
                Double maxAdjacentScore = Math.max(cnv.getLeftAdjacentPhenogramScore(), cnv.getRightAdjacentPhenogramScore());
                if (maxAdjacentScore > cnv.getOverlapPhenogramScore() && cnv.getOverlapPhenogramScore() > 0) {
                    cnv.setEffectMechanism("EAlowG", "EAlowG");
                } else {
                    cnv.setEffectMechanism("EAlowG", "Mixed");
                }
            } else {
                if (cnv.getOverlapPhenogramScore() > 0) {
                    cnv.setEffectMechanism("EAlowG", "GDE");
                } else {
                    cnv.setEffectMechanism("EAlowG", "NoData");
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
    public static void annotateTDBD(GenomicSet<CNV> cnvs, HashMap<Term, HashSet<String>> term2genes) {
        
        for (CNV cnv : cnvs.values()) {
            
            // get all genes relevant to the phenotype observed in the patient
            HashSet<String> allTargetGenes = term2genes.get(cnv.getTargetTerm());
            
            // test for matching enahncer and gene signature adjacent to the CNV
            boolean isTDBD = cnv.hasBoundaryOverlap() && matchingEnhancerAndGene(cnv, allTargetGenes);

            if (isTDBD) {
                // get maximal score
                Double maxAdjacentScore = Math.max(
                        cnv.getLeftAdjacentPhenogramScore(), 
                        cnv.getRightAdjacentPhenogramScore()
                );
                if (maxAdjacentScore > cnv.getOverlapPhenogramScore()) {
                    cnv.setEffectMechanism("TDBD", "TDBD");
                } else {
                    cnv.setEffectMechanism("TDBD", "Mixed");
                }
            } else {
                if (cnv.getOverlapPhenogramScore() > 0) {
                    cnv.setEffectMechanism("TDBD", "GDE");
                } else {
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
     * @param term2genes
     */
    public static void annotateEA(GenomicSet<CNV> cnvs, HashMap<Term, HashSet<String>> term2genes) {
        for (CNV cnv : cnvs.values()) {
            HashSet<String> allTargetGenes = term2genes.get(cnv.getTargetTerm());
            if (matchingEnhancerAndGene(cnv, allTargetGenes)) {
                Double maxAdjacentScore = Math.max(cnv.getLeftAdjacentPhenogramScore(), cnv.getRightAdjacentPhenogramScore());
                if (maxAdjacentScore > cnv.getOverlapPhenogramScore()) {
                    cnv.setEffectMechanism("EA", "EA");
                } else {
                    cnv.setEffectMechanism("EA", "Mixed");
                }
            } else {
                if (cnv.getOverlapPhenogramScore() > 0) {
                    cnv.setEffectMechanism("EA", "GDE");
                } else {
                    cnv.setEffectMechanism("EA", "NoData");
                }
            }
        }
    }
    
    /**
     * add annotation of {@link InteractionChange}s to the input {@link CNV} for 
     * each input {@link Gene}.
     * 
     * @param cnv
     * @param enhancers
     * @param genes
     * @param phenotypeData
     * @param method
     * @param changeDirection 
     */
    private static void addInteractionChanges(CNV cnv, 
            GenomicSet<GenomicElement> enhancers, GenomicSet<Gene> genes, 
            PhenotypeData phenotypeData, 
            String changeDirection, String method) {
        
        // handle each input gene
        for (Gene gene: genes.values()){
            
            // get phenomatch score of the gene
            Double score = phenotypeData.phenoMatchScore(cnv.getPhenotypes(), gene);
 
            // if the gene is relevant to the phenotypes observed in the patient
            if (score > 0){
            
                // create an interactionChange annotation object
                InteractionChange intChange = new InteractionChange(
                        method, gene, score, enhancers, 
                        cnv.getBoundaryOverlap(), changeDirection
                );
                
                // add the interactionChange annotation to the effect annotation list
                cnv.addInteractionChangeAnnotaion(intChange);
            }
        }
        
    }
    
    /**
     * add annotation of {@link GeneDosage}s to the input {@link CNV} for 
     * each input {@link Gene}.
     * 
     * @param cnv
     * @param genes
     * @param phenotypedata
     * @param method 
     */
    private static void addGeneDosageAnnotation(CNV cnv, 
            GenomicSet<Gene> genes, PhenotypeData phenotypeData,
            String method){
        
        // handle each input gene
        for (Gene gene: genes.values()){
            
            // get phenomatch score of the gene
            Double score = phenotypeData.phenoMatchScore(cnv.getPhenotypes(), gene);
            
            // if the gene is relevant to the phenotypes observed in the patient
            if (score > 0){
                
                // get "loss" or "gain" from the 
                String dosageChange = cnv.getType();

                // create a new GeneDosage annotaiton object
                GeneDosage gde = new GeneDosage(method, gene, score, dosageChange);

                // add it to the CNV's annotaion list
                cnv.addGeneDosageAnnotaion(gde);
            }
        }
        
    }
    
//    /**
//     * From the list of {@link EffectAnnotation}s of the input {@link CNV}, 
//     * return the one that best explaines the phenotypes of the patient. 
//     * That is the effect mechanism with the highest score (pheontmatch score of the invoved gene).
//     * 
//     * @param cnv
//     * @return 
//     */
//    public EffectAnnotation getBestExplaination(CNV cnv){
//        
//        
//    }
    
}
