/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package phenotypeontology;

import genomicregions.Gene;
import ontologizer.go.Term;
import java.util.HashSet;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Compute the phenogram score of a set of genes.
 * This is just the maximum phenoMatch score over all input genes.
 * 
 * @author jonas
 */
public class phenogramScore {

//    /**
//     * Compute the phenogram score of a set of genes.
//     * This is just the maximum phenoMatch score over all input genes.
//     * @param ontology ontologyWrapper instance
//     * @param patientTerms  set of terms used to annotate the patient
//     * @param genes set of genes
//     * @return 
//     */
//    public static double phenoGramScore(OntologyWrapper ontology, HashSet<Term> patientTerms, HashSet<Gene> genes ){
//                
//        ArrayList<Double> phenoMatchScores = new ArrayList<Double>(genes.size());
//        
//        for (Gene g: genes){
//            phenoMatchScores.add(ontology.phenoMatchScore(patientTerms, g));
//        }
//        
//        
//        return Collections.max(phenoMatchScores);
//    }

//     
//    public double phenoMatchScore(HashSet<Term> terms, Gene gene, double lambda, int k){
//        
//        
//
//        
//        //HashSet<Term> annot = geneId2annotations.get(geneId);
//        //ArrayList<Term> annotAl = new ArrayList<Term>(annot);
//        // annotAl ~> geneTerms 
//        
//        // get ArrayList of all Terms associated with the input gene
//        ArrayList<Term> geneTerms = new ArrayList<Term>(gene.phenotypeTerms);
//        ArrayList<Term> patientTerms = new ArrayList<Term>(terms);
//        
//        double similarity = 0;
//
//        // DMM-paper forumula
//        // iterate over all terms  with the input gene
//        for (Term t_g : geneTerms) {
//            
//            // initialize maximum over all patient's terms
//            double bestGeneTermScore = 0;
//            
//            // iterate over all term t_p  with the patient
//            for (Term t_p : patientTerms) {
//                
//                // compute pairwise similarity 
//                double termSim = sim.computeSimilarity(t_p, t_g);
//                
//                // take it as max if sim is larger
//                if (termSim > bestGeneTermScore)
//                    bestGeneTermScore = termSim;
//            }
//
//            if (bestGeneTermScore >= lambda) {
//                similarity = similarity + Math.pow(bestGeneTermScore, k);
//            }
//        }
//
//        return similarity;
//    }
    
}
