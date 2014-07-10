/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package phenotypeontology;

import genomicregions.Gene;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import ontologizer.go.Term;
import java.util.HashSet;
import ontologizer.go.OBOParser;
import ontologizer.go.OBOParserException;
import ontologizer.go.Ontology;
import ontologizer.go.TermContainer;

import similarity.SimilarityUtilities;
import similarity.concepts.ResnikSimilarity;
import similarity.objects.InformationContentObjectSimilarity;
import sonumina.math.graph.SlimDirectedGraphView;


/**
 *
 * @author jonas
 */
public class phenogramScore {

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
