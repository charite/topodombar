/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package phenotypeontology;

import genomicregions.Gene;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import ontologizer.go.OBOParser;
import ontologizer.go.OBOParserException;
import ontologizer.go.Ontology;
import ontologizer.go.Term;
import ontologizer.go.TermContainer;
import similarity.SimilarityUtilities;
import similarity.concepts.ResnikSimilarity;
import sonumina.math.graph.SlimDirectedGraphView;

/**
 *
 * @author jonas
 */
public class OntologyWrapper {
    
    public Ontology ontology;
    public HashMap<Term, Double> term2ic;
    public ResnikSimilarity sim;
    
    public OntologyWrapper(String oboFilePath, String annotationFilePath){
        
        // build oboParser object for the input ontology file
        OBOParser oboParser = new OBOParser(oboFilePath);
        
        // parse obo file
        try {
            String parseInfo = oboParser.doParse();
            System.out.println(parseInfo);
        } catch (OBOParserException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        // get the complete hpo
        TermContainer termContainer = new TermContainer(oboParser.getTermMap(), oboParser.getFormatVersion(), oboParser.getDate());
        this.ontology = new Ontology(termContainer);
        
        //SlimDirectedGraphView<Term> ontologySlim = ontology.getSlimGraphView();

        HashMap<String, HashSet<Term>> geneId2annotations = readAnnotations(annotationFilePath);
        this.term2ic = getTerm2InformationContent(geneId2annotations);

        this.sim = new ResnikSimilarity(ontology, term2ic);

    }
    
    
    /**
     * Computes the pheno match score between a gene and a set of phenotype terms.
     * The phenomatch score is a similarity between phenotypes that are associated
     * with a single gene and another set of phenotypes, which might be symptoms of a patient
     * that have genetic variations associated with the gene.
     * The calculation is described in Ibn-Salem et al. (2014) Genome Biology.
     * 
     * @param terms a set of phenotype terms
     * @param gene a {@link Gene} object
     * @return phenomatch score
     */
    public double phenoMatchScore(HashSet<Term> terms, Gene gene){

        // initialize similarity to zero
        double similarity = 0;

        // DMM-paper forumula
        // iterate over all terms  with the input gene
        for (Term t_g : gene.phenotypeTerms) {
            
            // initialize maximum over all patient's terms
            double bestGeneTermScore = 0;
            
            // iterate over all term t_p  with the patient
            for (Term t_p : terms) {
                
                // compute pairwise similarity 
                double termSim = sim.computeSimilarity(t_p, t_g);
                
                // take it as max if sim is larger
                if (termSim > bestGeneTermScore)
                    bestGeneTermScore = termSim;
            }

            //if (bestGeneTermScore >= lambda) {
            //    similarity = similarity + Math.pow(bestGeneTermScore, k);
            //}
            similarity += bestGeneTermScore;
        }

        return similarity;
    }

    /**
     * Reads the mapping form genes to set of phenotypes. Genes associated to phenotpyes 
     * are parsed by reading the annotation tables provided by the HPO or 
     * Uberpheno project to build a mapping from genes to terms.
     * @author adopted form Sebastian Koehler
     * 
     * @param annotationFilePath path to the annotation file that maps genes to 
     * HPO Terms
     * @param ontology an {@link Ontology} object
     * @return a mapping from EntrezGene IDs to a set of corresponding HPO {@link Term}s.
     */
    private HashMap<String, HashSet<Term>> readAnnotations(String annotationFilePath) {
	
            // define some constant patterns
            Pattern semicolon = Pattern.compile(";");
            Pattern tabstopp = Pattern.compile("\t");
            Pattern uphenoTermidPattern = Pattern.compile("\\(([HZM]P:\\d{7})\\)$");

            HashMap<String, HashSet<Term>> gene2annot = new HashMap<String, HashSet<Term>>();
            try {

                    BufferedReader in = new BufferedReader(new FileReader(annotationFilePath));

                    String firstLine = in.readLine();
                    boolean isHpo = false;
                    boolean isUpheno = false;
                    if (firstLine.startsWith("#Entrez Gene ID of human gene ; Gene symbol ; Annotated Uberpheno")) {
                            isUpheno = true;
                    } else if (firstLine.startsWith("#Format: entrez-gene-id<tab>entrez-gene-symbol<tab>HPO-Term-Name<tab>HPO-Ter")) {
                            isHpo = true;
                    } else {
                            throw new RuntimeException("Can't handle annotation-file format!");
                    }

                    String line = null;
                    while ((line = in.readLine()) != null) {

                            if (line.startsWith("#"))
                                    continue;

                            String[] split = null;
                            String entrezId = null;
                            String annotatedTermId = null;
                            if (isUpheno) {
                                    split = semicolon.split(line);
                                    entrezId = split[0];

                                    Matcher matcher = uphenoTermidPattern.matcher(split[2]);
                                    if (!matcher.find()) {
                                            System.out.println("pattern not matched in: " + line);
                                    }
                                    annotatedTermId = matcher.group(1);

                            }
                            if (isHpo) {
                                    split = tabstopp.split(line);
                                    entrezId = split[0];
                                    annotatedTermId = split[3];
                            }

                            Term t = ontology.getTermIncludingAlternatives(annotatedTermId);
                            if (t == null) {
                                    System.err.println("Could not find term for ID:" + annotatedTermId + " parsed from line: " + line);
                                    continue;
                            }

                            HashSet<Term> annotationsOfGene;
                            if (gene2annot.containsKey(entrezId))
                                    annotationsOfGene = gene2annot.get(entrezId);
                            else
                                    annotationsOfGene = new HashSet<Term>();
                            annotationsOfGene.add(t);
                            gene2annot.put(entrezId, annotationsOfGene);

                    }
                    in.close();

            } catch (IOException e) {
                    e.printStackTrace();
            }

            return gene2annot;
    }

    /**
     * Builds a mapping from each {@Term} to its information content. Thereby the
     * information content is calculated as the negative logarithm of the frequency
     * of gene annotations to that term.
     * 
     * @param geneId2annotations    a mapping of genes to set of {@Term}s to which
     * they are annotated
     * @return mapping from term to information content
     */
    private HashMap<Term, Double> getTerm2InformationContent(HashMap<String, HashSet<Term>> geneId2annotations) {
            
            // Build a SlimDirectedGraphView object from the input ontology
            SlimDirectedGraphView<Term> ontologySlim = ontology.getSlimGraphView();
            
            HashMap<Term, HashSet<String>> term2geneids = new HashMap<Term, HashSet<String>>();
            HashMap<String, HashSet<Term>> gene2terms = geneId2annotations;

            for (String geneId : gene2terms.keySet()) {

                    HashSet<Term> termsOfGene = gene2terms.get(geneId);
                    for (Term annotated : termsOfGene) {

                            // System.out.println("annot : "+annotated);
                            if (!ontologySlim.vertex2Index.containsKey(annotated)) {
                                    annotated = ontology.getTermIncludingAlternatives(annotated.getIDAsString());
                                    // System.out.println(" now annot : "+annotated);
                            }

                            for (Term annotatedAndInduced : ontologySlim.getAncestors(annotated)) {

                                    HashSet<String> genesForTerm;
                                    if (term2geneids.containsKey(annotatedAndInduced))
                                            genesForTerm = term2geneids.get(annotatedAndInduced);
                                    else
                                            genesForTerm = new HashSet<String>();

                                    genesForTerm.add(geneId);
                                    term2geneids.put(annotatedAndInduced, genesForTerm);
                            }
                    }
            }
            return calculateTermIC(term2geneids);
    }
    
    private HashMap<Term, Double> calculateTermIC(HashMap<Term, HashSet<String>> term2objectIds) {

            Term root = ontology.getRootTerm();
            HashMap<Term, Integer> term2frequency = new HashMap<Term, Integer>();
            for (Term t : term2objectIds.keySet()) {
                    term2frequency.put(t, term2objectIds.get(t).size());
            }
            int maxFreq = term2frequency.get(root);
            HashMap<Term, Double> term2informationContent = SimilarityUtilities.caculateInformationContent(maxFreq, term2frequency);

            int frequencyZeroCounter = 0;
            double ICzeroCountTerms = -1 * (Math.log(1 / (double) maxFreq));

            for (Term t : ontology) {
                    if (!term2frequency.containsKey(t)) {
                            ++frequencyZeroCounter;
                            term2informationContent.put(t, ICzeroCountTerms);
                    }
            }

            System.out.println("WARN: Frequency of " + frequencyZeroCounter
                            + " terms was zero!! Calculated by -1 * (Math.log(1/(double)maxFreq)) = -1 * (Math.log(1/(double)" + maxFreq + ")))");
            System.out.println("Set IC of these to : " + ICzeroCountTerms);
            return term2informationContent;
    }

}
    
