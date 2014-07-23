/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package phenotypeontology;

import genomicregions.Gene;
import genomicregions.GenomicSet;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
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
 * Implements several functionalities of a phenotype ontology and corresponding 
 * annotation tables such as information content access and similarity calculations
 * between tmers and sets of terms.
 * 
 * @author jonas
 */
public class OntologyWrapper  {
    
    /**
     * The phenotype ontology as {@link Ontology} object form the {@link ontologizer.go} package.
     */
    private Ontology ontology;
    
    /**
     * Mapping of gene IDs to phenotype terms.
     */
    private HashMap<String, HashSet<Term>>  gene2Terms;
    
    /**
     * Mapping from each {@link Term} to its information content (IC).
     */
    private HashMap<Term, Double> term2ic;
    
    /**
     * Similarity object from the {@link similarity.concepts.ResnikSimilarity} class
     * in the {@code de.sonumina.javautil} project.
     */
    private ResnikSimilarity sim;
    
    /**
     * Construct a new ontologyWrapper instance.
     * 
     * @param oboFilePath   
     *          path to the ontology file in .obo format
     * @param annotationFilePath 
     *          path to the file with phenotype annotation of genes
     */
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

        this.gene2Terms = readAnnotations(annotationFilePath);
        this.term2ic = getTerm2InformationContent(gene2Terms);

        this.sim = new ResnikSimilarity(getOntology(), term2ic);

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
        for (Term t_g : gene.getPhenotypeTerms()) {
            
            // initialize maximum over all patient's terms
            double bestGeneTermScore = 0;
            
            // iterate over all term t_p  with the patient
            for (Term t_p : terms) {
                
                // compute pairwise similarity 
                double termSim = this.sim.computeSimilarity(t_p, t_g);
                
                // take it as max if sim is larger
                if (termSim > bestGeneTermScore)
                    bestGeneTermScore = termSim;
            }
            
            // TODO: implement the additonal parameters lambda and k
            //if (bestGeneTermScore >= lambda) {
            //    similarity = similarity + Math.pow(bestGeneTermScore, k);
            //}
            
            similarity += bestGeneTermScore;
        }

        return similarity;
    }

    /**
     * Compute the phenogram score of a set of genes.
     * Here, this is just the maximum phenoMatch score over all input genes.
     * 
     * @param patientTerms  set of terms used to annotate the patient
     * @param genes set of genes that are compared to the phenotyes
     * @return the phenogram score for a set of genes
     */
    public double phenoGramScore(HashSet<Term> patientTerms, GenomicSet<Gene> genes ){
        
        // if gene set is empty, return zero
        if (genes.isEmpty()){
            return 0.0;
        }else{
            // initialize list for all phenoMatch scores
            ArrayList<Double> phenoMatchScores = new ArrayList(genes.size());

            // compute for each gene the phenoMatch score
            for (Gene g: genes.values()){
                phenoMatchScores.add(phenoMatchScore(patientTerms, g));
            }

            // return the maximal phenoMatch score
            return Collections.max(phenoMatchScores);   
        }
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

                            Term t = getOntology().getTermIncludingAlternatives(annotatedTermId);
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
            SlimDirectedGraphView<Term> ontologySlim = getOntology().getSlimGraphView();
            
            HashMap<Term, HashSet<String>> term2geneids = new HashMap<Term, HashSet<String>>();
            HashMap<String, HashSet<Term>> gene2terms = geneId2annotations;

            for (String geneId : gene2terms.keySet()) {

                    HashSet<Term> termsOfGene = gene2terms.get(geneId);
                    for (Term annotated : termsOfGene) {

                            // System.out.println("annot : "+annotated);
                            if (!ontologySlim.vertex2Index.containsKey(annotated)) {
                                    annotated = getOntology().getTermIncludingAlternatives(annotated.getIDAsString());
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

            Term root = getOntology().getRootTerm();
            HashMap<Term, Integer> term2frequency = new HashMap<Term, Integer>();
            for (Term t : term2objectIds.keySet()) {
                    term2frequency.put(t, term2objectIds.get(t).size());
            }
            int maxFreq = term2frequency.get(root);
            HashMap<Term, Double> term2informationContent = SimilarityUtilities.caculateInformationContent(maxFreq, term2frequency);

            int frequencyZeroCounter = 0;
            double ICzeroCountTerms = -1 * (Math.log(1 / (double) maxFreq));

            for (Term t : getOntology()) {
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

    /**
     * The phenotype ontology as {@link Ontology} object form the {@link ontologizer.go} package.
     * @return the ontology
     */
    public Ontology getOntology() {
        return ontology;
    }
    
    /**
     * Returns the term to a given term ID string or null.
     * 
     * @param term
     * @return 
     */
    public Term getTerm(String term){
        return this.getOntology().getTerm(term);
    }    

    /**
     * Tests if input gene ID is contained. If so, there might be phenotype
     * associations for this gene.
     * 
     * @param gID gene ID to search for
     * @return ture if the ID is found else false
     */
    public boolean containsGene(String gID) {
        return this.gene2Terms.containsKey(gID);
    }

    /**
     * Returns a {@link HashSet} of phenotype {@link Term}s to which the input 
     * gene is associated.
     * 
     * @param gID ID of the gene for which the phenotypes should be retrieved
     * @return {@link HashSet} of phenotype {@link Term}s to which the input 
     * gene is associated.
     */
    public HashSet<Term> getGenePhenotypes(String gID) {
        return this.gene2Terms.get(gID);
    }

    /**
     * get all gene IDs of genes that are mapped to phenotype annotations.
     * 
     * @return {@link HashSet} of all gene IDs that are mapped to phenotype 
     * annotations
     */
    public Set<String> getAllGenesIDs() {
        return this.gene2Terms.keySet();
    }

    /**
     * Returns the information content (IC) for an input {@link Term}.
     * The information content is calculated as the negative logarithm of the
     * frequency p_t of annotations to term t.
     * <br><br>
     * <center>
     * IC(t) = -log(p_t)
     * </center>
     * 
     * @param term
     * @return 
     */
    public Double getIC(Term term) {
        return this.term2ic.get(term);
    }

    /**
     * Returns an iterator to iterate over all terms.
     * @return an iterator
     */
    public Iterator<Term> iterator() {
        return this.ontology.iterator();
    }

}
    
