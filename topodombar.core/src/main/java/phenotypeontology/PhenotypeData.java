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
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import ontologizer.go.OBOParser;
import ontologizer.go.OBOParserException;
import ontologizer.go.OBOParserFileInput;
import ontologizer.go.Ontology;
import ontologizer.go.Term;
import ontologizer.go.TermContainer;
import ontologizer.go.TermID;
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
public class PhenotypeData  implements Cloneable{
    
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
    private final HashMap<Term, Double> term2ic;
    
    /**
     * Similarity object from the {@link similarity.concepts.ResnikSimilarity} class
     * in the {@code de.sonumina.javautil} project.
     */
    public final ResnikSimilarity sim;
    
    /** maps all phenotype terms to its more general ancestor terms */
    private final HashMap<Term, HashSet<Term>> term2ancestors;
    
    /** maps all phenotype terms to its more specific descendant terms */
    private final HashMap<Term, HashSet<Term>> term2descendants;
    
    /**
     * List of terms in topological order
     */
    private final ArrayList<TermID> orderedTermIDs;
    
    /**
     * Constructor with all members as variables.
     * 
     * @param ontology
     * @param gene2Terms
     * @param term2ic
     * @param sim
     * @param term2ancestors
     * @param term2descendants 
     */
    public PhenotypeData(Ontology ontology, 
            HashMap<String, HashSet<Term>>  gene2Terms,
            HashMap<Term, Double> term2ic,
            ResnikSimilarity sim,
            HashMap<Term, HashSet<Term>> term2ancestors,
            HashMap<Term, HashSet<Term>> term2descendants){
        
        this.ontology = ontology;
        this.gene2Terms = gene2Terms;
        this.term2ic = term2ic;
        this.sim = sim;
        this.term2ancestors = term2ancestors;
        this.term2descendants = term2descendants;
        ArrayList<Term> orderedTerms = ontology.getTermsInTopologicalOrder();
        this.orderedTermIDs = new ArrayList<TermID>();
        for (Term t: orderedTerms){
            this.orderedTermIDs.add(t.getID());
        }
    }
    
    /**
     * Construct a new phenotypeData instance.
     * 
     * @param oboFilePath   
     *          path to the ontology file in .obo format
     * @param annotationFilePath 
     *          path to the file with phenotype annotation of genes
     */
    public PhenotypeData(String oboFilePath, String annotationFilePath) throws IOException{
        
        // build oboParser object for the input ontology file
        OBOParser oboParser = new OBOParser(new OBOParserFileInput(oboFilePath));
        
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
        
        SlimDirectedGraphView<Term> ontologySlim = ontology.getSlimGraphView();

        this.gene2Terms = readAnnotations(annotationFilePath);
        this.term2ic = getTerm2InformationContent(gene2Terms);

        this.sim = new ResnikSimilarity(getOntology(), term2ic);

        // build mapping from each term to ancestors
        term2ancestors = new HashMap<>();
        term2descendants = new HashMap<>();
        for(Term t : ontology){
            
            // get ancestors from ontologySlim object
            term2ancestors.put(t, new HashSet<Term>(ontologySlim.getAncestors(t)));

            // get descendants form ontologyslim object
            term2descendants.put(t, new HashSet<Term>(ontologySlim.getDescendants(t)));

        }
        ArrayList<Term> orderedTerms = ontology.getTermsInTopologicalOrder();
        this.orderedTermIDs = new ArrayList<TermID>();
        for (Term t: orderedTerms){
            this.orderedTermIDs.add(t.getID());
        }
    }
    
    /**
     * returns a shallow copy of this {@link PhenotypeData} objecekt.
     * @return 
     */
    public PhenotypeData shallowCopy(){
        
        // make a new PhenotypeData object
        PhenotypeData copy = new PhenotypeData(
                this.ontology, this.gene2Terms, this.term2ic, this.sim, 
                this.term2ancestors, this.term2descendants);
        return copy;
    }
    
    /**
     * Coputes the Resnik Similarity score of two input terms.
     * @param t1
     * @param t2
     * @return 
     */
    public double resnikSim(Term t1, Term t2){
        
        // get all shared parent terms
        ArrayList<TermID> par = (ArrayList<TermID>) ontology.getSharedParents(t1.getID(), t2.getID());
                
        // initialize maximum score
        Double maxScore = 0.0;
        
        if (par.size() > 0){
            // take the first term, since it is the first common parent in the path from t2 to root.
            maxScore = getIC(ontology.getTerm(par.get(0)));
        }
        return(maxScore);
    }
    
    /**
     * Coputes the Resnik Similarity score of two input terms and returns the score together with the term id.
     * @param t1
     * @param t2
     * @return {@link TermPair} object with matching score and responsibel common parent term. 
     */
    public TermPair resnikSimWithTerm(Term t1, Term t2){
        
        // get all shared parent terms
        ArrayList<TermID> par = (ArrayList<TermID>) ontology.getSharedParents(t1.getID(), t2.getID());
        
        // initialize maximum term and score
        Term maxTerm = null;
        Double maxScore = 0.0;
        
        if (par.size() > 0){
            // take first term in order (with lowest term in ontology level, e.g. highest resnik score)
            maxTerm = this.ontology.getTerm(par.get(0));
            maxScore = getIC(maxTerm);
        }
        // construct a new pair of the maximum value with its term to retunr both
        TermPair tp = new TermPair(t1, t2, maxScore, maxTerm);
        
        return(tp);
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
//                double termSim = resnikSim(t_p, t_g);

//                TermPair tp = resnikSimWithTerm(t_p, t_g);
//                double termSim = tp.getS();
                
                
                // take it as max if sim is larger
                if (termSim > bestGeneTermScore)
                    bestGeneTermScore = termSim;
            }
            
            // TODO: implement the additonal parameters lambda and k
            //if (bestGeneTermScore >= lambda) {
            //    similarity = similarity + Math.pow(bestGeneTermScore, k);
            //}
        
            // take max across gene terms
//            System.out.println("DEBUG: return max across gene terms!");
            //similarity += bestGeneTermScore;
            if (bestGeneTermScore >= similarity){
                similarity = bestGeneTermScore;
            }
        }

        return similarity;
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
    public ArrayList<TermPair> phenoMatchScoreWithMatching(HashSet<Term> terms, Gene gene){
        
        // initialize matching
        ArrayList<TermPair> matching = new ArrayList<>();
        
        // iterate over all terms  with the input gene
        for (Term t_g : gene.getPhenotypeTerms()) {
            
            ArrayList<TermPair> geneTermPairs = new ArrayList<>();
            
            // iterate over all term t_p  with the patient
            for (Term t_p : terms) {
                
                // compute pairwise similarity and get responsible term                
                TermPair tp = resnikSimWithTerm(t_p, t_g);
                
                // add term pair to list pair to all patient terms
                geneTermPairs.add(tp);
                    
            }
            
            // get maxium over all patient terms for this gene
            TermPair maxPair = Collections.max(geneTermPairs, TermPair.TERM_PAIR_SCORE_ORDER);
            
            // check if maxPair has score larger than zeoro
            // check that lowest common ancester is not the root
            if (maxPair.getS() > 0.0){
                // add pair with highest score to output matching list
                matching.add(maxPair);
            }
        }
        
//        // take max across gene terms
//        ArrayList<TermPair> maxMatching = new ArrayList<TermPair>();
//        if (matching.size() > 0){
//            maxMatching.add(Collections.max(matching, TermPair.TERM_PAIR_SCORE_ORDER));
//        }
//        return maxMatching;
        return matching;
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
            //TODO: with double >...
            // compute for each gene the phenoMatch score
            for (Gene g: genes.values()){
                phenoMatchScores.add(phenoMatchScore(patientTerms, g));
            }

            // return the maximal phenoMatch score
            return Collections.max(phenoMatchScores);   
        }
    }
        
    /**
     * builds a mapping for phenotpye terms to the set of genes that are associated
     * with the phenotype or more specific descendants of the phenotype.
     * 
     * @param targetTerms the input terms for which the mapping should be computed.
     * @return a mapping for each input term to a set of genes associated to the term.
     */
    public HashMap<Term, HashSet<String>> mapTargetTermToGenes(ArrayList<TargetTerm> targetTerms){
        
        HashMap<Term, HashSet<String>> term2genes = new HashMap();
        // for all target terms initialize empty set:
        for (TargetTerm tT: targetTerms){
            term2genes.put(tT.getTerm(), new HashSet<String>());
        }
        
        // iterate over all genes that have phenotype associations
        for (String geneID : this.gene2Terms.keySet()){
            
            HashSet<Term> geneAncestorSet = new HashSet<Term>();
            
            // iterate over all phenotype terms associated with this gene:
            for (Term t : this.gene2Terms.get(geneID)){
                
                // test if gene term t is an ancester of target terms:
                for (TargetTerm tT: targetTerms){
                    
                    // check if targetTerm tT equals gene term t
                    // or if tT is an ancesotr of t.
                    if ( this.isAncestorOrEqual(t, tT.getTerm()) ){

                        // add gene to set of targetTerm associated genes
                        term2genes.get(tT.getTerm()).add(geneID);

                    }
                }
            }
        }                
        return term2genes;
    }

    /**
     * Tests if tow input terms are equal or have an ancestor/descendant relation.
     * @param term  input term tested to be descendant
     * @param ancestor input term tested to be ancestor 
     * @return ture if {@code ancestor} is equal to or an ancestor of {@code term}. 
     */
    public boolean isAncestorOrEqual(Term term, Term ancestor){
        
        return this.term2ancestors.get(term).contains(ancestor);
    }
    
    /**
     * Reads the mapping form genes to set of phenotypes. Genes associated to phenotypes 
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
     * TODO: use term2ancestor mapping, not slimOntology
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
     * returns all terms in the ontology
     * @return all terms
     */
    public HashSet<Term> getAllTerms(){
        HashSet<Term> terms = new HashSet<Term>();
        for (Term t: this.ontology){
            terms.add(t);
        }
        return terms;
    }
    
    /**
     * Returns the term to a given term ID string or null.
     * 
     * @param term
     * @return 
     */
    public Term getTermIncludingAlternatives(String term) throws IOException{
        Term t = this.getOntology().getTermIncludingAlternatives(term);
        if (t==null){
            System.out.println("WARNING term is null:");
            System.out.println(term);
            throw new IOException(String.format(
                        "[ERROR] this term '%s'. is null. Probably it "
                                + "was not found in the provied phenotype "
                                + "ontology.", term));
        }
        return t;
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

    /**
     * returns a set of all terms that are ancestors of the input term t.
     * @param t input term
     * @return all ancestors of t (including t)
     */
    public HashSet<Term> getAncestors(Term t) {
        return this.term2ancestors.get(t);
    }
    
    /**
     * returns a set of all terms that are descendants of the input term t.
     * @param t input term
     * @return all descendants of t (including t itself)
     */
    public HashSet<Term> getDescendants(Term t){
        return this.term2descendants.get(t);
    }

    /**
     * Mapping of gene IDs to phenotype terms.
     * @param gene2Terms the gene2Terms to set
     */
    public void setGene2Terms(HashMap<String, HashSet<Term>> gene2Terms) {
        this.gene2Terms = gene2Terms;
    }

}
    
