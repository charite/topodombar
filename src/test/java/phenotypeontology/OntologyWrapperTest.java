/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package phenotypeontology;

import genomicregions.Gene;
import genomicregions.GenomicSet;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import ontologizer.go.Term;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Unit tests for the OntologyWrapper class.
 * See description of the example data in the test class {@link TopodombarTest}.
 * 
 <pre>
 {@code 
    Toy example dataset for testing:
    ================================
                        10        20        30        40    
    Coord:    01234567890123456789012345678901234567890
    Domain    <---------->   <-------------------->        (0,12), (15,37)
    Boundary              ###                              (12,15)
    GeneC,B,D,A -->     --->    -->     ----->             C:(2,5), B:(10,14), D:(18,20), A:(26,32)
    Enhacer     **                      **                 (2,4), (26,28) 
    cnv1               ==========                          (9,19)
    cnv2              =========================            (8,33)
    cnv3                  =======                          (12,19)
    cnv4                    ==                             (14,16)
                        10        20        30        40    
              01234567890123456789012345678901234567890


    Example phenotype ontology (EPO) data for testing:

                        EP00
                       /    \       
                    EP01     EP02        
                   /       /    \       
                EP03    EP04    EP05
                           \    /  \
                            EP06    EP07

    Term    p   IC=-log(p)
    ========================
    EP01    1   0
    EP02    .25 1.39
    EP03    .75 0.29
    EP04    .25 1.39
    EP05    .75 0.29
    EP06    0   1.39 (Inf)
    EP07    .25 1.39

    CNV phenotypes: EP06
    CNV targetTerm: EP06

    Gene    Annotation  PhenoMatchScore (to EP:06)
    ==================================
    GeneA   EP04,EP05   1.68
    GeneB   EP07        0.29
    GeneC   EP03        0
    GeneD   EP05        0.29
  }
 </pre>
 * 
 * 
 * @author jonas
 */
public class OntologyWrapperTest {
    
    private static OntologyWrapper ontologyWrapper;
    
    public OntologyWrapperTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {

        // create ontlolgyWrapper object form example data:
        java.net.URL oboURL = OntologyWrapperTest.class.getResource("/example_ontology.obo");
        String oboPath = oboURL.getPath();
        
        String annotPath = OntologyWrapperTest.class.getResource("/example_genes_to_penotype.txt").getPath();
        
        // parse ontology and create wrapper object
        ontologyWrapper = new OntologyWrapper(oboPath, annotPath);        
    }
    
    @Test
    public void testOntologyParsing() {

        // test if all terms are read
        System.out.println("TEST: Numer of terms in example data");
        assertEquals(8, ontologyWrapper.getOntology().getNumberOfTerms());
        
        // test the ic of all terms
        // print it to output
        System.out.println("TEST: print IC of all Terms in example dataset:");
        for (Iterator<Term> it = ontologyWrapper.iterator(); it.hasNext(); ){
            Term t = it.next();
            System.out.println("TEST: term and IC:" + t.toString() + ontologyWrapper.getIC(t));
        }
        
    }
    
    @Test
    public void testInformationContentCalculation() {
        // test IC calculation of example term EP:01 with frequency p=.25
        Term t1 = ontologyWrapper.getTerm("EP:01");
        Double expIC = -Math.log(0.25);
        
        System.out.println("TEST: IC of EP:01");
        assertEquals(expIC, ontologyWrapper.getIC(t1));

        // test i all genes are read:
        System.out.println("TEST: Number of genes with annotation in example dataset.");
        assertEquals(4, ontologyWrapper.getAllGenesIDs().size());
        
        // test if geneA is containted in term2ic dict
        HashSet<String> genes = new HashSet<String>();
        genes.add("geneA");
        genes.add("geneB");
        genes.add("geneC");
        genes.add("geneD");
        System.out.println("TEST: Test all gene IDs in example dataset.");
        assertEquals(genes, ontologyWrapper.getAllGenesIDs());
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of phenoMatchScore method, of class OntologyWrapper.
     */
    @Test
    public void testPhenoMatchScore() {
        System.out.println("phenoMatchScore");
        HashSet<Term> terms = new HashSet<Term>();
        terms.add(ontologyWrapper.getTerm("EP:06"));
        
        // build gene A of example data set
        ArrayList<String> genePhenotypes = new ArrayList<String>();
        genePhenotypes.add("EP:04");
        genePhenotypes.add("EP:05");        
        Gene geneA = new Gene("chr1", 26, 32, "geneA", genePhenotypes);
        geneA.setPhenotypeTerms( new HashSet<Term>() );
        geneA.addPhenotypeTerm(ontologyWrapper.getTerm("EP:04"));
        geneA.addPhenotypeTerm(ontologyWrapper.getTerm("EP:05"));
        
        // expect sum IC of most specific common terms for each gene phenotype
        double expResult = ontologyWrapper.getIC(ontologyWrapper.getTerm("EP:04"));
        expResult += ontologyWrapper.getIC(ontologyWrapper.getTerm("EP:05"));
        double result = ontologyWrapper.phenoMatchScore(terms, geneA);
        
        assertEquals(expResult, result, 0.001);
        
    }

    /**
     * Test of phenoGramScore method, of class OntologyWrapper.
     */
    @Test
    public void testPhenoGramScore() {
        System.out.println("phenoGramScore");
        
        // build set of patient terms with only EP:06 from the example dataset
        HashSet<Term> patientTerms = new HashSet<Term>();
        patientTerms.add(ontologyWrapper.getTerm("EP:06"));
        
        // build gene A of example data set
        ArrayList<String> genePhenotypes = new ArrayList<String>();
        genePhenotypes.add("EP:04");
        genePhenotypes.add("EP:05");        
        Gene geneA = new Gene("chr1", 26, 32, "geneA", genePhenotypes);
        geneA.setPhenotypeTerms( new HashSet<Term>() );
        geneA.addPhenotypeTerm(ontologyWrapper.getTerm("EP:04"));
        geneA.addPhenotypeTerm(ontologyWrapper.getTerm("EP:05"));
        
        // build gene B of example data set
        ArrayList<String> geneBPhenotypes = new ArrayList<String>();
        geneBPhenotypes.add("EP:07");
        Gene geneB = new Gene("chr1", 26, 32, "geneB", geneBPhenotypes);
        geneB.setPhenotypeTerms( new HashSet<Term>() );
        geneB.addPhenotypeTerm(ontologyWrapper.getTerm("EP:07"));

        // add geneA and geneB to a genomic set of genes
        GenomicSet<Gene> genes = new GenomicSet<Gene>();
        genes.put("geneA", geneA);
        genes.put("geneB", geneB);
        
        // since geneA have higher phenomatchScore than geneB, we expect the
        // phenomatch score of geneB here again.
        double expResult = ontologyWrapper.getIC(ontologyWrapper.getTerm("EP:04"));
        expResult += ontologyWrapper.getIC(ontologyWrapper.getTerm("EP:05"));

        // calculate phenoMatch score for the patient temrs (EP:06) and the gene phenotypes
        double result = ontologyWrapper.phenoGramScore(patientTerms, genes);

        // assert eauality with a tolerance of 0.001
        assertEquals(expResult, result, 0.001);
        
    }
    
}
