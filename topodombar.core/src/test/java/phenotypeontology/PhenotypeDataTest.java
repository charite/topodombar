/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package phenotypeontology;

import genomicregions.Gene;
import genomicregions.GenomicSet;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import ontologizer.go.Ontology;
import ontologizer.go.Term;
import ontologizer.go.TermRelation;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import toyexampledata.ExampleData;

/**
 * Unit tests for the PhenotypeData class.
 * See description of the toy example data in the class {@link ExampleData}.
 * 
 * @see ExampleData
 * 
 * @author jonas
 */
public class PhenotypeDataTest {
    
    private static PhenotypeData phenotypeData;
    
    private static ExampleData exampleData;
    
    public PhenotypeDataTest() {
    }
    
    @BeforeClass
    public static void setUpClass() throws IOException {

        // create PhenotypeData object form example data:
        java.net.URL oboURL = PhenotypeDataTest.class.getResource("/example_ontology.obo");
        String oboPath = oboURL.getPath();
        
        String annotPath = PhenotypeDataTest.class.getResource("/example_genes_to_penotype.txt").getPath();
        
        // parse ontology and create PhenotypeData object
        phenotypeData = new PhenotypeData(oboPath, annotPath);        
        
        exampleData = new ExampleData();
    }
    
    @Test
    public void testOntologyParsing() {

        // test if all terms are read
        System.out.println("TEST: Numer of terms in example data");
        assertEquals(8, phenotypeData.getOntology().getNumberOfTerms());
        
        // test the ic of all terms
        // print it to output
        System.out.println("TEST: print IC of all Terms in example dataset:");
        for (Iterator<Term> it = phenotypeData.iterator(); it.hasNext(); ){
            Term t = it.next();
            System.out.println("TEST: term and IC:" + t.toString() + phenotypeData.getIC(t));
        }
        
    }
    
    @Test
    public void testInformationContentCalculation() {
        // test IC calculation of example term EP:01 with frequency p=.25
        Term t1 = phenotypeData.getTermIncludingAlternatives("EP:01");
        Double expIC = -Math.log(0.25);
        
        System.out.println("TEST: IC of EP:01");
        assertEquals(expIC, phenotypeData.getIC(t1));

        // test i all genes are read:
        System.out.println("TEST: Number of genes with annotation in example dataset.");
        assertEquals(4, phenotypeData.getAllGenesIDs().size());
        
        // test if geneA is containted in term2ic dict
        HashSet<String> genes = new HashSet<String>();
        genes.add("geneA");
        genes.add("geneB");
        genes.add("geneC");
        genes.add("geneD");
        System.out.println("TEST: Test all gene IDs in example dataset.");
        assertEquals(genes, phenotypeData.getAllGenesIDs());
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
     * Test of phenoMatchScore method, of class PhenotypeData.
     */
    @Test
    public void testPhenoMatchScore() {
        System.out.println("phenoMatchScore");
        HashSet<Term> terms = new HashSet<Term>();
        terms.add(phenotypeData.getTermIncludingAlternatives("EP:06"));
        
        // build gene A of example data set
        ArrayList<String> genePhenotypes = new ArrayList<String>();
        genePhenotypes.add("EP:04");
        genePhenotypes.add("EP:05");        
        Gene geneA = new Gene("chr1", 26, 32, "geneA", genePhenotypes);
        geneA.setPhenotypeTerms( new HashSet<Term>() );
        geneA.addPhenotypeTerm(phenotypeData.getTermIncludingAlternatives("EP:04"));
        geneA.addPhenotypeTerm(phenotypeData.getTermIncludingAlternatives("EP:05"));
        
        // expect sum IC of most specific common terms for each gene phenotype
        double expResult = phenotypeData.getIC(phenotypeData.getTermIncludingAlternatives("EP:04"));
        expResult += phenotypeData.getIC(phenotypeData.getTermIncludingAlternatives("EP:05"));
        double result = phenotypeData.phenoMatchScore(terms, geneA);
        
        assertEquals(expResult, result, 0.001);
        
    }

    /**
     * Test of phenoGramScore method, of class PhenotypeData.
     */
    @Test
    public void testPhenoGramScore() {
        System.out.println("phenoGramScore");
        
        // build set of patient terms with only EP:06 from the example dataset
        HashSet<Term> patientTerms = new HashSet<Term>();
        patientTerms.add(phenotypeData.getTermIncludingAlternatives("EP:06"));
        
        // build gene A of example data set
        ArrayList<String> genePhenotypes = new ArrayList<String>();
        genePhenotypes.add("EP:04");
        genePhenotypes.add("EP:05");        
        Gene geneA = new Gene("chr1", 26, 32, "geneA", genePhenotypes);
        geneA.setPhenotypeTerms( new HashSet<Term>() );
        geneA.addPhenotypeTerm(phenotypeData.getTermIncludingAlternatives("EP:04"));
        geneA.addPhenotypeTerm(phenotypeData.getTermIncludingAlternatives("EP:05"));
        
        // build gene B of example data set
        ArrayList<String> geneBPhenotypes = new ArrayList<String>();
        geneBPhenotypes.add("EP:07");
        Gene geneB = new Gene("chr1", 26, 32, "geneB", geneBPhenotypes);
        geneB.setPhenotypeTerms( new HashSet<Term>() );
        geneB.addPhenotypeTerm(phenotypeData.getTermIncludingAlternatives("EP:07"));

        // add geneA and geneB to a genomic set of genes
        GenomicSet<Gene> genes = new GenomicSet<Gene>();
        genes.put("geneA", geneA);
        genes.put("geneB", geneB);
        
        // since geneA have higher phenomatchScore than geneB, we expect the
        // phenomatch score of geneB here again.
        double expResult = phenotypeData.getIC(phenotypeData.getTermIncludingAlternatives("EP:04"));
        expResult += phenotypeData.getIC(phenotypeData.getTermIncludingAlternatives("EP:05"));

        // calculate phenoMatch score for the patient temrs (EP:06) and the gene phenotypes
        double result = phenotypeData.phenoGramScore(patientTerms, genes);

        // assert eauality with a tolerance of 0.001
        assertEquals(expResult, result, 0.001);
        
    }

    /**
     * Test of mapTargetTermToGenes method, of class PhenotypeData.
     */
    @Test
    public void testMapTargetTermToGenes() {
        System.out.println("mapTargetTermToGenes");
         ArrayList<TargetTerm> targetTerms = exampleData.getTargetTermList();

        HashSet<String> genes = new HashSet();
        genes.add("geneA");
        genes.add("geneB");
        genes.add("geneD");
        
        HashMap<Term, HashSet<String>> expResult = new HashMap();
        expResult.put(phenotypeData.getTermIncludingAlternatives("EP:05"), genes);
        
        HashMap<Term, HashSet<String>> result = phenotypeData.mapTargetTermToGenes(targetTerms);
        assertEquals(expResult, result);
    }
    
    /**
     * test the getDirectRelation method from the Ontology class from the ontologizer package
     */
    @Test
    public void testTermRelation() {
    
        Term ep5 = phenotypeData.getTermIncludingAlternatives("EP:05");
        Term ep2 = phenotypeData.getTermIncludingAlternatives("EP:02");
        Term ep0 = phenotypeData.getTermIncludingAlternatives("EP:00");

        Ontology ep = phenotypeData.getOntology();
        
        // direct ancesotr
        assertTrue(TermRelation.IS_A.equals(ep.getDirectRelation(ep2.getID(), ep5.getID() )) );

        // indirect ancesotr
        assertFalse(TermRelation.IS_A.equals(ep.getDirectRelation(ep0.getID(), ep5.getID() )) );

    }

    /**
     * Test of isAncestorOrEqual method, of class PhenotypeData.
     */
    @Test
    public void testIsAncestorOrEqual() {
        System.out.println("isAncestorOrEqual");

        Term ep5 = phenotypeData.getTermIncludingAlternatives("EP:05");
        Term ep2 = phenotypeData.getTermIncludingAlternatives("EP:02");
        Term ep0 = phenotypeData.getTermIncludingAlternatives("EP:00");

        Ontology ep = phenotypeData.getOntology();
        assertTrue(phenotypeData.isAncestorOrEqual(ep5, ep2));
        assertTrue(phenotypeData.isAncestorOrEqual(ep5, ep0));
    }

    
}

