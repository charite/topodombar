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

package io;

import genomicregions.CNV;
import genomicregions.Gene;
import genomicregions.GenomicElement;
import genomicregions.GenomicSet;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import ontologizer.go.Term;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import phenotypeontology.PhenotypeData;
import phenotypeontology.PhenotypeDataTest;
import toyexampledata.ExampleData;

/**
 *
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class TabFileParserTest {
    
    private static TabFileParser parser;
    private static TabFileParser cnvParser;
    private static TabFileParser geneParser;
    private static TabFileParser domainParser;
    private static TabFileParser boundaryParser;
    private static PhenotypeData phenotypeData;
    
    public TabFileParserTest() {
    }
    
    @BeforeClass
    public static void setUpClass() throws IOException {
	java.net.URL url = TabFileParserTest.class.getResource("/sample_CNV_chr22.tab");
        String path = url.getPath();
	parser = new TabFileParser(path);
        
        // create parser for example CNV dataset
        String cnvPath = TabFileParserTest.class.getResource("/example_CNV.tab").getPath();
        cnvParser = new TabFileParser(cnvPath);

        // create parser for example gene dataset
        String genePath = TabFileParserTest.class.getResource("/example_genes.tab").getPath();
        geneParser = new TabFileParser(genePath);
        
        // create parser for domain example dataset
        String domainPath = TabFileParserTest.class.getResource("/example_domains.tab").getPath();
        domainParser = new TabFileParser(domainPath);

        // create parser for domain example dataset
        String boundaryPath = TabFileParserTest.class.getResource("/example_boundary.tab").getPath();
        boundaryParser = new TabFileParser(boundaryPath);

        // create an phenotypeData object from the example dataset:
        java.net.URL oboURL = PhenotypeDataTest.class.getResource("/example_ontology.obo");
        String oboPath = oboURL.getPath();
        String annotPath = PhenotypeDataTest.class.getResource("/example_genes_to_penotype.txt").getPath();
       
        // parse ontology and create phenotype data object
        phenotypeData = new PhenotypeData(oboPath, annotPath);        
        
    
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
     * Test of parse method, of class TabFileParser.
     */
    @Test
    public void testParse() throws Exception {
        System.out.println("parse");
        GenomicElement containedElem = new GenomicElement("chr22", 49932021, 51187844, "132");
        GenomicSet<GenomicElement> elements = parser.parse();
        
        assertTrue("contained element is parsed", containedElem.equals(elements.get("132")));
        assertTrue("number of CNVs is 53, like the lines in the input file", elements.size() == 53);
    
    }

    /**
     * Test of parseCNV method, of class TabFileParser.
     */
    @Test
    public void testParseCNV() throws Exception {
        System.out.println("parseCNV");
        
        // fist CNV in sample file:
        // chr22	49932021	51187844	132	loss	HP:0001249;HP:0000717;HP:0001252	HP:0003011
                
        CNV firstCNV = new CNV("chr22", 49932021, 51187844, "132", "loss");
        GenomicSet<CNV> cnvs = parser.parseCNV();
        
        // parse the example CNV form the CNV set
        CNV example =  cnvs.get("132");
        
        assertTrue("contained element is parsed", firstCNV.equals(cnvs.get("132")));
        assertTrue("number of CNVs is 53, like the lines in the input file", cnvs.size() == 53);
    }
    /**
     * Test of parseCNV method, of class TabFileParser.
     */
    @Test
    public void testParseCNVwithTerms() throws Exception {
        System.out.println("parseCNV");
        
        ExampleData exampleData = new ExampleData();
        HashSet<Term> terms = exampleData.getPhenotypeData().getAllTerms();
        Term targetTerm = exampleData.getPhenotypeData().getTermIncludingAlternatives("EP:05");                

        // create cnv1 from examle dataset from scratch
        Term t5 = phenotypeData.getTermIncludingAlternatives("EP:05");
        Term t6 = phenotypeData.getTermIncludingAlternatives("EP:06");
        
        CNV cnv1 = new CNV("chr1", 9, 19, "cnv1", "loss", new HashSet<Term>(Arrays.asList(t6)), t5);
                
        // parse CNVs from example data with phenotype terms
        GenomicSet<CNV> cnvs = cnvParser.parseCNVwithPhenotypeAnnotation(phenotypeData);
        
        // fetch cnv1
        CNV parsedCnv1 =  cnvs.get("cnv1");
        
        assertEquals(parsedCnv1.getPhenotypes(), cnv1.getPhenotypes());
        assertEquals(parsedCnv1.getTargetTerm(), t5);
        assertTrue("contained element is parsed", cnv1.equals(parsedCnv1));
        assertTrue("number of CNVs is 4, like the lines in the input file of example dataset", cnvs.size() == 4);        
    }

    /**
     * Test of parseGene method, of class TabFileParser.
     * @throws java.io.IOException
     */
    @Test
    public void testParseGene() throws IOException {
        System.out.println("parseGene");
        
        GenomicSet<Gene> genes = geneParser.parseGene();
        
        // build gene A of example data set
        ArrayList<String> genePhenotypes = new ArrayList<String>();
        genePhenotypes.add("EP:04");
        genePhenotypes.add("EP:05");        
        Gene geneA = new Gene("chr1", 26, 32, "geneA", genePhenotypes);
        geneA.setPhenotypeTerms( new HashSet<Term>() );
        geneA.addPhenotypeTerm(phenotypeData.getTermIncludingAlternatives("EP:04"));
        geneA.addPhenotypeTerm(phenotypeData.getTermIncludingAlternatives("EP:05"));

        assertTrue(geneA.equals( genes.get("geneA") ));
    }

    /**
     * Test of parseGeneWithTerms method, of class TabFileParser.
     * @throws java.io.IOException
     */
    @Test
    public void testParseGeneWithTerms() throws IOException {
        System.out.println("parseGeneWithTerms");
        GenomicSet<Gene> genes = geneParser.parseGeneWithTerms(phenotypeData);
        
        // build gene A of example data set
        ArrayList<String> genePhenotypes = new ArrayList<String>();
        genePhenotypes.add("EP:04");
        genePhenotypes.add("EP:05");        
        Gene geneA = new Gene("chr1", 26, 32, "geneA", genePhenotypes);
        geneA.setPhenotypeTerms( new HashSet<Term>() );
        geneA.addPhenotypeTerm(phenotypeData.getTermIncludingAlternatives("EP:04"));
        geneA.addPhenotypeTerm(phenotypeData.getTermIncludingAlternatives("EP:05"));
        
        assertTrue(geneA.equals( genes.get("geneA") ));
        assertEquals(geneA.getPhenotypeTerms(), genes.get("geneA").getPhenotypeTerms());
    }

    /**
     * Test of parseBoundariesFromDomains method, of class TabFileParser.
     * @throws java.io.IOException
     */
    @Test
    public void testParseBoundariesFromDomains() throws IOException{
        System.out.println("parseBoundariesFromDomains");
        GenomicSet<GenomicElement> expResult = boundaryParser.parse();
        GenomicSet<GenomicElement> result = domainParser.parseBoundariesFromDomains();
        assertEquals(expResult.get("b_1").toString(), result.get("b_1").toString());
    }

    /**
     * Test of parseTargetTermSet method, of class TabFileParser.
     */
    @Test
    public void testParseTargetTermSet() throws Exception {
        System.out.println("parseTargetTermSet");
        
        HashSet<Term> expResult = new HashSet<Term>();
        expResult.add(phenotypeData.getTermIncludingAlternatives("EP:05"));

        HashSet<Term> result = cnvParser.parseTargetTermSet(phenotypeData);

        assertEquals(expResult, result);
        
    }
    
}
