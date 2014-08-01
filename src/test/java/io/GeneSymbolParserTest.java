/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package io;

import java.util.HashMap;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import phenotypeontology.PhenotypeDataTest;

/**
 *
 * @author jonas
 */
public class GeneSymbolParserTest {
    
    private static GeneSymbolParser parser;
    
    public GeneSymbolParserTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
        
        // create parser for example data set
        String annotPath = PhenotypeDataTest.class.getResource("/example_genes_to_penotype.txt").getPath();        
        parser = new GeneSymbolParser(annotPath);
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
     * Test of parseEntrezToSymbol method, of class GeneSymbolParser.
     */
    @Test
    public void testParseEntrezToSymbol() throws Exception {
        System.out.println("parseEntrezToSymbol");
        GeneSymbolParser instance = this.parser;
        
        HashMap<String, String> expResult = new HashMap<String, String>();
        expResult.put("geneA", "GeneA");
        expResult.put("geneB", "GeneB");
        expResult.put("geneC", "GeneC");
        expResult.put("geneD", "GeneD");

        HashMap<String, String> result = instance.parseEntrezToSymbol();
        assertEquals(expResult, result);

    }
    
}
