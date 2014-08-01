/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package genomicregions;

import io.GeneSymbolParser;
import java.io.IOException;
import java.util.HashMap;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import toyexampledata.ExampleData;

/**
 *
 * @author jonas
 */
public class AnnotateGenesTest {
    
    private static ExampleData exampleData;
    
    public AnnotateGenesTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() throws IOException {
        this.exampleData = new ExampleData();
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of addGeneSymbol method, of class AnnotateGenes.
     */
    @Test
    public void testAddGeneSymbol() throws IOException {
        System.out.println("addGeneSymbol");
        GenomicSet<Gene> genes = exampleData.getGenes();
        String annotPath = AnnotateGenesTest.class.getResource("/example_genes_to_penotype.txt").getPath();
        HashMap<String, String> idToSymbol = new GeneSymbolParser(annotPath).parseEntrezToSymbol();
        AnnotateGenes.addGeneSymbol(genes, idToSymbol);
        
        assertEquals("GeneA", genes.get("geneA").getSymbol());
        assertEquals("GeneD", genes.get("geneD").getSymbol());
    }
    
}
