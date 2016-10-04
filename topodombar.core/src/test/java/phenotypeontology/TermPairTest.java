/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phenotypeontology;

import java.io.IOException;
import ontologizer.go.Term;
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
public class TermPairTest {
    
    private static TermPair termPair;

    private static ExampleData exampleData;
    
    private static Term p1;
    private static Term p2;
    private static Double s;
    

    public TermPairTest() throws IOException {
    }
    
    @BeforeClass
    public static void setUpClass() {

        // create TermPair object form example data:
        java.net.URL oboURL = PhenotypeDataTest.class.getResource("/example_ontology.obo");
        String oboPath = oboURL.getPath();
        
        String annotPath = PhenotypeDataTest.class.getResource("/example_genes_to_penotype.txt").getPath();
        
        // parse ontology and create PhenotypeData object
        PhenotypeData phenotypeData = new PhenotypeData(oboPath, annotPath);        
        
        p1 = phenotypeData.getOntology().getTerm("EP:04");
        p2 = phenotypeData.getOntology().getTerm("EP:05");
        s = 1.39;
        
        termPair = new TermPair(p1, p2, s);
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
     * Test of getPp method, of class TermPair.
     */
    @Test
    public void testGetPp() {
        System.out.println("getPp");
        TermPair instance = termPair;
        Term expResult = p1;
        Term result = instance.getPp();
        assertEquals(expResult, result);
    }

    /**
     * Test of getGp method, of class TermPair.
     */
    @Test
    public void testGetGp() {
        System.out.println("getGp");

        TermPair instance = this.termPair;
        Term expResult = p2;
        Term result = instance.getGp();
        assertEquals(expResult, result);

    }
    /**
     * Test of getS method, of class TermPair.
     */
    @Test
    public void testGetS() {
        System.out.println("getS");
        TermPair instance = termPair;
        Double expResult = s;
        Double result = instance.getS();
        assertEquals(expResult, result);
    }
    
}
