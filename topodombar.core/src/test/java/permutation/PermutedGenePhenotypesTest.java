/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package permutation;

import java.io.IOException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import phenotypeontology.PhenotypeData;
import toyexampledata.ExampleData;

/**
 *
 * @author jonas
 */
public class PermutedGenePhenotypesTest {
    
    private static ExampleData exampleData;
    
    public PermutedGenePhenotypesTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() throws IOException {
        exampleData = new ExampleData();
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of permuteGenePhenotypes method, of class PermutedGenePhenotypes.
     */
    @Test
    public void testPermuteGenePhenotypes() {
        System.out.println("permuteGenePhenotypes");
        PhenotypeData orgPhenotypeData = exampleData.getPhenotypeData();

        PhenotypeData permPhenotypeData = PermutedGenePhenotypes.permuteGenePhenotypes(orgPhenotypeData);
        
        // check that not changed members are still the same
        assertEquals(orgPhenotypeData.getOntology(), permPhenotypeData.getOntology());
        assertEquals(orgPhenotypeData.getAllGenesIDs(), permPhenotypeData.getAllGenesIDs());
        
        for (String geneID : orgPhenotypeData.getAllGenesIDs()){
            
            System.out.println("TEST OUT: gene: " + geneID);
            System.out.println("TEST OUT: org phentoype: " + orgPhenotypeData.getGenePhenotypes(geneID));
            System.out.println("TEST OUT: perm phentoype: " + permPhenotypeData.getGenePhenotypes(geneID));
        }
    }
    
}
