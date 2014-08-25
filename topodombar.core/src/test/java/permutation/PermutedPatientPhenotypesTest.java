/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package permutation;

import genomicregions.CNV;
import genomicregions.GenomicSet;
import java.io.IOException;
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
public class PermutedPatientPhenotypesTest {
    
    private static ExampleData exampleData;
    
    public PermutedPatientPhenotypesTest() {
    }
    
    @BeforeClass
    public static void setUpClass(){
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
     * Test of permutatePatientPhenotypes method, of class PermutedPatientPhenotypes.
     */
    @Test
    public void testPermutatePatientPhenotypes() {
        System.out.println("permutatePatientPhenotypes");
        GenomicSet<CNV> cnvs = exampleData.getCnvs();

        GenomicSet<CNV> permCNVs = PermutedPatientPhenotypes.permutatePatientPhenotypes(cnvs);
        
        assertEquals(cnvs.size(), permCNVs.size());
                
        assertEquals(cnvs.get("cnv1").getName(), permCNVs.get("cnv1").getName());
        assertEquals(cnvs.get("cnv1").getStart(), permCNVs.get("cnv1").getStart());
        assertEquals(cnvs.get("cnv1").getEnd(), permCNVs.get("cnv1").getEnd());
    }
    
}
