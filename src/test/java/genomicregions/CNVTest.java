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

package genomicregions;

import java.util.Arrays;
import java.util.List;
import org.junit.After;
import org.junit.AfterClass;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 *
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class CNVTest {
    
    List<String> terms;
    CNV cnvB;
    
    public CNVTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
        // construct CNV with all annotations:
        terms = Arrays.asList("HP:0001249", "HP:0000717", "HP:0001252");
        cnvB = new CNV("chr1", 10, 101, "cnvB", "loss", terms, "HP:0003011");
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testSomeMethod() {
        
        // create CVN with the super class constructor:
        CNV cnvA = new CNV("chr1", 10, 101, "cnvA");
        
        assertTrue("Default CNV type", cnvA.getType() == ".");
        assertTrue("Default CNV phenotypes", cnvA.getPhenotpyes().isEmpty());
        assertTrue("Default CNV targetTerm", cnvA.getTargetTerm() == ".");
        

        assertEquals("CNV type", cnvB.getType(), "loss");
        assertEquals("phenotypes", cnvB.getPhenotpyes(), terms);
        assertEquals("CNV targetTerm", cnvB.getTargetTerm(), "HP:0003011");
    }

    /**
     * Test of toOutputLine method, of class CNV.
     * 
     */
    @Test
    public void testToOutputLine() {
        System.out.println("toOutputLine");
//        terms = Arrays.asList("HP:0001249", "HP:0000717", "HP:0001252");
//        cnvB = new CNV("chr1", 10, 101, "cnvB", "loss", terms, "HP:0003011");
        String expResult = "chr1\t10\t101\tcnvB\tloss\tHP:0001249;HP:0000717;HP:0001252\tHP:0003011\t.\t.\t.";
        String result = cnvB.toOutputLine();
        assertEquals(expResult, result);
    }

    /**
     * Test of getOutputHeaderLine method, of class CNV.
     *         return super.getOutputHeaderLine()
                + "\t" 
                + StringUtils.join(new String[]{
                    "type", 
                    "phenotypes", 
                    "targetTerm",
                    "boundaryOverlap",
                    "geneOverlap",
                    "overlapPhenogramScore"
                }, '\t');
     */
    @Test
    public void testGetOutputHeaderLine() {
        System.out.println("getOutputHeaderLine");
        CNV cnvB = new CNV("chr1", 10, 101, "cnvB", "loss", terms, "HP:0003011");
        String expResult = "#chr\tstart\tend\tname\ttype\tphenotypes\ttargetTerm\tboundaryOverlap\tgeneOverlap\toverlapPhenogramScore";
        String result = cnvB.getOutputHeaderLine();
        System.out.println(expResult);
        System.out.println(result);
        assertEquals(expResult, result);
    }
    
}
