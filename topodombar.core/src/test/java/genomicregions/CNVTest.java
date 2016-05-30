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

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import ontologizer.go.Term;
import org.junit.After;
import org.junit.AfterClass;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import toyexampledata.ExampleData;

/**
 *
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class CNVTest {
    
    private static ExampleData exampleData;
    private static HashSet<Term> terms;
    private static Term targetTerm;
    private static CNV cnv1;
    private static CNV cnv2;
    private static CNV cnv3;
    private static CNV cnv4;
    
    public CNVTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() throws IOException {
        // construct CNV with all annotations:
        exampleData = new ExampleData();
        terms = exampleData.getPhenotypeData().getAllTerms();
        targetTerm = exampleData.getPhenotypeData().getTermIncludingAlternatives("EP:05");
        //cnvB = new CNV("chr1", 10, 101, "cnvB", "loss", terms, "HP:0003011");
        cnv1 = exampleData.getCnvs().get("cnv1");
        cnv2 = exampleData.getCnvs().get("cnv2");
        cnv3 = exampleData.getCnvs().get("cnv3");
        cnv4 = exampleData.getCnvs().get("cnv4");
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testSomeMethod() {
        
        // create CVN with the super class constructor:
        CNV cnvA = new CNV("chr1", 10, 101, "cnvA", "loss");
        CNV cnvB = new CNV("chr1", 10, 101, "cnvA", "loss", terms, cnv1.getTargetTerm());
        
        assertTrue("Default CNV phenotypes", cnvA.getPhenotypes().isEmpty());
        assertTrue("Default CNV targetTerm", cnvA.getTargetTerm() == null);
        

        assertEquals("CNV type", cnvA.getType(), "loss");
        assertEquals("phenotypes", cnvB.getPhenotypes(), terms);
        assertEquals("CNV targetTerm", cnvB.getTargetTerm(), cnv1.getTargetTerm());
    }

    /**
     * Test of toOutputLine method, of class CNV.
     * 
     */
    @Test
    public void testToOutputLine() {
        System.out.println("toOutputLine");
        
        CNV cnvB = new CNV("chr1", 10, 101, "cnvB", "loss", cnv1.getPhenotypes(), cnv1.getTargetTerm());
        String expResult = "chr1\t10\t101\tcnvB\tloss\tEP:0000006\tEP:0000005\t0\t.\t.\t.\t.\t.\t.\t.\t.\t.";
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
//    @Test
//    public void testGetOutputHeaderLine() {
//        System.out.println("getOutputHeaderLine");
//        String expResult = "#chr\tstart\tend\tname\ttype\tphenotypes\ttargetTerm\tboundaryOverlap\tgeneOverlap\toverlapPhenogramScore\tTDBD";
//        String result = cnv1.getOutputHeaderLine();
//        System.out.println(expResult);
//        System.out.println(result);
//        assertEquals(expResult, result);
//    }

    /**
     * Test of getType method, of class CNV.
     */
    @Test
    public void testGetType() {
        System.out.println("getType");
        String expResult = "loss";
        String result = cnv1.getType();
        assertEquals(expResult, result);
    }

    /**
     * Test of getPhenotypes method, of class CNV.
     */
    @Test
    public void testGetPhenotypes() throws IOException {
        System.out.println("getPhenotypes");
        HashSet<Term> expResult = new HashSet<Term>();
        expResult.add(exampleData.getPhenotypeData().getTermIncludingAlternatives("EP:06"));
        HashSet<Term> result = cnv1.getPhenotypes();
        assertEquals(expResult, result);
    }

    /**
     * Test of getTargetTerm method, of class CNV.
     */
    @Test
    public void testGetTargetTerm() {
        System.out.println("getTargetTerm");
        CNV instance = cnv1;
        Term expResult = targetTerm;
        Term result = instance.getTargetTerm();
        assertEquals(expResult, result);
    }

    /**
     * Test of getBoundaryOverlap method, of class CNV.
     */
    @Test
    public void testGetBoundaryOverlap() {
        System.out.println("getBoundaryOverlap");
        CNV instance = cnv1;
        GenomicSet<GenomicElement> result = instance.getBoundaryOverlap();
        assertEquals(new GenomicSet<GenomicElement>(), result);

    }

    
}
