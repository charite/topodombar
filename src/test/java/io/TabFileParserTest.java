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
import genomicregions.GenomicElement;
import genomicregions.GenomicSet;
import java.util.ArrayList;
import java.util.Arrays;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class TabFileParserTest {
    
    private static TabFileParser parser;
    
    public TabFileParserTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
	java.net.URL url = TabFileParserTest.class.getResource("/sample_CNV_chr22.tab");
        String path = url.getPath();
	parser = new TabFileParser(path);
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
        
        ArrayList<String> phenotypes = new ArrayList(Arrays.asList("HP:0001249", "HP:0000717", "HP:0001252"));
        
        GenomicElement firstCNV = new CNV("chr22", 49932021, 51187844, "132", "loss", phenotypes, "HP:0003011");
        GenomicSet<CNV> cnvs = parser.parseCNV();
        
        // parse the example CNV form the CNV set
        CNV example =  (CNV) cnvs.get("132");
        
        assertEquals(example.phenotpyes, phenotypes);
        assertEquals(example.targetTerm, "HP:0003011");
        assertTrue("contained element is parsed", firstCNV.equals(cnvs.get("132")));
        assertTrue("number of CNVs is 53, like the lines in the input file", cnvs.size() == 53);
    }
    
}
