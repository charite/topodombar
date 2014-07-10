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

import java.util.Collections;
import java.util.List;
import java.util.ArrayList;
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
public class GenomicSetTest {
    
    public GenomicSetTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
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
     * Test of anyOverlap method, of class GenomicSet.
     */
    @Test
    public void testAnyOverlap() {
        System.out.println("anyOverlap");
        GenomicElement query = new GenomicElement("chr1", 3, 6, "q");
                
         // fill set with two elements a and b
        GenomicSet<GenomicElement> instance = new GenomicSet<GenomicElement>();
        instance.put("a", new GenomicElement("chr1", 3, 5, "a"));
        instance.put("b", new GenomicElement("chr1", 5, 6, "b"));
        
       
        List<GenomicElement> expResult = new ArrayList<GenomicElement>(instance.values());
        Collections.sort(expResult);
        
        // add other elements that shoudl not be found
        instance.put("f", new GenomicElement("chr1", 1, 3, "f"));
        instance.put("g", new GenomicElement("chr1", 6, 8, "g"));
        instance.put("h", new GenomicElement("chr2", 3, 6, "h"));
        
        // serach for overlap
        List<GenomicElement> result = new ArrayList<GenomicElement>(instance.anyOverlap(query).values());
        Collections.sort(result);

        assertEquals(expResult, result);
    }

    /**
     * Test of completeOverlap method, of class GenomicSet.
     */
    @Test
    public void testCompleteOverlap() {
        System.out.println("completeOverlap");
        GenomicElement query = new GenomicElement("chr1", 3, 6, "q");
        
        // fill set with two elements a and b
        GenomicSet<GenomicElement> instance = new GenomicSet<GenomicElement>();
        instance.put("a", new GenomicElement("chr1", 3, 5, "a"));
        instance.put("b", new GenomicElement("chr1", 5, 6, "b"));
        
        List<GenomicElement> expResult = new ArrayList<GenomicElement>(instance.values());
        Collections.sort(expResult);
        
        // add other elements that shoudl not be found
        instance.put("c", new GenomicElement("chr1", 1, 4, "c"));
        instance.put("d", new GenomicElement("chr1", 5, 7, "d"));
        instance.put("e", new GenomicElement("chr1", 7, 10, "e"));
        instance.put("g", new GenomicElement("chr1", 1, 9, "h"));
        
        
        List<GenomicElement> result = new ArrayList<GenomicElement>(instance.completeOverlap(query).values());
        Collections.sort(result);
                
        assertEquals(expResult, result);
    
        GenomicSet<GenomicElement> emptySet = new GenomicSet<GenomicElement>();
        List<GenomicElement> emptyResult = new ArrayList<GenomicElement>( emptySet.completeOverlap(query).values());
        
        assertTrue("empty set yield empty result", emptyResult.isEmpty());
    
    }
    
}
