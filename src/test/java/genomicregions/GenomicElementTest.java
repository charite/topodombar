/*
 * Copyright (C) 2014 Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package genomicregions;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import jannovar.interval.Interval;
import jannovar.interval.IntervalTree;
import java.util.ArrayList;

/**
 *
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class GenomicElementTest {
    
    public GenomicElementTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
        GenomicElement gi5_16 = new GenomicElement("chr1", 5, 16, "instance");
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testGenmicInterval() {
        
        GenomicElement gi = new GenomicElement("chr1", 1, 10, "aName");
        assertTrue("Start smaller or equal end?", gi.getStart() <= gi.getEnd());
                

    }

    @Test(expected=IllegalArgumentException.class)
    public void testGenmicIntervalNegativCoordinate() {
            GenomicElement wrongGI = new GenomicElement("chr1", -1, 1, "wrongGI");
    }
    
    @Test(expected=IllegalArgumentException.class)
    public void testGenmicIntervalInvertedInterval() {
            GenomicElement wrongGI = new GenomicElement("chr1", 10, 1, "wrongGI");
    }

    /**
     * Test of hasOverlap method, of class GenomicElement.
     */
    @Test
    public void testHasOverlap() {
        System.out.println("hasOverlap");
        GenomicElement instance = new GenomicElement("chr1", 5, 16, "instance");
        GenomicElement giPositive = new GenomicElement("chr1", 10, 20, "pos");
        GenomicElement giNegative = new GenomicElement("chr1", 20, 30, "neg");
        
        boolean posResult = instance.hasOverlap(giPositive);
        boolean negResult = instance.hasOverlap(giNegative);
        assertTrue("has overlap", posResult);
        assertFalse("no overlap", negResult);
    }

    /**
     * Test of getChr method, of class GenomicElement.
     */
    @Test
    public void testGetChr() {
        System.out.println("getChr");
        GenomicElement instance = new GenomicElement("chr1", 0, 1, "name");
        String expResult = "chr1";
        String result = instance.getChr();
        assertEquals(expResult, result);
    }

    /**
     * Test of getStart method, of class GenomicElement.
     */
    @Test
    public void testGetStart() {
        System.out.println("getStart");
        GenomicElement instance = new GenomicElement("chr1", 0, 1, "name");
        int expResult = 0;
        int result = instance.getStart();
        assertEquals(expResult, result);
    }

    /**
     * Test of getEnd method, of class GenomicElement.
     */
    @Test
    public void testGetEnd() {
        System.out.println("getEnd");
        GenomicElement instance = new GenomicElement("chr1", 0, 1, "name");
        int expResult = 1;
        int result = instance.getEnd();
        assertEquals(expResult, result);
    }

    /**
     * Test of getName method, of class GenomicElement.
     */
    @Test
    public void testGetName() {
        System.out.println("getName");
        GenomicElement instance = new GenomicElement("chr1", 0, 1, "name");
        String expResult = "name";
        String result = instance.getName();
        assertEquals(expResult, result);
    }

    /**
     * Test of toString method, of class GenomicElement.
     */
    @Test
    public void testToString() {
        System.out.println("toString");
        GenomicElement instance = new GenomicElement("chr1", 0, 1, "name");
        //String expResult = "name|chr1:0-0";
        String expResult = "name:chr1:[0,1)";
        String result = instance.toString();
        assertEquals(expResult, result);
    }

    /**
     * Test of equals method, of class GenomicElement.
     */
    @Test
    public void testEquals() {
        System.out.println("equals");
        GenomicElement instance = new GenomicElement("chr1", 0, 1, "name");
        GenomicElement other = new GenomicElement("chr1", 0, 1, "name");
        boolean expResult = true;
        boolean result = instance.equals(other);
        assertEquals(expResult, result);
    }
    
    /**
     * Test the import and functionallity of the jannovar.interval.Interval and IntervalTree classes
     * 
     */
    @Test
    public void testIntervalTree(){
        
        // a: ----
        // b:     ----
        // c:   ----
        /*
        GenomicElement a = new GenomicElement("chr1", 0, 100, "a");
        GenomicElement b = new GenomicElement("chr1", 100, 201, "b");
        GenomicElement c = new GenomicElement("chr1", 50, 150, "c");
        
        Interval<GenomicInterval> aIV = new Interval(a.getStart(), a.getEnd(), a);
        Interval<GenomicInterval> bIV = new Interval(b.getStart(), b.getEnd(), b);
        */
        Interval<String> a = new Interval(0, 100, "a");
        Interval<String> b = new Interval(100, 200, "b");

        ArrayList<Interval<String>> list = new ArrayList();
        list.add(a);
        
        IntervalTree<String> tree = new IntervalTree(list);
        ArrayList<String> foundList = tree.search(100, 200);
        
        System.out.println(foundList);
        
        boolean foundA = foundList.contains("a");
        assertTrue(foundA);
        
    }

    /**
     * Test of toOutputLine method, of class GenomicElement.
     */
    @Test
    public void testToOutputLine() {
        System.out.println("toOutputLine");
        GenomicElement instance = new GenomicElement("chr1", 0, 1, "name");
        String expResult = "chr1\t0\t1\tname";
        String result = instance.toOutputLine();
        assertEquals(expResult, result);
    }

    /**
     * Test of getOutputHeaderLine method, of class GenomicElement.
     */
    @Test
    public void testGetOutputHeaderLine() {
        System.out.println("getOutputHeaderLine");
        GenomicElement instance =  new GenomicElement("chr1", 0, 1, "name");
        String expResult = "#chr\tstart\tend\tname";
        String result = instance.getOutputHeaderLine();
        assertEquals(expResult, result);
    }

    /**
     * Test of toInterval method, of class GenomicElement.
     * @throws java.lang.Exception
     */
    @Test
    public void testToInterval() throws Exception {
        System.out.println("toInterval");
        GenomicElement instance =  new GenomicElement("chr1", 0, 10, "name");
        Interval expResult = new Interval(0, 9, instance);
        Interval result = instance.toInterval();
        assertEquals(expResult.toString(), result.toString());
    }

    /**
     * Test of toInterval method, of class GenomicElement with zero-length element.
     * @throws java.lang.Exception
     */
    @Test(expected = Exception.class)
    public void testToIntervalZeroLenght() throws Exception {
        System.out.println("toInterval with lengtz zero");
        GenomicElement zeroLength = new GenomicElement("chr1", 10, 10, "name");
        
        // toInterval() method should now throw an Excpeiont
        Interval zeroIV = zeroLength.toInterval();

    }

    /**
     * Test of completeOverlaped method, of class GenomicElement.
     */
    @Test
    public void testCompleteOverlaped() {
        System.out.println("completeOverlaped");
        GenomicElement other1 = new GenomicElement("chr1", 0, 100, "name");
        GenomicElement other2 = new GenomicElement("chr1", 15, 100, "name");
        GenomicElement instance = new GenomicElement("chr1", 10, 20, "name");
        boolean result1 = instance.completeOverlaped(other1);
        boolean result2 = instance.completeOverlaped(other2);
        assertEquals(true, result1);
        assertEquals(false, result2);
    }

    /**
     * Test of compareTo method, of class GenomicElement.
     */
    @Test
    public void testCompareTo() {
        System.out.println("compareTo");
        GenomicElement instance = new GenomicElement("chr1", 0, 100, "name");
        GenomicElement same = new GenomicElement("chr1", 0, 100, "name");
        GenomicElement other1 = new GenomicElement("chr1", 0, 200, "name");
        GenomicElement other2 = new GenomicElement("chr1", 0, 100, "name2");
        GenomicElement other3 = new GenomicElement("chr2", 0, 100, "name");
        GenomicElement other4 = new GenomicElement("chr0", 0, 100, "name");
        assertEquals(0, instance.compareTo(same));
        assertEquals(-1, instance.compareTo(other1));
        assertTrue(instance.compareTo(other2) > 0);
        assertTrue(instance.compareTo(other3) < 0 );
        assertTrue(instance.compareTo(other4) > 0);
    }

}
