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

import jannovar.interval.Interval;
import java.util.Comparator;
import java.util.Objects;



/**
 * Implements a genomic element or genomic interval.
 * 
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class GenomicElement implements Comparable<GenomicElement>{
    
    // Genomic location in zero-based half-open BED-like format:
    private final String chr;
    private final int start;
    private final int end;
    private final String name;    
    
    /**
     * A Comparator that orders {@link GenomicElement} objects by there 
     * start coordinates.
     */
    public static final Comparator<GenomicElement> START_COORDINATE_ORDER = new StartCoordinateComparator();
    
    /**
     * Constructor for the {@link GenomicElement} class.
     * 
     * @param chr   Chromosome identifier
     * @param start the zero-based start coordinate
     * @param end   the end coordinate (zero-based half open, BED-like format)
     * @param name  a name for this genomic element/interval
     * @throws IllegalArgumentException     in case of negative coordinates or if end coordinate is smaller or equal start coordinate
     */
    public GenomicElement(String chr, int start, int end, String name)  throws IllegalArgumentException {
        
        this.chr = chr;
        this.name = name;
        
        // test for negative coordinates
        if (start < 0 | end < 0 ){

            throw new IllegalArgumentException(
                    "GenomicInterval constructor: start and end coordinates have to be positive" 
                    + start + " end=" + end 
                );
        }
        
        // test for proper start coordinates
        if ( start > end ){

            throw new IllegalArgumentException(
                    "GenomicInterval constructor: Start coordinate has to be "
                            + "smaller (or equal to) the end coordintae. start=" 
                            + start + " end=" + end 
                );
        }
        
        
        this.start = start;
        this.end = end;
        
    }

    
    /**
     * Convert element to {@link String} representation in the format {@code "name:chr:[start,end)"}.
     * 
     * @return String representation of the genomic interval object 
     */
    @Override
    public String toString(){
        return name + ":" + chr + ":[" + start + "," + end + ")";
    }
    
    /**
     * This function constructs {@link String} that represents an output line
     * for a TAB separated file like BED files. The line  contains the following
     * columns: chromosome, start, end, name. 
     * This function is overwritten by subclasses like {@link CNV} to output
     * additional columns for each element.
     * 
     * @return a TAB-separated output line to write BED like files.
     */
    public String toOutputLine(){
        return chr + "\t" + start + "\t" + end + "\t" + name;
    }
    
    /**
     * This functions returns a header line for a TAB-separated output file.
     * This function is overwritten by subclasses like {@link CNV}.
     * 
     * @return header line for tab separated output file 
     */
    public static String getOutputHeaderLine(){
        return "#chr\tstart\tend\tname";
    }
    
    
    /**
     * Convert {@link GenomicElement} into an {@link Interval} object.
     * Since the {@link Interval} and {@link IntervalTree} classes from the 
     * {@link jannovar} project is implemented for 1-based coordinates and 
     * assumes closed intervals (end position included), we subtracted 1 from 
     * the end coordinate.
     * 
     * @return {@link Interval} object for the {@link GenomicElement}
     * @throws java.lang.Exception
     */
    public Interval toInterval() throws Exception{
        
        // if interval has length zero, throw exception
        if (this.length() == 0){
            throw new Exception ("Interval of length zero are not supported.");
        }else{
            Interval iv = new Interval(start, end-1, this);
            return iv;
        }
    }
    
    /**
     * Test if another {@link GenomicElement} object is equal to this by
     * comparing the String representation of the objects.
     * 
     * @param other 
     * @return true if other {@link GenomicElement} object is equal to this 
     */
    @Override
    public boolean equals(Object other){
        
        if(this == other) return true;
        
        if(other == null || (this.getClass() != other.getClass())){
           return false;
        }
        
        GenomicElement elem = (GenomicElement) other;
        // compaire elements by there hash codes
        return this.hashCode() == elem.hashCode();
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 29 * hash + Objects.hashCode(this.chr);
        hash = 29 * hash + this.start;
        hash = 29 * hash + this.end;
        hash = 29 * hash + Objects.hashCode(this.name);
        return hash;
    }

    /**
     * Compares this {@link GenomicElement} with the specified {@link GenomicElement} for order. 
     * Thereby, it compares the String representation {@code "name:chr:[start,end)"} of the 
     * {@link GenomicElement} objects for order.
     * Returns a negative integer, zero, or a positive integer as this object 
     * is less than, equal to, or greater than the specified object.
     * 
     * @param o the {@link GenomicElement} to be compared. 
     * @return a negative integer, zero, or a positive integer as this object is less than, equal to, or greater than the specified object. 
     */
    @Override
    public int compareTo(GenomicElement o) {
        return this.toString().compareTo(o.toString());
    }
    

    /**
     * Test if the another {@link GenomicElement} has any overlap with this element.
     * 
     * @param other  An {@link GenomicElement} object that is tested for overlap
     * @return ture if any overlap else false
     */
    public boolean hasOverlap(GenomicElement other){
        
        // check if chromsomses for the two elements are equal
        boolean chrEquals = chr.equals(other.getChr());
        
        // check for any overlap: s1 < e2 and s2 < e1
        boolean overlap = (start < other.getEnd()) && (other.getStart() < end);
        
        return chrEquals &&  overlap ;
    }
    
    /**
     * Test if another {@link GenomicElement} overlaps this element completely.
     * 
     * @param other  
     *      An {@link GenomicElement} object that is tested for complete overlap
     * @return ture if input element overlaps this completely.
     */
    public boolean completeOverlaped(GenomicElement other){
        
        // check if chromsomses for the two elements are equal
        if ( ! chr.equals(other.getChr()) ){
            return false;
        }else{
        
            // check for complete overlap: s1 >= s2 and e1 <= e2
            return (this.start >= other.getStart()) && (this.end <= other.getEnd());

        }
    }

    /**
     * Test if another {@link GenomicElement} has reciprocal overlap >= a given 
     * fraction with this element.
     * 
     * @param other  
     *      An {@link GenomicElement} object that is tested for reciprocal overlap
     * @param fraction reciprocal ovelrap parameter (a double between 0 and 1)
     * @return ture if input element has reciprocal overlap >= {@code fraction} with this element.
     */
    public boolean reciprocalOverlap(GenomicElement other, double fraction){
        
        // check if chromsomses for the two elements are equal
        if ( ! chr.equals(other.getChr()) ){
            return false;
        }else{
        
            // get the length of the two elements in bp
            int n = this.length();
            int m = other.length();
            
            // get the number of bp in the overlp as min(end1, end2) - max(start1, start2)
            int overlap = Math.min(this.getEnd(), other.getEnd()) - Math.max(this.getStart(), other.getStart());
            
            // calculate the reciprocal overlap fraction:
            double reciprocalOverlap = 1.0 * overlap / Math.max(n, m);
            
            // test if overlap is greater or equal to the given input threshold
            return reciprocalOverlap >= fraction ;

        }
    }

    /**
     * Returns the length of this {@link GenomicElement} in base pairs (bp).
     * @return length in bp
     */
    public int length() {
        return this.end - this.start;
    }

    
    /**
     * @return the chr
     */
    public String getChr() {
        return chr;
    }

    /**
     * @return the start
     */
    public int getStart() {
        return start;
    }

    /**
     * @return the end
     */
    public int getEnd() {
        return end;
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * A comparison function, which imposes a total ordering on some collection 
     * of {@link GenomicElement}s by there start coordinates.
     * 
     * See {@see Comparator}.
     */
    private static class StartCoordinateComparator implements Comparator<GenomicElement> {

        @Override
        public int compare(GenomicElement e1, GenomicElement e2) {
            
            // get start coordinate
            int s1 = e1.getStart();
            int s2 = e2.getStart();
            
            // Returns a negative integer, zero, or a positive integer as the 
            // s1 is less than, equal to, or greater than s2.
            if (s1 < s2){
                return -1;
            }else if (s1 == s2) {
                return 0;
            }else{
                return 1;
            }
        }
    }
        
}
