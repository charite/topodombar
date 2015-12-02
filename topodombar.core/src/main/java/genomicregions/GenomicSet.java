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

import com.google.common.collect.ImmutableList;
import java.util.ArrayList;
import java.util.HashMap;
import de.charite.compbio.jannovar.impl.intervals.IntervalArray;
import java.util.Arrays;
import org.apache.commons.lang3.StringUtils;



/**
 * Implements a set of {@link GenomicElement}s as an extension to the {@link HashMap} class.
 * The key is the name or ID of the element.
 * The value is a {@link GenomicElement} or any of its subclasses like {@link CNV}.
 * This class provides several functionalities to compute overlaps with other 
 * genomic features. Internally it holds an {@link IntervalTree} form the {@link jannovar} package
 * for fast overlap computation.
 * 
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 * @param <T> {@link GenomicElement} or any of its subclasses
 */
public class GenomicSet<T extends GenomicElement> extends HashMap<String, T>{
    
    /**
     * maps each chromosome to an intervalTree consisting of 
     * the GenomicElements of that chromosome
     */
    private HashMap<String, IntervalArray<T>> chr2IntervalArray;
    
    
    /**
     * Build an IntervalArray form the list of {@link GenomicElement}s for faster overlap search.
     * This method is called from the search methods if the interval array is not already build.
     */
    private void buildIntervalArray (){
        
        /**
         * hold an {@IntervalArray} for each chromosome
         */
        this.chr2IntervalArray = new HashMap();
                
        // create hash map for sorting all elements by thier chromsoms
        HashMap<String, ArrayList<T>> chr2ArrayList = new HashMap<>();

        // distribute elements to chrom specific lists
        for (T e : this.values() ) {
            String chr = e.getChr();
            
            // test if chr2IntervalArray has already a chr as key
            if (!chr2ArrayList.containsKey(chr)) {
                chr2ArrayList.put(chr, new ArrayList<T>());
            }
            
            chr2ArrayList.get(chr).add(e);
        }

        // construct an IntervalArray for each chromosome
        for (String chr : chr2ArrayList.keySet()) {
            
            IntervalArray<T> iTree = new IntervalArray<>(chr2ArrayList.get(chr),
					new GenomicElementEndExtractor<T>());
            this.chr2IntervalArray.put(chr, iTree);
        }
                
    }
    
    /**
     * Searches for all elements that have any overlap with input element {@code e}.
     * 
     * @param e GenomicElement for which the overlapping elements are searched for.
     * 
     * @return List of all overlapping GenomicElements
     */
    public GenomicSet<T> anyOverlap(GenomicElement e){

        GenomicSet<T> result = new GenomicSet();
        
        // if interval tree is not build already, build it
        if (this.chr2IntervalArray == null){
            buildIntervalArray();
        }
        // Check if chromosm of input element is contained in theis set
        if (! this.chr2IntervalArray.containsKey(e.getChr())){
            
            // return empty list
            return result;
        
        }else{            
            // search for overlapping intervals. 
            ImmutableList<T> resultList = this.chr2IntervalArray.get(e.getChr()).findOverlappingWithInterval(e.getStart(), e.getEnd()).getEntries();
            
            // convert the elements form List to GenomicSet
            for (T elem : resultList){
                result.put(elem.getName(), elem);
            }
            return result;
            
        }
    }

    /**
     * Searches for all elements in the this set of elements that are completely 
     * overlapped by the input element {@code e}.
     * 
     * @param e GenomicElement for which the overlapping elements are searched for.
     * 
     * @return A List of {@link GenomicElement}s that are overlapped completely by {@code e}.
     */
    public GenomicSet<T> completeOverlap(GenomicElement e){

        // initialize empty GenomicSet
        GenomicSet<T> result = new GenomicSet<T>();
        
        // first, get all elements that have ANY overlap as candidates for complete overlap
        for (T cand : this.anyOverlap(e).values()){
            
            // test for complete overlap with input element e:
            if (cand.completeOverlapped(e)){
                
                // add candidate to result set
                result.put(cand.getName(), cand);
            }
        }
        
        return result;
    }

    /**
     * Searches for all elements in the this set of elements that have 
     * reciprocal overlap greater than a fraction with the input element {@code e}.
     * 
     * @param e GenomicElement for which the overlapping elements are searched for.
     * @param fraction fraction of minimal reciprocal overlap (between 0 and 1).
     * 
     * @return A List of {@link GenomicElement}s that reciprocal overlap >= {@code fraction} wiht {@code e}.
     */
    public GenomicSet<T> reciprocalOverlap(GenomicElement e, double fraction){

        // initialize empty GenomicSet
        GenomicSet<T> result = new GenomicSet<T>();
        
        // first, get all elements that have ANY overlap as candidates for complete overlap
        for (T cand : anyOverlap(e).values()){
            
            // test for reciprocal overlap with input element e:
            if (cand.reciprocalOverlap(e, fraction)){
                
                // add candidate to result set
                result.put(cand.getName(), cand);
            }
        }
        
        return result;
    }

    /**
     * Filters this {@link GenomicSet} for only those elements that do not overlap
     * with one of the input elements.
     * Thereby a user specific overlap function is used:
     * <ul>
     * <li> "any" filters out elements that have any overlap with one of the input elements.
     * <li> "complete" filters out elements that overlap one of the input elements completely.
     * <li> "reciprocal" filter out elements that have >= {@code overlapFraction} reciprocal 
     * overlap with one input element.
     * 
     * @param others a set of elements
     * @param overlapFunc String representation of the overlap function, should be one of "any", "complete", or "recliprocal50"
     * @param overlapFraction minimal overlap fraction for reciprocal overlap. This will be ignored in case of "complete" or "any" overlap functions.
     * @return a subset of this set with only those elements that do not overlap with the input elements
     */
    public GenomicSet<T> filterOutOverlap(GenomicSet<T> others, String overlapFunc, Double overlapFraction){
        
        // check for proper overlapFunction argument:
        if (! Arrays.asList(new String [] {"any", "complete", "reciprocal"}).contains(overlapFunc)){
            System.err.printf("[ERROR] Wrong overlap fanction argument for filterOutOverlap function: '%s'%n", overlapFunc);
            System.exit(1);
        }
        
        GenomicSet<T> filteredSet = new GenomicSet<T>();
        
        for (T elem : this.values()){
            
            boolean hasNoOverlap = false;
            
            // decide which overlap function should be used:
            if ("any".equals(overlapFunc)){
                
                // check if the current cnv has ANY overlap with one of the input elements
                hasNoOverlap = others.anyOverlap(elem).isEmpty();
                
            }

            if ("complete".equals(overlapFunc)){
                
                // check if the current cnv COMPLETELY overlaps one of the input elements
                hasNoOverlap = others.completeOverlap(elem).isEmpty();
                
            }

            // decide which overlap function should be used:
            if ("reciprocal".equals(overlapFunc)){
                
                // check if the current cnv has RECIPROCAL overlap >= 50% with one of the input elements
                hasNoOverlap = others.reciprocalOverlap(elem, overlapFraction).isEmpty();
                
            }
            
            // if the element has no overlap with the input elemets, add it to the
            // filtered output set
            if (hasNoOverlap){
                filteredSet.put(elem.getName(), elem);
            }
            
        }
        
        return filteredSet;
    }
    
    /**
     * Removes all input elements from this {@link GenomicSet}.
     * 
     * @param others set of elements that will be removed 
     */
    private void removeAll(GenomicSet<T> others) {
        for (String k : others.keySet()){
            this.remove(k);
        }
    }
    
    /**
     * Joins all element names as String separate by semicolons ';'.
     * If this element is empty, it returns a dot ".".
     * 
     * @return all member names separated by ';' or "." if set is empty  
     */
    public String allNamesAsString(){
        
        // if the this set is empty return only a dot ".":
        if (this.isEmpty()){
            
            return ".";
        
        }else{
            
            return StringUtils.join(this.keySet(), ';');
        
        }
    }
    
    public ArrayList<String> getOutputLines(){
        
        // convert elements to list of output lines
        ArrayList<String> lines = new ArrayList<String>();
        // iterate over each element and add a line for it to the output lines
        for ( T e : this.values() ){
            
            // call the memberfunction toOutputLine to convert each element to 
            // one output line in the appropriate format. 
            lines.add(e.toOutputLine());
        
        }
                
        return lines;
    }
    
}
