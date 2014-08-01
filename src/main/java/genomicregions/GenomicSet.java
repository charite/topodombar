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

import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import jannovar.interval.Interval;
import jannovar.interval.IntervalTree;
import org.apache.commons.lang3.StringUtils;



/**
 * Implements a set of {@link GenomicElement}s as an extension to the {@link HashMap} class.
 * The key is the name or ID of the element.
 * The value is a {@link GenomicElement} or any of its subclasses like {@link CNV}.
 * 
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 * @param <T> {@link GenomicElement} or any of its subclasses
 */
public class GenomicSet<T extends GenomicElement> extends HashMap<String, T>{
    
    /**
     * maps each chromosome to an intervalTree consisting of 
     * the GenomicElements of that chromosome
     */
    private HashMap<String, IntervalTree<T>> chr2tree;
    
    
    /**
     * Build an IntervalTree form the list of {@link GenomicElement}s for faster overlap search.
     * This method is called from the search methods if the interval tree is not already build.
     */
    private void buildIntervalTree (){

        this.chr2tree = new HashMap();
                
        // sort all elements by thier chromsoms
        HashMap<String, ArrayList<Interval<T>>> chr2intervalList = new HashMap();
        
        for (T e : this.values()){
            
            String chr = e.getChr();
            
            // test if chr is not already in chr2tree map:
            if (! chr2intervalList.containsKey(chr)){
                chr2intervalList.put(chr, new ArrayList<Interval<T>>());
            }
            
            // fill map:
            try {
                Interval<T> iv = e.toInterval();
                chr2intervalList.get(chr).add(iv);
            }
            catch (Exception excep){
                
                assert e.length() <= 0;
                
                System.err.println("WARNING: Failed to build Interval object form element with length: " + e.length());
                System.err.println("WARNING: Interval tree will not contain element with name: " + e.getName());
    
            }
        }
        
        // for each chromosome build an separate interval tree:
        for (String chr : chr2intervalList.keySet()){
            
            // fill map with IntervalTree that is build form the interval list
            this.chr2tree.put(chr, new IntervalTree(chr2intervalList.get(chr)));
            
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
        if (this.chr2tree == null){
            buildIntervalTree();
        }
        
        // Check if chromosm of input element is contained in theis set
        if (! this.chr2tree.containsKey(e.getChr())){
            
            // return empty list
            return result;
        
        }else{
            
            // search for overlapping intervals. Thereby convert from 0-based 
            // half open GenomicElement to 1-based closed Interval object
            List<T> resultList = this.chr2tree.get(e.getChr()).search(e.getStart(), e.getEnd()-1);
            
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
        for (T cand : anyOverlap(e).values()){
            
            // test for complete overlap with input element e:
            if (cand.completeOverlaped(e)){
                
                // add candidate to result set
                result.put(cand.getName(), cand);
            }
        }
        
        return result;
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
