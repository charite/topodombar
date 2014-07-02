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



/**
 * Implements a set of {@link GenomicElement}s as an extension to the {@link HashMap} class.
 * The key is the name or ID of the element.
 * The value is a {@link GenomicElement} or any of its subclasses like {@link CNV}.
 * 
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 * @param <T> {@link GenomicElement} or any of its subclasses
 */
public class GenomicElementSet<T extends GenomicElement> extends HashMap<String, T>{
    
    // maps each crhomsome to and intervalTree consisting of the GenomicElements of that chromosome
    private HashMap<String, IntervalTree> chr2tree = null;
    
    
    /**
     * Build an IntervalTree form the list of {@link GenomicElement}s for faster overlap search.
     * This method is called from the search methods if the interval tree is not already build.
     */
    private void buildIntervalTree (){

        //TODO build separate Tree for each chromosome!
        chr2tree = new HashMap();
                
        // sort all elements by thier chromsoms
        HashMap<String, List<Interval>> chr2intervalList = new HashMap();
        
        for (T e : this.values()){
            
            String chr = e.getChr();
            
            // test if chr is not already in chr2tree map:
            if (! chr2intervalList.containsKey(chr)){
                chr2intervalList.put(chr, new ArrayList());
            }
            
            // fill map:
            chr2intervalList.get(chr).add(e.toInterval());
                        
        }
        
        // for each chromosome build an separate interval tree:
        for (String chr : chr2intervalList.keySet()){
            
            // fill map with IntervalTree that is build form the interval list
            chr2tree.put(chr, new IntervalTree(chr2intervalList.get(chr)));
            
        }
        
    }
    
    /**
     * Searches for all elements that have any overlap with input element {@code e}.
     * 
     * @param e GenomicElement for which the overlapping elements are searched for.
     * 
     * @return List of all overlapping GenomicElements
     */
    public List<GenomicElement> anyOverlap(GenomicElement e){

        List<GenomicElement> result = new ArrayList();
        
        // if interval tree is not build already, build it
        if (chr2tree == null){
            System.out.print("DEBUG null chr2tree:" + chr2tree);
            buildIntervalTree();
            System.out.print("DEBUG build chr2tree:" + chr2tree);
        }
        
        // Check if chromosm of input element is contained in theis set
        if (! chr2tree.containsKey(e.getChr())){
            
            // return empty list
            return result;
        
        }else{
            
            // search for overlapping intervals. Thereby convert from 0-based 
            // half open GenomicElement to 1-based closed Interval object
            return chr2tree.get(e.getChr()).search(e.getStart(), e.getEnd()-1);
            
        }
        
        
    }

    /**
     * Searches for all elements in the this set of elements that are completely overlapped by the input element {@code e}.
     * 
     * @param e GenomicElement for which the overlapping elements are searched for.
     * 
     * @return A List of {@link GenomicElement}s that are overlapped completely by {@code e}.
     */
    public List<GenomicElement> completeOverlap(GenomicElement e){
        //TODO implement same chromosom check!
        
        if (chr2tree == null){
            buildIntervalTree();
        }
        // convert from 0-based half open GenomicElement to 1-based closed Interval object
        List<GenomicElement> result =  anyOverlap(e);
        GenomicElement leftQuery = new GenomicElement(e.getChr(), e.getStart()-1, e.getStart(), "left");
        GenomicElement rightQuery = new GenomicElement(e.getChr(), e.getEnd(), e.getEnd()+1, "right");

        List<GenomicElement> halfOverlappingLeft = anyOverlap(leftQuery);
        List<GenomicElement> halfOverlappingRight = anyOverlap(rightQuery);
                    
        result.removeAll(halfOverlappingLeft);
        result.removeAll(halfOverlappingRight);
        
        return result;
    }
    
}
