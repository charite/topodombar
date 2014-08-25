/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package io;

import genomicregions.CNV;
import genomicregions.GenomicSet;
import java.util.HashMap;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 *
 * @author jonas
 */
public class CountStatistics {
    
    /**
     * Counts for all CNV effect mechanism classes the number of annotated
     * mechanisms in the input {@link GenomicSet} of {@link CNV}s.
     * Note, the effect mechanism classes are defined in the CNV class.
     * @see CNV
     * @see SimpleStatsWriter
     * 
     * @param cnvs
     * @return 
     */
    public static HashMap<String, HashMap<String, Integer>> getEffectMechanismCounts(GenomicSet<CNV> cnvs){
        
        HashMap<String, HashMap<String, Integer>> counts = new HashMap<String, HashMap<String, Integer>>();
        
        // iterate over all effect mechanism classes defined in the {@link CNV}
        // class
        for (String mechanismClass : CNV.getEffectMechanismClasses()){
            
            // counts the annotated mechanism for the gieven class
            counts.put(mechanismClass, countMechanisms(cnvs, mechanismClass));
        }
        
        return counts;
    }

    /**
     * Counts the occurances of effect mechanisms for a given class in the input 
     * set of CNVs.
     * 
     * @param cnvs {@link GenomicSet} of {@CNV}s with annotated effect mechanism.
     * @param mechanismClass A {@link String} representing the effect mechanism class.
     * @return a {@link HashMap} that maps each found mechanism to its number of
     * occurances.
     */
    protected static HashMap<String, Integer> countMechanisms(GenomicSet<CNV> cnvs, String mechanismClass){
        
        HashMap<String, Integer> counts = new HashMap<String, Integer>();
        
        // initialize zero counts for all possible effect mechanisms:
        for (String effect : CNV.possibleEeffectAnnotations(mechanismClass)){
            counts.put(effect, 0);
        }
        
        for (CNV cnv : cnvs.values()){
            
            String mechanism = cnv.getEffectMechanism(mechanismClass);
            
            // count the occurrences
            int freq = counts.containsKey(mechanism) ? counts.get(mechanism) : 0;
            counts.put(mechanism, freq + 1);
        }
        
        return counts;
        
    }

    /** 
     * For each subset defined by effect mechanism, calculate the mean and sd of
     * the CNV sizes.
     * 
     * @param cnvs set of CNVs
     * @param mechanismClass effect mechanism class
     * @param keySet set of possible mechanism
     * @return a String in the format "mean(SD)" of the CNV sizes in the subsets
     */
    protected static String meanSdSizeOfSubsets(GenomicSet<CNV> cnvs, String mechanismClass, String mechanism) {

        
        DescriptiveStatistics sizes = new DescriptiveStatistics();

        for (CNV cnv: cnvs.values()){
            if (cnv.getEffectMechanism(mechanismClass).equals(mechanism)){
                sizes.addValue( cnv.length());                      
            }
        }
        // calculate mean and sd 
        Double mean = sizes.getMean();
        Double sd = sizes.getStandardDeviation();
        String meanSdStr = String.format("%s(+/-%s)", Utils.roundToString(mean), Utils.roundToString(sd));
            
        return meanSdStr;
    }

    /**
     * For the entire input set of CNVs, calculate the mean and sd of
     * the CNV sizes.
     * @param cnvSet
     * @return a String in the format "mean(SD)" of the CNV sizes
     */
    protected static String meanSdSize(GenomicSet<CNV> cnvSet) {

        
        DescriptiveStatistics sizes = new DescriptiveStatistics();

        for (CNV cnv: cnvSet.values()){
            sizes.addValue( cnv.length());                      
        }
        // calculate mean and sd 
        Double mean = sizes.getMean();
        Double sd = sizes.getStandardDeviation();
        String meanSdStr = String.format("%s(+/-%s)", Utils.roundToString(mean), Utils.roundToString(sd));
            
        return meanSdStr;
    }

}
