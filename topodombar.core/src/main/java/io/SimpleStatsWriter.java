/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package io;

import genomicregions.CNV;
import genomicregions.GenomicSet;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 * Provides functionality to write a simple statistics about the number and percent of
 * explainable CNVs to an output file.
 * 
 * @author jonas
 */
public class SimpleStatsWriter {
    
    // define charset for nio.Files.write function
    private final static Charset charset = Charset.forName("utf-8");
    
    /**
     * Calculates simple count statistics about the number and percent of 
     * explainable CNVs and writes it to the output file.
     * 
     * @param outFile
     * @param cnvs 
     */
    public static void calcAndWriteStats(String outFile, GenomicSet<CNV> cnvs) throws IOException{
        
        ArrayList<String> outLines = new ArrayList<String>();
        outLines.add(StringUtils.join(new String[]{"#CNVs", "N", "percent", "size_mean(sd)"}, '\t'));
        
        // entire set of all cnvs
        Integer n = cnvs.size();
        String allMeanSD = meanSdSize(cnvs);
        outLines.add(StringUtils.join(new String[]{"All", n.toString(), "100.00", allMeanSD}, '\t'));
        
        // create subset of CNVs with boundary overlap:
        GenomicSet<CNV> withBoundaryOverlap = new GenomicSet<CNV>();
        for (CNV cnv : cnvs.values()){
            if (cnv.hasBoundaryOverlap()){
                withBoundaryOverlap.put(cnv.getName(), cnv);
            }
        }
        // write output line for CVNs with at leaste one boundary overlap
        Integer nBoundary = withBoundaryOverlap.size();
        String percentBoundary = Utils.roundToString( 100.0 * (double)nBoundary / (double)n);
        String sizeBoundary = meanSdSize(withBoundaryOverlap);
        outLines.add(StringUtils.join(new String[]{"Boundary", nBoundary.toString(), percentBoundary, sizeBoundary}, '\t'));
        
        for (String mechanismClass : new String [] {"TDBD", "newTDBD", "EA", "EAlowG"}){

            HashMap<String, Integer> counts = countMechanisms(cnvs, mechanismClass);
            
            // add intermediate header line for the new class
            outLines.add("#Analysis for effect mechanism class: " + mechanismClass );
            for (String mechanism : counts.keySet()){
                
                Integer nSubset = counts.get(mechanism);
                String percent = Utils.roundToString( 100.0 * (double)nSubset / (double)n );
                String meanSdStr = meanSdSizeOfSubsets(cnvs, mechanismClass, mechanism);
                
                // add vlaues to new output line
                outLines.add(StringUtils.join(new String[]{mechanism, nSubset.toString(), percent, meanSdStr}, '\t'));
                
            }
        }
        
        Path outPath =  Paths.get(outFile);
        java.nio.file.Files.write(outPath, outLines, charset);

    }
    
    private static HashMap<String, Integer> countMechanisms(GenomicSet<CNV> cnvs, String mechanismClass){
        
        HashMap<String, Integer> counts = new HashMap<String, Integer>();
        
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
    private static String meanSdSizeOfSubsets(GenomicSet<CNV> cnvs, String mechanismClass, String mechanism) {

        
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
    private static String meanSdSize(GenomicSet<CNV> cnvSet) {

        
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
