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
import java.util.Arrays;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
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
        String allMeanSD = CountStatistics.meanSdSize(cnvs);
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
        String sizeBoundary = CountStatistics.meanSdSize(withBoundaryOverlap);
        outLines.add(StringUtils.join(new String[]{"Boundary", nBoundary.toString(), percentBoundary, sizeBoundary}, '\t'));
        
        for (String mechanismClass : CNV.getEffectMechanismClasses()){

            HashMap<String, Integer> counts = CountStatistics.countMechanisms(cnvs, mechanismClass);
            
            // add intermediate header line for the new class
            outLines.add("#Analysis for effect mechanism class: " + mechanismClass );
            for (String mechanism : counts.keySet()){
                
                Integer nSubset = counts.get(mechanism);
                String percent = Utils.roundToString( 100.0 * (double)nSubset / (double)n );
                String meanSdStr = CountStatistics.meanSdSizeOfSubsets(cnvs, mechanismClass, mechanism);
                
                // add vlaues to new output line
                outLines.add(StringUtils.join(new String[]{mechanism, nSubset.toString(), percent, meanSdStr}, '\t'));
                
            }
        }
        
        Path outPath =  Paths.get(outFile);
        java.nio.file.Files.write(outPath, outLines, charset);

    }
    
    /**
     * writes a tab separated file with statistics for all effect mechanisms
     * for the actual and permuted data. 
     * @param actualCounts the counts of effect mechanisms for the real data
     * @param permutedCounts counts for the permuted data
     * @param permutations number of permutations
     * @param outFile path to the output file
     */ 
    public static void writePermutationStatisticsOLD(HashMap<String, HashMap<String, Integer>> actualCounts, HashMap<String, HashMap<String, Integer []>> permutedCounts, Integer permutations, String outFile){
        
        /*
        HeaderLine:     class1_Effect1  Class1_Effect2, ... Class2_Effect1....
        actual data     count           count               count
        permuted_Mean   mean            mean                mean
        permuted_SD     sd              sd                  sd
        permuted_p(>actual)
        
        */
        
        // initialize list with output lines
        ArrayList<String> outLines = new ArrayList<String>();
        outLines.add("#DataType");
        outLines.add("real_data");
        outLines.add("permutation_mean");
        outLines.add("permutation_SD");
        outLines.add("permutation_p-value");
        
        // iterate over all effect classes and possible effect and add vlaues to the output lines:
        for (String effectClass: CNV.getEffectMechanismClasses()){

            for (String effect: CNV.possibleEeffectAnnotations(effectClass)){
                
                Integer actualCount = actualCounts.get(effectClass).get(effect);
                
                // get stats oject
                double [] doubleArray = Utils.toDoubleArray(permutedCounts.get(effectClass).get(effect));
                DescriptiveStatistics permValues = new DescriptiveStatistics(doubleArray);
                Double pVal = Utils.fractionGraterEqualX(permValues.getValues(), (double)actualCount );
                Double smallestP = 1.0/permValues.getN();
                
                // fill lines with addional column
                outLines.set(0, outLines.get(0) + "\t" +  effectClass + "_" + effect);
                outLines.set(1, outLines.get(1) + "\t" +  actualCount);                        
                outLines.set(2, outLines.get(2) + "\t" +  permValues.getMean());
                outLines.set(3, outLines.get(3) + "\t" +  permValues.getStandardDeviation());
                outLines.set(4, outLines.get(4) + "\t" +  (pVal > 0.0 ? pVal : "<" + smallestP) );
                
            }
            
        }
        
        // write lines to output file
        Path outPath =  Paths.get(outFile);
        
        try {
            java.nio.file.Files.write(outPath, outLines, charset);
        } catch (IOException ex) {
            Logger.getLogger(SimpleStatsWriter.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    /**
     * writes a tab separated file with statistics for all effect mechanisms
     * for the actual and permuted data. 
     * @param actualCounts the counts of effect mechanisms for the real data
     * @param permutedCounts counts for the permuted data
     * @param permutations number of permutations
     * @param outFile path to the output file
     */ 
    public static void writePermutationStatistics(HashMap<String, HashMap<String, Integer>> actualCounts, HashMap<String, HashMap<String, Integer []>> permutedCounts, Integer permutations, String outFile){
        
        /*
        HeaderLine:     #mechanismClass effectMechanism     actual  permuted_mean   permuted_sd pVal    
        actual data     count           count               count
        permuted_Mean   mean            mean                mean
        permuted_SD     sd              sd                  sd
        permuted_p(>actual)
        
        */
        
        // initialize list with output lines
        ArrayList<String> outLines = new ArrayList<String>();
        outLines.add(StringUtils.join(new String[]{
            "#mechanismClass", "effectMechanism", "real_data", 
            "permutation_mean", "permutation_SD", "permutation_p-value"}, '\t'));

        // iterate over all effect classes and possible effect and add vlaues to the output lines:
        for (String effectClass: CNV.getEffectMechanismClasses()){

            for (String effect: CNV.possibleEeffectAnnotations(effectClass)){
                
                Integer actualCount = actualCounts.get(effectClass).get(effect);
                
                // get stats oject
                double [] doubleArray = Utils.toDoubleArray(permutedCounts.get(effectClass).get(effect));
                DescriptiveStatistics permValues = new DescriptiveStatistics(doubleArray);
                Double pVal = Utils.fractionGraterEqualX(permValues.getValues(), (double)actualCount );
                Double smallestP = 1.0/permValues.getN();
                
                // create output line for this mechanism
                String line = StringUtils.join(new String[]{
                        effectClass, effect, actualCount.toString(), 
                        Double.toString(permValues.getMean()), Double.toString(permValues.getStandardDeviation()),
                        (pVal > 0.0 ? Double.toString(pVal) : "<" + smallestP)
                    }, '\t');
                // append line to output
                outLines.add(line);
                
            }
            
        }
        
        // write lines to output file
        Path outPath =  Paths.get(outFile);
        
        try {
            java.nio.file.Files.write(outPath, outLines, charset);
        } catch (IOException ex) {
            Logger.getLogger(SimpleStatsWriter.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
}
