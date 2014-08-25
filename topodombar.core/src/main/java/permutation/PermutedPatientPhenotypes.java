/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package permutation;

import genomicregions.CNV;
import genomicregions.GenomicSet;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import ontologizer.go.Term;

/**
 * This class provides functionality to permute the phenotype annotation of the 
 * patients as a background control experiment.
 * 
 * @author jonas
 */
public class PermutedPatientPhenotypes {
    
    
    /**
     * Returns a new {@link GenomicSet} of {@link CNV}s with randomly permuted
     * phenotype annotations. Thereby each CNV x in the input set get the 
     * phenoypes annotations (target term and set of individual term annotations)
     * of an randomly choosen other CNV y in the input set.
     * 
     * @param cnvs input set of CNVs with phenotpye annotations
     * @return an new {@link GenomicSet} of {@link CNV}s with randomly permuted
     * phenotype annotations.
     * 
     */
    public static GenomicSet<CNV> permutatePatientPhenotypes(GenomicSet<CNV> cnvs){


        // get all phenotypes annotations as list
        ArrayList<PairTuple<Term, HashSet<Term>>> allPhenotypes = new ArrayList<PairTuple<Term, HashSet<Term>>>();

        for (CNV cnv : cnvs.values()){
            allPhenotypes.add(new PairTuple(cnv.getTargetTerm(), cnv.getPhenotypes()));
        }

        // shuffle phenotype annotations randomly 
        Collections.shuffle(allPhenotypes);

        // create new genomic element set of CNVs with permuted phenotypes.
        GenomicSet<CNV> permutedCNVs = new GenomicSet<CNV>();

        // get all CNV names as array
        String [] cnvNames = new String [cnvs.size()];
        cnvs.keySet().toArray(cnvNames);

        // iterate over all input cnvs and copy them wiht permuted phenotype
        for (int i=0 ; i < cnvs.size(); i++){

            String orgName = cnvNames[i];
            CNV orgCNV = cnvs.get(orgName);

            // get random phenotypes from shuffeld list
            Term newTargetTerm = allPhenotypes.get(i).x;
            HashSet<Term> newPhenotypes = allPhenotypes.get(i).y;

            // create a new CNV object an assign the permuted phenotypes
            CNV newCNV = new CNV(
                    orgCNV.getChr(), orgCNV.getStart(), orgCNV.getEnd(), orgName, 
                    orgCNV.getType(), newPhenotypes, newTargetTerm
                 );

            // add new CNV to set of permuted CNVs
            permutedCNVs.put(orgName, newCNV);

        }

        return permutedCNVs;
    }
    
    
}
