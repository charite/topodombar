/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package permutation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import ontologizer.go.Term;
import phenotypeontology.PhenotypeData;

/**
 * This class provides functionality to create background dataset with permuted
 * mapping of genes to phenotype for control analysis.
 * 
 * @author jonas
 */
public class PermutedGenePhenotypes {
    
    
    /**
     * returns a shallow copy of the input {@link PhenotypeData} object, but with 
     * permuted gene to phenotype mapping.
     * 
     * @param orgPhenotypeData
     * @return 
     */
    public static PhenotypeData permuteGenePhenotypes(PhenotypeData orgPhenotypeData){
        
        // get array of all genen IDs
        String [] orgGeneIDs = new String [orgPhenotypeData.getAllGenesIDs().size()];
        orgPhenotypeData.getAllGenesIDs().toArray(orgGeneIDs);

        // get list of correspoonding set of pheontype terms
        ArrayList<HashSet<Term>> orgGenePhenotypes = new ArrayList<HashSet<Term>>();
        for (String geneID : orgGeneIDs){
            orgGenePhenotypes.add(orgPhenotypeData.getGenePhenotypes(geneID));
        }
        
        // shuffle the list of phenotype terms
        Collections.shuffle(orgGenePhenotypes);
        
        // initialize permuted mapping of genes to phenotypes
        HashMap<String,HashSet<Term>> permutedGene2Phenotypes = new HashMap<String,HashSet<Term>>();
        
        // iterate over all genes wiht phenotype annotaiton:
        for (int i=0; i<orgGeneIDs.length; i++){
            
            String geneID = orgGeneIDs[i];
            
            // fill mapping with permuted phenotypes
            permutedGene2Phenotypes.put(geneID, orgGenePhenotypes.get(i));
        }
        
        // create a shallow copy of the original PhenotypeData object
        PhenotypeData permPhenotypeData = orgPhenotypeData.shallowCopy();
        
        // replace the original mapping of genes to phenotypes with the new permuted one
        permPhenotypeData.setGene2Terms(permutedGene2Phenotypes);
        
        return permPhenotypeData;

    }

}
