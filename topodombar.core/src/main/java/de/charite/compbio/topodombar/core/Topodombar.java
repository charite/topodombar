/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package de.charite.compbio.topodombar.core;

import genomicregions.AnnotateCNVs;
import static genomicregions.AnnotateGenes.addGeneSymbol;
import genomicregions.CNV;
import genomicregions.Gene;
import genomicregions.GenomicElement;
import genomicregions.GenomicSet;
import io.GeneSymbolParser;
import io.SimpleStatsWriter;
import io.TabFileParser;
import io.TabFileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import ontologizer.go.Term;
import phenotypeontology.PhenotypeData;
import phenotypeontology.TargetTerm;


/**
 * This is the main program class of the topodombar project.
 * 
 * TODOs:
 * - Parse genes with respect to different transcripts (use TranscriptModel from 
 * the Jannovar project).
 * 
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class Topodombar {
    
        
    /**
     * The phenotype ontology instance (can be HPO or Uberpheno).
     */
    
    /** access to the phenotype ontology and gene to phenotype associations. */
    private PhenotypeData phenotypeData;

    /** Input copy number variations (CNVs) to annotate with effect mechanism*/ 
    private GenomicSet<CNV> cnvs;

    /** Benign control copy number variations (CNVs) used to filter input CNVs*/ 
    private GenomicSet<CNV> controlCNVs;

    /** unspecific phenotype terms used to group patients. */
    private ArrayList<TargetTerm> targetTerms;

    /** mapping of {@link TargetTerms} to the subsets of CNVs annotated to the term. */
    private HashMap<TargetTerm, GenomicSet<CNV>> targetTerm2cnvSubset;
    
    /** mapping of some general targetTerms to the genes associated with 
     * corresponding targetTerm2targetGenes phenotypes.*/
    private HashMap<Term, HashSet<String>> targetTerm2targetGenes;
        
    /** Topological domains regions */
    private GenomicSet<GenomicElement> domains;

    /** Topological domain boundary regions*/
    private GenomicSet<GenomicElement> boundaries;
    
    /** Genes */
    private GenomicSet<Gene> genes;

//    /** Enhancers (or other regulatory elements) */
//    private GenomicSet<GenomicElement> enhancers;

    // some parameters:
    /** Size of adjacent regions, in case they are defined by size */
    private int regionSize;
    
    private String outputPath;
    
    public Topodombar(Map<String, Object> argMap) throws IOException{
        
        // get the individual values
        String ontologyPath = (String) argMap.get("phenotype_ontology");
        String annotationPath = (String) argMap.get("annotation_file");
        
        String cnvPath = (String) argMap.get("input_file");
        String domainPath = (String) argMap.get("domains");
        String genesPath = (String) argMap.get("genes");
        String enhancersPath = (String) argMap.get("enhancers");
        String targetTermPath = (String) argMap.get("target_terms");
        this. outputPath = (String) argMap.get("output_file");
        
        // parse optional arguments:
        this.regionSize = (Integer) argMap.get("adjacent_region_size");
        String globalPhenotype = (String) argMap.get("phenotype");
        String controlPath = (String) argMap.get("filter_out");
        String overlapFunc = (String) argMap.get("overlap_function");
        Double overlapFraction = (Double) argMap.get("overlap_fraction");

        // read the phenotype ontology
        this.phenotypeData = new PhenotypeData(ontologyPath, annotationPath);
        System.out.println("[INFO] Topodombar: Ontology and annotation table were parsed.");

        ////////////////////////////////////////////////////////////////////////
        //  CNVs
        ////////////////////////////////////////////////////////////////////////

        // read CNV data from input file:
        TabFileParser cnvParser = new TabFileParser(cnvPath);

        // if global phenotype for the entire cohort is given:
        if (globalPhenotype != null){

            // check if a singel enhancer file is given by -e
            Term globalTerm = this.phenotypeData.getTermIncludingAlternatives(globalPhenotype);
            if (globalTerm == null){
                System.err.printf("[ERROR] Could not find the input phenotype "
                        + "term with ID >%s< in the ontology that was read from "
                        + "this file >%s<. Exit now.%n",globalPhenotype, ontologyPath);
                System.exit(1);
            }

            System.out.printf("[LOGGING] Use the phenotype term '%s' (%s) for "
                    + "all input CNVs.%n"
                    , globalTerm.getNamespaceAsAbbrevString(), globalTerm.getIDAsString());

            // parse CNVs by using only one global term as phenotype annotation
            cnvs = cnvParser.parseCNVwithGlobalTerm(globalTerm);

        // assume individual phenotype annotation per patient in input CNV file
        }else{
            cnvs = cnvParser.parseCNVwithPhenotypeAnnotation(phenotypeData);
        }
        
        
        // check if only a single enhancer file is given.
        if (enhancersPath != null){

            GenomicSet<GenomicElement> enhancers = new TabFileParser(enhancersPath).parse();

            // build the list of target terms by using only the global term
            targetTerms = new ArrayList<TargetTerm>();
            
            // if global phenotype for the entire cohort is given:
            if (globalPhenotype != null){

                // build the list of target terms by using only the global term
                Term globalTerm = this.phenotypeData.getTermIncludingAlternatives(globalPhenotype);
                targetTerms.add(new TargetTerm(globalTerm, globalTerm.getNamespaceAsAbbrevString(), enhancers));

            // assume individual phenotype annotation per patient in input CNV file
            }else{

                
                // parse set of all target terms in the input CNV data set
                HashSet<Term> patientTargetTerms = cnvParser.parseTargetTermSet(phenotypeData);
                
                // for all target temrs defined in the patitens add a target Term to the targetTerm list
                // and use the single set of enhancers for all
                for (Term t : patientTargetTerms){
                    targetTerms.add(new TargetTerm(t, t.getNamespaceAsAbbrevString(), enhancers));
                }
            }
        
        // if ther are enhancers specific for each target term (by -t option):
        }else{
            
            // read target terms and corresponding enhancer data
            TabFileParser targetTermParser = new TabFileParser(targetTermPath);
            targetTerms = targetTermParser.parseTargetTerms(phenotypeData);            
            
            // if global phenotype for the entire cohort is given:
            if (globalPhenotype != null){
                
                // make set of all input term IDs
                HashSet<String> targetTermIDs = new HashSet<String>();
                for (TargetTerm t : targetTerms){
                    targetTermIDs.add(t.getIDAsString());
                }
                
                // check if single input term for patients is contained in the file with target terms
                if (! targetTermIDs.contains(globalPhenotype)){
                    System.err.printf("[ERROR] Could not find the input phenotype "
                            + "term with ID >%s< in the given list of target terms "
                            + "in file >%s<. Exit now.%n",globalPhenotype, targetTermPath);
                    System.exit(1);
                    
                }

            }
        }

        ////////////////////////////////////////////////////////////////////////        
        // filter out CNVs that overlap with the control set, if given
        ////////////////////////////////////////////////////////////////////////
        if (controlPath != null){
            // parse control CNVs:
            TabFileParser controlCnvParser = new TabFileParser(controlPath);
            controlCNVs = controlCnvParser.parseCNV();
            cnvs = cnvs.filterOutOverlap(controlCNVs, overlapFunc, overlapFraction);
            
            // TODO: report in log file (or any stats file) how much CNVs are filterd out!

        }
        
        // get the mapping of target terms to the subset of CNVs corresponding 
        // to these terms
        targetTerm2cnvSubset = getTargetTermSubset(targetTerms, cnvs);
        
        ////////////////////////////////////////////////////////////////////////
        // Build mapping of targetTerms to targetGenes
        ////////////////////////////////////////////////////////////////////////
        // TODO: maybe better to add the set of target genes directly as member to the TargetTerm class

        this.targetTerm2targetGenes = phenotypeData.mapTargetTermToGenes(targetTerms);
        
        // TODO: report number of target terms and target genes in any log or stats file
        System.out.println("[LOGGING] number of targetGenes:");

        ////////////////////////////////////////////////////////////////////////
        //  Domains and Boundaries
        ////////////////////////////////////////////////////////////////////////

        // read topological domain regions and compute boundaries from it.
        TabFileParser domainParser = new TabFileParser(domainPath);
        // parse topological domains 
        domains = domainParser.parse();

        // TODO: compute boundaries directely form input GenomicSet of domains
        // read topological domain regions and compute boundaries from it.
        boundaries = domainParser.parseBoundariesFromDomains();
                
        ////////////////////////////////////////////////////////////////////////
        //  Genes
        ////////////////////////////////////////////////////////////////////////
        // read genes and compute overlap with genes
        genes = new TabFileParser(genesPath).parseGeneWithTerms(phenotypeData);        
        
        // add GeneSymbol to genes
        HashMap<String, String> entrezToSymbol = new GeneSymbolParser(annotationPath).parseEntrezToSymbol();
        addGeneSymbol(genes, entrezToSymbol);
        
        ////////////////////////////////////////////////////////////////////////
        //  Enhancers
        ////////////////////////////////////////////////////////////////////////
        
        // read enhancers 
//        enhancers = new TabFileParser(enhancersPath).parse();
    }
    
    /**
     * Runs the entire analysis.
     */
    public void runAnalysis(){
        
        // iterate over all target term and corresponding CNV subsets:
        for (TargetTerm targetTerm : targetTerm2cnvSubset.keySet()){
            
            // get subset of CNVs annotated to the target term
            GenomicSet<CNV> cnvSubset =  targetTerm2cnvSubset.get(targetTerm);
            
            // get set of enhancers corresponding to the target term (tissue class)
            GenomicSet<GenomicElement> enhancers = targetTerm.getEnhancers();
            
            // TODO use a logger class with log file!
            // output some numbers:
            Term tT = targetTerm.getTerm();
            System.out.printf("[LOGGING] Run analysis for "
                    + "target term '%s' with %d CNVs.%n", targetTerm.getName(), cnvSubset.size());
            System.out.printf(
                "[LOGGING] Corresponding target term '%s'(%s) has %d sub-terms and %d target genes.%n"
                ,tT.getName(), tT.getIDAsString(),  
                phenotypeData.getDescendants(tT).size(), targetTerm2targetGenes.get(tT).size() 
            );
            
            ////////////////////////////////////////////////////////////////////////
            // Annotate CNVs with all other annotation sets and compute phenogram scores
            ////////////////////////////////////////////////////////////////////////
            // annotate CNVs for overlap with boundary elements 
            // annotate CNVs with genes that are completely overlapped by the CNV
            // annotate CNVs with genes that lie in the adjacent regions
            // annotate CNVs with adjacent enhancers
            // compute phenogram score for overlapped and adjacent genes:
            AnnotateCNVs.annoateCNVsForTDBD(cnvSubset, domains, boundaries, genes, enhancers, phenotypeData);

            // annotate Overlap region of each CNV with boundaries, genes, and enhancers and phenogram score:
            AnnotateCNVs.annoateOverlap(cnvSubset, boundaries, genes, enhancers, phenotypeData);

            ////////////////////////////////////////////////////////////////////////
            //  Toplological domain boundary disruption (TDBD)
            ////////////////////////////////////////////////////////////////////////

            // define adjacent regions as the interval form CNV breakpoint to the end
            // of the underlying toplological domain.
            AnnotateCNVs.defineAdjacentRegionsByDomains(cnvSubset, domains);

            // annotate the adjacent regions wiht enhancers, genes, and pheogram scores
            AnnotateCNVs.annoateAdjacentRegions(cnvSubset, boundaries, genes, enhancers, phenotypeData);

            // annotate CNVs as toplological domain boundary disruption (TDBD)
            AnnotateCNVs.annotateTDBD(cnvSubset, targetTerm2targetGenes);

            // annotate CNVs with new definition of TDBD (just by score, without the need of target phenotype)
            AnnotateCNVs.annotateTDBDjustByScore(cnvSubset);

            ////////////////////////////////////////////////////////////////////////
            // Enhancer adoption (EA)mechansim based on fixed size adjacent regions 
            // and without boundary and domain data
            ////////////////////////////////////////////////////////////////////////

            // define adjacent regions as the 400kb windows flanking the CNV.
            AnnotateCNVs.defineAdjacentRegionsByDistance(cnvSubset, regionSize);

            // overwrite the annotation of adjacent regions wiht enhancers, genes, and pheogram scores
            AnnotateCNVs.annoateAdjacentRegions(cnvSubset, boundaries, genes, enhancers, phenotypeData);

            // annotate CNVs as enhancer adoption mechanism (EA)
            AnnotateCNVs.annotateEA(cnvSubset, targetTerm2targetGenes);

            // annotate CNVs as enhancer adoption mechanism with low similarity of 
            // overlapped genes (EAlowG)
            AnnotateCNVs.annotateEAlowG(cnvSubset, targetTerm2targetGenes);


            ////////////////////////////////////////////////////////////////////////
            // Test duplications for Enhancer adoption effect by assuming tandem dups
            ////////////////////////////////////////////////////////////////////////
            // define the regions of interest
            AnnotateCNVs.defineOverlappedDomainRegions(cnvSubset, domains);
            // test effect mechanism
            AnnotateCNVs.tandemDuplicationEnhancerAdoption(cnvSubset, genes, enhancers, phenotypeData);
        }
    }
    
    /**
     * writes the results to the output file(s).
     * 
     * @throws IOException 
     */
    public void writeOutput() throws IOException{
        ////////////////////////////////////////////////////////////////////////
        // write annotated CNVs to output file
        ////////////////////////////////////////////////////////////////////////
        
        // TODO: Fix writing function to write the adjacent phenogram score fo the adjacent regions 
        // defined by the domain structure and not by the overwritten distance!!!
        
        TabFileWriter<CNV> outWriter = new TabFileWriter<CNV>(this.outputPath);
        outWriter.writeCNVs(this.cnvs, this.phenotypeData);
        System.out.println("[INFO] Topodombar: Wrote annotated CNVs to output file '"+this.outputPath+"'.");

        // calculate simple stats and write them to an additional output file:
        SimpleStatsWriter.calcAndWriteStats(this.outputPath + ".simple_stats.txt", cnvs);
        System.out.println("[INFO] Topodombar: Wrote simple count statistics to output file '"+this.outputPath + ".simple_stats.txt"+"'.");
        
    }
    
    /**
     * Helper function to get the subset of CNVs belonging to a given target term.
     * 
     * @param targetTerms an {@link ArrayList} of {@link TargetTerm}s
     * @param cnvs the entire {@link GenomicSet} of {@link CNV}s
     * @return a mapping from each target term to the subset of CNVs corresponding to this term.
     */
    private static HashMap<TargetTerm, GenomicSet<CNV>> getTargetTermSubset(ArrayList<TargetTerm> targetTerms, GenomicSet<CNV> cnvs){
   
        HashMap<TargetTerm, GenomicSet<CNV>> term2cnvs = new HashMap<TargetTerm, GenomicSet<CNV>>();
        
        // iterate over all terms and get subset:
        for(TargetTerm t: targetTerms){
            
            GenomicSet<CNV> subset = new GenomicSet<CNV>();
            // iterate over all cnvs and test if they belong to theis subset:
            for (Map.Entry<String, CNV> entry : cnvs.entrySet()){
                CNV cnv = entry.getValue();
                
                if (cnv.getTargetTerm().equals(t.getTerm())){
                    subset.put(entry.getKey(), cnv);
                }
            }
            
            // put subset to mapping
            term2cnvs.put(t, subset);
        }
        
        return term2cnvs;
    }
}
