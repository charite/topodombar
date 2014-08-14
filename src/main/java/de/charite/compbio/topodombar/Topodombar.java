/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package de.charite.compbio.topodombar;

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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import ontologizer.go.Term;
import org.apache.commons.cli.ParseException;
import phenotypeontology.PhenotypeData;


/**
 * This is the main program class of the topodombar project.
 * 
 * TODOs:
 * - Handle non-unique IDs in input regions (CNVs, Genes, Boundaries)
 * - Catch wrong cmd-line args by proper help message.
 * - Parse genes with respect to different transcripts (use TranscriptModel from 
 * the Jannovar project.
 * 
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class Topodombar {
    
    /**
     * The phenotype ontology instance (can be HPO or Uberpheno).
     */
    
    /** access to the phenotype ontology and gene to phenotype associations. */
    private static PhenotypeData phenotypeData;

    /** Input copy number variations (CNVs) to annotate with effect mechanism*/ 
    private static GenomicSet<CNV> cnvs;

    /** Benign control copy number variations (CNVs) used to filter input CNVs*/ 
    private static GenomicSet<CNV> controlCNVs;

    /** unspecific phenotype terms used to group patients. */
    private static HashSet<Term> targetTerms;

    /** mapping of some general targetTerms to the genes associated with 
     * corresponding targetTerm2targetGenes phenotypes.*/
    private static HashMap<Term, HashSet<String>> targetTerm2targetGenes;
        
    /** Topological domains regions */
    private static GenomicSet<GenomicElement> domains;

    /** Topological domain boundary regions*/
    private static GenomicSet<GenomicElement> boundaries;
    
    /** Genes */
    private static GenomicSet<Gene> genes;

    /** Enhancers (or other regulatory elements) */
    private static GenomicSet<GenomicElement> enhancers;

    // some parameters:
    /** Size of adjacent regions, in case they are defined by size */
    private static int regionSize;
  
    public static void main(String[] args) throws ParseException, IOException{
        
        // build and argument parser and set its properties
        ArgumentParser argsParser = ArgumentParsers.newArgumentParser("Topodombar")
                .description("Phenotypic analysis of microdeletions and topological "
                        + "chromosome domain boundaries. These scripts are meant"
                        + " to document the analysis performed in Ibn-Salem J "
                        + "et al., \"Deletions of Chromosomal Regulatory Boundaries "
                        + "are Associated with Congenital Disease\", Genome Biology, 2014.")
                .epilog("2014 by Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>")
                .defaultHelp(true)
                .version("${prog} 0.0.1");    
        //TODO remove hard-coded version. e.g. by this approach:http://stackoverflow.com/questions/2469922/generate-a-version-java-file-in-maven
        argsParser.addArgument("-i", "--input-file").required(true)
                .help("input file with copy number variations (CNVs) in TAB separated file format");
        argsParser.addArgument("-d", "--domains").required(true)
                .help("topological domains in BED file format. Non-overlapping domain regions are assumed");
        
        argsParser.addArgument("-g", "--genes").required(true).help("Genes in BED like format");
        argsParser.addArgument("-e", "--enhancers").required(true).help("Enhancers in BED like format");
        argsParser.addArgument("-O", "--phenotype-ontology").required(true)
               .help("the phenotype ontology in OBO file format");
        argsParser.addArgument("-a", "--annotation-file").required(true)
                .help("phenotype annotation file that maps genes to phenotpye terms");
        argsParser.addArgument("-p", "--phenotype")
                .help("the term ID of a phenotpye that will be used to annotate the entire cohort of input CNVs");
        argsParser.addArgument("-o", "--output-file").required(true)
                .help("output file to which the annotated CNVs will be written");
        // add optional parameters
        argsParser.addArgument("--adjacent-region-size").type(Integer.class)
                .setDefault(400000).help("size in base pairs (bp) of adjacent regions");
        argsParser.addArgument("--filter-out").help("benign control CNVs used to filter out input CNVs that overlap with this controls.");
        argsParser.addArgument("--overlap-function").choices(new String [] {"any", "complete", "reciprocal"} )
                .setDefault("reciprocal").help("overlap function used to filter out CNVs that overlap with the contorle data set.");
        argsParser.addArgument("--overlap-fraction").type(Double.class).choices(Arguments.range(0, 1))
                .setDefault(0.5).help("minimal fraction of reciprocal overlap. This will be ignored for 'complete' or 'any' overlap functions. ");
        argsParser.addArgument("-v", "--version").action(Arguments.version());

        // build objects to parse the commandline to
        Namespace ns = null;
        Map<String, Object> argMap = new HashMap<String, Object>(); 
        
        // parse arguments and handle errors
        try{
            ns = argsParser.parseArgs(args);
        }catch (ArgumentParserException e) {
            argsParser.handleError(e);
            System.exit(1);
        }
        
        argMap = ns.getAttrs();
        // get the individual values
        String ontologyPath = ns.get("phenotype_ontology");
        String annotationPath = (String) argMap.get("annotation_file");
        
        String cnvPath = (String) argMap.get("input_file");
        String domainPath = (String) argMap.get("domains");
        String genesPath = (String) argMap.get("genes");
        String enhancersPath = (String) argMap.get("enhancers");
        String outputPath = (String) argMap.get("output_file");
        
        // parse optional arguments:
        regionSize = (Integer) argMap.get("adjacent_region_size");
        String globalPhenotype = (String) argMap.get("phenotype");
        String controlPath = (String) argMap.get("filter_out");
        String overlapFunc = (String) argMap.get("overlap_function");
        Double overlapFraction = (Double) argMap.get("overlap_fraction");

        // read the phenotype ontology
        phenotypeData = new PhenotypeData(ontologyPath, annotationPath);
        System.out.println("[INFO] Topodombar: Ontology and annotation table were parsed.");

        ////////////////////////////////////////////////////////////////////////
        //  CNVs
        ////////////////////////////////////////////////////////////////////////

        // read CNV data from input file:
        TabFileParser cnvParser = new TabFileParser(cnvPath);
        
        // if global phenotype for the entire cohort is given:
        if (globalPhenotype != null){
            Term globalTerm = phenotypeData.getTermIncludingAlternatives(globalPhenotype);
            if (globalTerm == null){
                System.err.printf("[ERROR] Could not find the input phenotype "
                        + "term with ID >%s< in the ontology that was read from "
                        + "this file >%s<. Exit now.%n",globalPhenotype, ontologyPath);
                System.exit(1);
            }
            cnvs = cnvParser.parseCNVwithGlobalTerm(globalTerm);
            targetTerms = new HashSet<Term>();
            targetTerms.add(globalTerm);
            
        // assume individual phenotype annotation per patient in input CNV file
        }else{
            cnvs = cnvParser.parseCNVwithPhenotypeAnnotation(phenotypeData);
            // parse set of all target terms
            targetTerms = cnvParser.parseTargetTermSet(phenotypeData);
        }
        
        // filter out CNVs that overlap with the control set, if given
        if (controlPath != null){
            // parse control CNVs:
            TabFileParser controlCnvParser = new TabFileParser(controlPath);
            controlCNVs = controlCnvParser.parseCNV();
            cnvs = cnvs.filterOutOverlap(controlCNVs, overlapFunc, overlapFraction);
            
            // TODO: report in log file (or any stats file) how much CNVs are filterd out!

        }
        // bild mapping of targetTerms to targetGenes
        targetTerm2targetGenes = phenotypeData.mapTargetTermToGenes(targetTerms);
        
        // TODO: report number of target terms and target genes in any log or stats file
        System.out.println("[LOGGING] number of targetGenes:");
        for (Term tT:targetTerm2targetGenes.keySet()){
            System.out.printf(
                "[LOGGING] Target term '%s'(%s) has %d sub-terms and %d target genes.%n"
                ,tT.getName(), tT.getIDAsString(),  
                phenotypeData.getDescendants(tT).size(), targetTerm2targetGenes.get(tT).size() 
            );
        }

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
        enhancers = new TabFileParser(enhancersPath).parse();
        
        ////////////////////////////////////////////////////////////////////////
        // Annotate CNVs with all other annotation sets and compute phenogram scores
        ////////////////////////////////////////////////////////////////////////
        // annotate CNVs for overlap with boundary elements
        // annotate CNVs with genes that are completely overlapped by the CNV
        // annotate CNVs with genes that lie in the adjacent regions
        // annotate CNVs with adjacent enhancers
        // compute phenogram score for overlaped and adjacent genes:
        AnnotateCNVs.annoateCNVsForTDBD(cnvs, domains, boundaries, genes, enhancers, phenotypeData);
               
        // annotate Overlap region of each CNV with boundaries, genes, and enhancers and phenogram score:
        AnnotateCNVs.annoateOverlap(cnvs, boundaries, genes, enhancers, phenotypeData);

        ////////////////////////////////////////////////////////////////////////
        //  Toplological domain boundary disruption (TDBD)
        ////////////////////////////////////////////////////////////////////////
        
        // define adjacent regions as the interval form CNV breakpoint to the end
        // of the underlying toplological domain.
        AnnotateCNVs.defineAdjacentRegionsByDomains(cnvs, domains);
        
        // annotate the adjacent regions wiht enhancers, genes, and pheogram scores
        AnnotateCNVs.annoateAdjacentRegions(cnvs, boundaries, genes, enhancers, phenotypeData);
        
        // annotate CNVs as toplological domain boundary disruption (TDBD)
        AnnotateCNVs.annotateTDBD(cnvs, targetTerm2targetGenes);

        // annotate CNVs with new definition of TDBD (just by score, without the need of target phenotype)
        AnnotateCNVs.annotateTDBDjustByScore(cnvs);
        
        ////////////////////////////////////////////////////////////////////////
        // Enhancer adoption (EA)mechansim based on fixed size adjacent regions 
        // and without boundary and domain data
        ////////////////////////////////////////////////////////////////////////
        
        // define adjacent regions as the 400kb windows flanking the CNV.
        AnnotateCNVs.defineAdjacentRegionsByDistance(cnvs, regionSize);

        // overwrite the annotation of adjacent regions wiht enhancers, genes, and pheogram scores
        AnnotateCNVs.annoateAdjacentRegions(cnvs, boundaries, genes, enhancers, phenotypeData);
            
        // annotate CNVs as enhancer adoption mechanism (EA)
        AnnotateCNVs.annotateEA(cnvs, targetTerm2targetGenes);

        // annotate CNVs as enhancer adoption mechanism with low similarity of 
        // overlaped genes (EAlowG)
        AnnotateCNVs.annotateEAlowG(cnvs, targetTerm2targetGenes);

        ////////////////////////////////////////////////////////////////////////
        // write annotated CNVs to output file
        ////////////////////////////////////////////////////////////////////////
        TabFileWriter<CNV> outWriter = new TabFileWriter<CNV>(outputPath);
        outWriter.writeCNVs(cnvs, phenotypeData);
        System.out.println("[INFO] Topodombar: Wrote annotated CNVs to output file '"+outputPath+"'.");

        // calculate simple stats and write them to an additional output file:
        SimpleStatsWriter.calcAndWriteStats(outputPath + ".simple_stats.txt", cnvs);
        System.out.println("[INFO] Topodombar: Wrote simple count statistics to output file '"+outputPath + ".simple_stats.txt"+"'.");
        
    }
    
}
