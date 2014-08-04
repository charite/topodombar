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
//        
//        // define commandline ns
//        // TODO this has to be moved to an separte class to enable to choose between GUI and CLI
//        final CommandLineParser commandLineParser = new BasicParser();
//	final Options ns = new Options();
//	Option help = new Option("h", "help", false, "print this (help-)message");
//	ns.addOption(help);
//	ns.addOption("i", "input-file", true, "input file with copy number variations (CNVs) in TAB separated file format");
//	ns.addOption("d", "domains", true, "topological domains in BED file format. Non-overlapping domain regions are assumed.");
//	//options.addOption("b", "boundaries", true, "topological domain boundaries in BED file format.");
//	ns.addOption("g", "genes", true, "Genes in BED like format.");
//	ns.addOption("e", "enhancers", true, "Enhancers in BED like format.");
//	ns.addOption("O", "phenotype-ontology", true, "the phenotype ontology in obo file format.");
//	ns.addOption("a", "annotation-file", true, "phenotype annotation file that maps genes to phenotpye terms.");
//	ns.addOption("o", "output-file", true, "output file to which the annotated CNVs will be written.");
//        	
//        
//        // TODO catch wrong arguments and print help massage
//
//        // parse commandline arguments
//        final CommandLine cmd = commandLineParser.parse(ns, args);
//        String cnvPath = cmd.getOptionValue("i");
//        String domainPath = cmd.getOptionValue("d");
//        String genesPath = cmd.getOptionValue("g");
//        String enhancersPath = cmd.getOptionValue("e");
//        String ontologyPath = cmd.getOptionValue("O");
//        String annotationPath = cmd.getOptionValue("a");
//        String outputPath = cmd.getOptionValue("o");

//        // if help option is set, print usage
//        HelpFormatter formatter = new HelpFormatter();
//        if (cmd.hasOption(help.getOpt()) || cmd.hasOption(help.getLongOpt())) {
//            formatter.printHelp(Topodombar.class.getSimpleName(), ns);
//            return;
//	}


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
        argsParser.addArgument("-o", "--output-file").required(true)
                .help("output file to which the annotated CNVs will be written");
        argsParser.addArgument("-s", "--adjacent-region-size").type(Integer.class)
                .setDefault(400000).help("size in base pairs (bp) of adjacent regions");
        argsParser.addArgument("-v", "--version").action(Arguments.version());

//        regionSize = Integer.parseInt(cmd.getOptionValue("s"));
        // TODO: add parameter with default value for size of adjacent regions:
        //regionSize = 400000;
        
        Namespace ns = null;
        Map<String, Object> argMap = new HashMap<String, Object>(); 
        
        try{
            ns = argsParser.parseArgs(args);
        }catch (ArgumentParserException e) {
            argsParser.handleError(e);
            System.exit(1);
        }
        
        argMap = ns.getAttrs();
        System.out.println("DEBUG args ns: " + ns);
        System.out.println("DEBUG args argMap: " + argMap);

        String ontologyPath = ns.get("phenotype_ontology");
        String annotationPath = (String) argMap.get("annotation_file");
        System.out.println("DEBUG args ontologyPath: " + ontologyPath);
        System.out.println("DEBUG args annotPaht: " + annotationPath);
        
        String cnvPath = (String) argMap.get("input_file");
        String domainPath = (String) argMap.get("domains");
        String genesPath = (String) argMap.get("genes");
        String enhancersPath = (String) argMap.get("enhancers");
        String outputPath = (String) argMap.get("output_file");
        regionSize = (Integer) argMap.get("adjacent_region_size");
        
        // read the phenotype ontology
//        phenotypeData = new PhenotypeData(ontologyPath, annotationPath);ns.
        
        // read the phenotype ontology
//        phenotypeData = new PhenotypeData(ontologyPath, annotationPath);
        phenotypeData = new PhenotypeData(ontologyPath, annotationPath);
        System.out.println("[INFO] Topodombar: Ontology and annotation table were parsed.");

        // read CNV data from input file:
        TabFileParser cnvParser = new TabFileParser(cnvPath);
        cnvs = cnvParser.parseCNVwithTerms(phenotypeData);
        targetTerms = cnvParser.parseTargetTermSet(phenotypeData);

        // bild mapping of targetTerms to targetGenes
        targetTerm2targetGenes = phenotypeData.mapTargetTermToGenes(targetTerms);
        
        System.out.println("DEBUG number of targetGenes:");
        for (Term tT:targetTerm2targetGenes.keySet()){
            System.out.println(tT.getIDAsString() + ": " + targetTerm2targetGenes.get(tT).size() + " target genes");
            System.out.println(tT.getIDAsString() + ": " + phenotypeData.getDescendants(tT).size() + " subterms" );
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
        
    }
    
}
