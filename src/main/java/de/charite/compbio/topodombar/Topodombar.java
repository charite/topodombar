/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package de.charite.compbio.topodombar;

import genomicregions.AnnotateCNVs;
import genomicregions.Gene;
import genomicregions.GenomicElement;
import genomicregions.GenomicSet;
import io.TabFileParser;
import io.TabFileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.logging.Logger;
import ontologizer.go.Term;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
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
 * - Parse only domain data (as bed) and implement computation of boundaries from 
 * the domain regions
 * 
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class Topodombar {
    
    /**
     * The phenotype ontology instance (can be HPO or Uberpheno).
     */
    
    /** access to the phenotype ontology and gene to phenotype associations. */
    private static PhenotypeData phenotypeData;

    /** unspecific phenotype terms used to group patients. */
    private static HashSet<Term> targetTerms;
    
    /** mapping of some general targetTerms to the genes associated with 
     * corresponding targetTerm2targetGenes phenotypes.*/
    private static HashMap<Term, HashSet<String>> targetTerm2targetGenes;
    
    public static void main(String[] args) throws ParseException, IOException{
        
        // define commandline options
        // TODO this has to be moved to an separte class to enable to choose between GUI and CLI
        final CommandLineParser commandLineParser = new BasicParser();
	final Options options = new Options();
	Option help = new Option("h", "help", false, "print this (help-)message");
	options.addOption(help);
	options.addOption("i", "input-file", true, "input file with copy number variations (CNVs) in TAB separated file format");
	options.addOption("d", "domains", true, "topological domains in BED file format. Non-overlapping domain regions are assumed.");
	//options.addOption("b", "boundaries", true, "topological domain boundaries in BED file format.");
	options.addOption("g", "genes", true, "Genes in BED like format.");
	options.addOption("e", "enhancers", true, "Enhancers in BED like format.");
	options.addOption("O", "phenotype-ontology", true, "the phenotype ontology in obo file format.");
	options.addOption("a", "annotation-file", true, "phenotype annotation file that maps genes to phenotpye terms.");
	options.addOption("o", "output-file", true, "output file to which the annotated CNVs will be written.");
        	
        
        // TODO catch wrong arguments and print help massage

        // parse commandline arguments
        final CommandLine cmd = commandLineParser.parse(options, args);
        String cnvPath = cmd.getOptionValue("i");
        String domainPath = cmd.getOptionValue("d");
//        String boundaryPath = cmd.getOptionValue("b");
        String genesPath = cmd.getOptionValue("g");
        String enhancersPath = cmd.getOptionValue("e");
        String ontologyPath = cmd.getOptionValue("O");
        String annotationPath = cmd.getOptionValue("a");
        String outputPath = cmd.getOptionValue("o");
        
        // if help option is set, print usage
        HelpFormatter formatter = new HelpFormatter();
        if (cmd.hasOption(help.getOpt()) || cmd.hasOption(help.getLongOpt())) {
            formatter.printHelp(Topodombar.class.getSimpleName(), options);
            return;
	}
        
        // read the phenotype ontology
        phenotypeData = new PhenotypeData(ontologyPath, annotationPath);
        System.out.println("[INFO] Topodombar: Ontology and annotation table were parsed.");

        // read CNV data from input file:
        TabFileParser cnvParser = new TabFileParser(cnvPath);
        GenomicSet cnvs = cnvParser.parseCNVwithTerms(phenotypeData);
        targetTerms = cnvParser.parseTargetTermSet(phenotypeData);

        // bild mapping of targetTerms to targetGenes
        targetTerm2targetGenes = phenotypeData.mapTargetTermToGenes(targetTerms);
        
        System.out.println("DEBUG number of targetGenes:");
        for (Term tT:targetTerm2targetGenes.keySet()){
            System.out.println(tT.getIDAsString() + ": " + targetTerm2targetGenes.get(tT).size() + " target genes");
            System.out.println(tT.getIDAsString() + ": " + phenotypeData.getDescendants(tT).size() + " subterms" );
        }
        
        // read topological domain regions and compute boundaries from it.
        TabFileParser domainParser = new TabFileParser(domainPath);
        // parse topological domains 
        GenomicSet domains = domainParser.parse();

        // TODO: compute boundaries directely form input GenomicSet of domains
        // read topological domain regions and compute boundaries from it.
        GenomicSet boundaries = domainParser.parseBoundariesFromDomains();
        
        // annotate CNVs for overlap with boundary elements
        AnnotateCNVs.boundaryOverlap(cnvs, boundaries);
        
        // output message
        System.out.println("[INFO] Topodombar: Boundary overlap computed.");
        
        // read genes and compute overlap with genes
        GenomicSet<Gene> genes = new TabFileParser(genesPath).parseGeneWithTerms(phenotypeData);

        
        AnnotateCNVs.geneOverlap(cnvs, genes);
        System.out.println("[INFO] Topodombar: Gene overlap computed.");
        
        // compute phenogram score for overlaped genes:
        AnnotateCNVs.phenogramScore(cnvs, phenotypeData);
        
        // read enhancers 
        GenomicSet<GenomicElement> enhancers = new TabFileParser(enhancersPath).parse();
        
        // annotate CNVs definition of adjacent regions:
        AnnotateCNVs.defineAdjacentRegionsByDomains(cnvs, domains);
        
        // calculate phenogram score of adjacent genes
        AnnotateCNVs.phenogramScoreAdjacentGenes(cnvs, genes, phenotypeData);
        
        // define adjacent regions of the CNVs as the the region from
        // the CNV breakpoint to the end of the underling domain.
        AnnotateCNVs.defineAdjacentRegionsByDomains(cnvs, domains);

        // Calcuclate phenogram score for genes in these adjacent regions
        AnnotateCNVs.phenogramScoreAdjacentGenes(cnvs, genes, phenotypeData);

        // annotate CNVs as toplological domain boundary disruption (TDBD)
//        AnnotateCNVs.annotateTDBD(cnvs, enhancers);
        AnnotateCNVs.annotateTDBD(cnvs, enhancers, genes, targetTerm2targetGenes, phenotypeData);
        
        
        // write annotated CNVs to output file
        TabFileWriter outWriter = new TabFileWriter(outputPath);
        outWriter.write(cnvs);
        System.out.println("[INFO] Topodombar: Wrote annotated CNVs to output file '"+outputPath+"'.");
        
    }
    
}
