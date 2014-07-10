/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package de.charite.compbio.topodombar;

import genomicregions.AnnotateCNVs;
import genomicregions.Gene;
import genomicregions.GenomicSet;
import io.TabFileParser;
import io.TabFileWriter;
import java.io.IOException;
import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;


/**
 *
 * @author Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 */
public class Topodombar {
    
    public static void main(String[] args) throws ParseException, IOException{
        
        // define commandline options
        // TODO this has to be moved to an separte class to enable to choose between GUI and CLI
        final CommandLineParser commandLineParser = new BasicParser();
	final Options options = new Options();
	Option help = new Option("h", "help", false, "print this (help-)message");
	options.addOption(help);
	options.addOption("i", "input-file", true, "input file with copy number variations (CNVs) in TAB separated file format");
	options.addOption("b", "boundaries", true, "topological domain boundaries in BED file format.");
	options.addOption("g", "genes", true, "Genes in BED like format.");
	options.addOption("o", "output-file", true, "output file to which the annotated CNVs will be written.");
        	
        
        // TODO catch wrong arguments and print help massage
        // parse commandline arguments
        final CommandLine cmd = commandLineParser.parse(options, args);
        String cnvPath = cmd.getOptionValue("i");
        String boundaryPath = cmd.getOptionValue("b");
        String genesPath = cmd.getOptionValue("g");
        String outputPath = cmd.getOptionValue("o");
        
        HelpFormatter formatter = new HelpFormatter();
        if (cmd.hasOption(help.getOpt()) || cmd.hasOption(help.getLongOpt())) {
            formatter.printHelp(Topodombar.class.getSimpleName(), options);
            return;
	}
        
        // read data from input files:
        TabFileParser cnvParser = new TabFileParser(cnvPath);
        GenomicSet cnvs = cnvParser.parseCNV();
        
        // read boundary data
        TabFileParser boundaryParser = new TabFileParser(boundaryPath);
        GenomicSet boundaries = boundaryParser.parse();
        
        // annotate CNVs for overlap with boundary elements
        AnnotateCNVs.boundaryOverlap(cnvs, boundaries);
        
        // output message
        System.out.println("[INFO] Topodombar: Boundary overlap computed.");
        

        // read genes and compute overlap with genes
        GenomicSet<Gene> genes = new TabFileParser(genesPath).parseGene();
        AnnotateCNVs.geneOverlap(cnvs, genes);
        System.out.println("[INFO] Topodombar: Gene overlap computed.");
        
        // write annotated CNVs to output file
        TabFileWriter outWriter = new TabFileWriter(outputPath);
        outWriter.write(cnvs);
        System.out.println("[INFO] Topodombar: Wrote annotated CNVs to output file.");
    }
    
}
