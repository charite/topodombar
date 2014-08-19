/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package commandline;

import java.util.HashMap;
import java.util.Map;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.MutuallyExclusiveGroup;
import net.sourceforge.argparse4j.inf.Namespace;

/**
 *
 * @author jonas
 */
public class ArgumentParser {
    
    
    public static Map<String, Object> parseCommnadLineArguments(String [] args) throws ArgumentParserException{
        
        // build and argument parser and set its properties
        net.sourceforge.argparse4j.inf.ArgumentParser argsParser = ArgumentParsers.newArgumentParser("Topodombar")
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
        
        // mutually exclusive arguments (-e |-t) for a single enhancer file or a list of target terms:
        MutuallyExclusiveGroup enhancerArgs = argsParser.addMutuallyExclusiveGroup("enhancers").required(true);
        enhancerArgs.addArgument("-e", "--enhancers").help("Enhancers in BED like format");
        
        enhancerArgs.addArgument("-t", "--target-terms").help("tab separated file holding each target term (tissue) per line. "
                + "Columns should hold the corresponding target term as HPO term"
                + " ID, a unique tissue name, and the file paht to the corresponding "
                + "enhancer data (the file path may be absolut or relative to this file)");     

        argsParser.addArgument("-p", "--phenotype")
                .help("the term ID of a phenotpye that will be used to annotate "
                        + "the entire cohort of input CNVs (may only be used "
                        + "with -e option but not with -t)");

        
        argsParser.addArgument("-O", "--phenotype-ontology").required(true)
               .help("the phenotype ontology in OBO file format");
        argsParser.addArgument("-a", "--annotation-file").required(true)
                .help("phenotype annotation file that maps genes to phenotpye terms");

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

        // convert argument into an HashMap
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

        return argMap;

    }
    
}
