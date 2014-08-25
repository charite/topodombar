/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package de.charite.compbio.topodombar.commandline;

import commandline.ArgumentParser;
import de.charite.compbio.topodombar.core.Topodombar;
import java.io.IOException;
import java.util.Map;
import net.sourceforge.argparse4j.inf.ArgumentParserException;

/**
 *
 * @author jonas
 */
public class TopodombarCMD {
    
        public static void main(String[] args) throws ArgumentParserException, IOException{
            
            // parse commandline arguments
            Map<String, Object> argMap = ArgumentParser.parseCommnadLineArguments(args);
            
            // run Topodombar tool with arguments:
            Topodombar topodombar = new Topodombar(argMap);
            
            topodombar.runAnalysis();
            
            topodombar.writeOutput();

            // run permutaion analysisi to get significance
            topodombar.runPermutations();
            
            
        }

}
