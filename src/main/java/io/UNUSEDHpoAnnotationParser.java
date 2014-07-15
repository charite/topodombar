/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import ontologizer.go.Ontology;
import ontologizer.go.Term;

/**
 *
 * @author jonas
 */
public class UNUSEDHpoAnnotationParser {
    
    // define some patterns
    private static final Pattern semicolon = Pattern.compile(";");
    private static final Pattern tabstopp = Pattern.compile("\t");

    private static final Pattern uphenoTermidPattern = Pattern.compile("\\(([HZM]P:\\d{7})\\)$");

    
    /**
     * Reads the HPO annotation file to map genes to phenotypic terms to which 
     * the gene is annotated.
     * 
     * @param annotationFilePath    annotation file from the HPO website
     * @param ontology  matching ontology as an {@link Ontology} object
     * @return a mapping from EntrezGene ID to a set of phenotpye terms to which
     * the gene is annotated to.
     */
    private static HashMap<String, HashSet<Term>> readAnnotations(String annotationFilePath, Ontology ontology) {

        HashMap<String, HashSet<Term>> gene2annot = new HashMap<String, HashSet<Term>>();
        try {

            BufferedReader in = new BufferedReader(new FileReader(annotationFilePath));

            String firstLine = in.readLine();
            boolean isHpo = false;
            boolean isUpheno = false;
            if (firstLine.startsWith("#Entrez Gene ID of human gene ; Gene symbol ; Annotated Uberpheno")) {
                    isUpheno = true;
            } else if (firstLine.startsWith("#Format: entrez-gene-id<tab>entrez-gene-symbol<tab>HPO-Term-Name<tab>HPO-Ter")) {
                    isHpo = true;
            } else {
                    throw new RuntimeException("Can't handle annotation-file format!");
            }

            String line = null;
            while ((line = in.readLine()) != null) {

                if (line.startsWith("#"))
                    continue;

                String[] split = null;
                String entrezId = null;
                String annotatedTermId = null;
                if (isUpheno) {
                    split = semicolon.split(line);
                    entrezId = split[0];

                    Matcher matcher = uphenoTermidPattern.matcher(split[2]);
                    if (!matcher.find()) {
                            System.out.println("pattern not matched in: " + line);
                    }
                    annotatedTermId = matcher.group(1);
                }

                if (isHpo) {
                    split = tabstopp.split(line);
                    entrezId = split[0];
                    annotatedTermId = split[3];
                }

                Term t = ontology.getTermIncludingAlternatives(annotatedTermId);
                if (t == null) {
                    System.err.println("Could not find term for ID:" + annotatedTermId + " parsed from line: " + line);
                    continue;
                }

                HashSet<Term> annotationsOfGene;
                if (gene2annot.containsKey(entrezId))
                    annotationsOfGene = gene2annot.get(entrezId);
                else
                    annotationsOfGene = new HashSet<Term>();
                    annotationsOfGene.add(t);
                    gene2annot.put(entrezId, annotationsOfGene);

            }
            in.close();

        } catch (IOException e) {
                e.printStackTrace();
        }

        return gene2annot;
    }

}

