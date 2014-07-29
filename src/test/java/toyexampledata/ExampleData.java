/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package toyexampledata;

import genomicregions.CNV;
import genomicregions.Gene;
import genomicregions.GenomicElement;
import genomicregions.GenomicSet;
import io.TabFileParser;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import ontologizer.go.Term;
import phenotypeontology.PhenotypeData;

/**
 * This class provides toy example data that are read form the resources folder.
 * Here is an overview of the toy example data set:
 * 
 * 
 <pre>
 {@code 
    ############################################################################
                        10        20        30        40    
    Coord:    01234567890123456789012345678901234567890
    Domain    <---------->   <-------------------->        (0,12), (15,37)
    Boundary              ###                              (12,15)
    GeneC,B,D,A -->     --->    -->     ----->             C:(2,5), B:(10,14), D:(18,20), A:(26,32)
    Enhacer     **                      **                 (2,4), (26,28) 
    cnv1               ==========                          (9,19)
    cnv2              =========================            (8,33)
    cnv3                  =======                          (12,19)
    cnv4                    ==                             (14,16)
                        10        20        30        40    
              01234567890123456789012345678901234567890


    Example phenotype ontology (EPO) data for testing:

                        EP:00
                       /    \       
                    EP:01   EP:02        
                   /       /    \       
                EP:03   EP:04   EP:05
                           \    /  \
                            EP:06  EP:07

    Term     p   IC=-log(p)
    ========================
    EP:01    1   0
    EP:02    .25 1.39
    EP:03    .75 0.29
    EP:04    .25 1.39
    EP:05    .75 0.29
    EP:06    0   1.39 (Inf)
    EP:07    .25 1.39

    CNV phenotypes: EP:06
    CNV targetTerm: EP:05

    Gene    Annotation  PhenoMatchScore (to EP:06)
    ==================================
    GeneA   EP:04,EP:05 1.68
    GeneB   EP:07       0.29
    GeneC   EP:03       0
    GeneD   EP:05       0.29
    
    ############################################################################

* }
 </pre>
 * 
 *  * 
 * @author jonas
 */
public class ExampleData {
    /**
     * Example {@link PhenotypeData} instance from the toy example data set.
     */
    private final PhenotypeData phenotypeData;
    /**
     * Example CNV data set.
     */
    private final GenomicSet<CNV> cnvs;
    /**
     * Example gene data set.
     */
    private final GenomicSet<Gene> genes;
    /**
     * Example topological domain data set.
     */
    private final GenomicSet<GenomicElement> domains;
    /**
     * Example boundaries.
     */
    private final GenomicSet<GenomicElement> boundaries;
    
    /**
     * Example enhancer data
     */
    private final GenomicSet<GenomicElement> enhancers;
    
    /** target phenotype terms */
    private final HashSet<Term> targetTerms;
    
    /** Mapping of target terms to target genes */
    private final HashMap<Term, HashSet<String>> targetTerm2targetGene;
    

    /**
     * This constructor parses the example data and build up the example data
     * structures provided in this class.
     * @throws java.io.IOException
     */
    public ExampleData() throws IOException{
        
        // create phenotypeData object form example data:
        java.net.URL oboURL = ExampleData.class.getResource("/example_ontology.obo");
        String oboPath = oboURL.getPath();    
        String annotPath = ExampleData.class.getResource("/example_genes_to_penotype.txt").getPath();
        
        // parse ontology and create phenotypeData object
        phenotypeData = new PhenotypeData(oboPath, annotPath);        

        // create parser for example CNV dataset
        String cnvPath = ExampleData.class.getResource("/example_CNV.tab").getPath();
        TabFileParser cnvParser = new TabFileParser(cnvPath);
        cnvs = cnvParser.parseCNVwithTerms(phenotypeData);
        targetTerms = cnvParser.parseTargetTermSet(phenotypeData);
        
        // create parser for example gene dataset
        String genePath = ExampleData.class.getResource("/example_genes.tab").getPath();
        TabFileParser geneParser = new TabFileParser(genePath);
        genes = geneParser.parseGeneWithTerms(phenotypeData);

        targetTerm2targetGene = phenotypeData.mapTargetTermToGenes(targetTerms);

        // create parser for domain example dataset
        String domainPath = ExampleData.class.getResource("/example_domains.tab").getPath();
        domains = new TabFileParser(domainPath).parse();
        
        // create parser for boundary example dataset
        String boundaryPath = ExampleData.class.getResource("/example_boundary.tab").getPath();
        boundaries = new TabFileParser(boundaryPath).parse();

        // parse for enahncer example dataset
        String enhancerPath = ExampleData.class.getResource("/example_enhancer.tab").getPath();
        enhancers = new TabFileParser(enhancerPath).parse();

        
    }

    /**
     * Example {@link PhenotypeData} instance from the toy example data set.
     * @return the phenotypeData
     */
    public PhenotypeData getphenotypeData() {
        return phenotypeData;
    }

    /**
     * Example CNV data set.
     * @return the cnvs
     */
    public GenomicSet<CNV> getCnvs() {
        return cnvs;
    }

    /**
     * Example gene data set.
     * @return the genes
     */
    public GenomicSet<Gene> getGenes() {
        return genes;
    }

    /**
     * Example topological domain data set.
     * @return the domains
     */
    public GenomicSet<GenomicElement> getDomains() {
        return domains;
    }

    /**
     * Example boundaries.
     * @return the boundaries
     */
    public GenomicSet<GenomicElement> getBoundaries() {
        return boundaries;
    }

    /**
     * Example enhancer data
     * @return the enhancers
     */
    public GenomicSet<GenomicElement> getEnhancer() {
        return enhancers;
    }

    /**
     * target phenotype terms
     * @return the targetTerms
     */
    public HashSet<Term> getTargetTerms() {
        return targetTerms;
    }

    /**
     * Mapping of target terms to target genes
     * @return the targetTerm2targetGene
     */
    public HashMap<Term, HashSet<String>> getTargetTerm2targetGene() {
        return targetTerm2targetGene;
    }
    
}
