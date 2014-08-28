/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package annotation;

import genomicregions.CNV;
import genomicregions.Gene;
import genomicregions.GenomicElement;
import genomicregions.GenomicSet;
import io.TabFileParser;
import io.TabFileParserTest;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import ontologizer.go.Term;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import static org.junit.Assert.*;
import phenotypeontology.PhenotypeData;
import toyexampledata.ExampleData;

/**
 *
 * @author jonas
 */
public class InterpretCNVsTest {
    
    // declare some variables to be knwon in this test
    private static GenomicSet<CNV> cnvs; 
    private static GenomicSet<CNV> exampleCNVs; 
    private static GenomicSet<GenomicElement> domains;
    private static GenomicSet<GenomicElement> boundaries;
    //private static GenomicElementSet<GenGenomicSet  public AnnotateCNVsTest() {
    private static GenomicSet<Gene> genes;
    
    private static ExampleData exampleData;

    public InterpretCNVsTest() {
    }
    
    @BeforeClass
    public static void setUpClass() throws IOException {
        // read sample CNVs:
        String cnvPath = TabFileParserTest.class.getResource("/sample_CNV_chr22.tab").getPath();
        TabFileParser cnvParser = new TabFileParser(cnvPath);
        cnvs = cnvParser.parseCNV();

        // read sample boundary elements
        String boundaryPath = TabFileParserTest.class.getResource("/sample_boundary_chr22.tab.addOne").getPath();
        TabFileParser boundaryParser = new TabFileParser(boundaryPath);
        boundaries = boundaryParser.parse();
                
        // read sample genes elements
        String genePath = TabFileParserTest.class.getResource("/sample_genes_chr22.tab").getPath();
        TabFileParser geneParser = new TabFileParser(genePath);
        genes = geneParser.parseGene();

        // parse toy example data set
        exampleData = new ExampleData();
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of annotateEA method, of class AnnotateCNVs.
     */
    public void testAnnotateEA() {
        System.out.println("annotateEA");
        GenomicSet<CNV> cnvs = this.exampleData.getCnvs();
        GenomicSet<Gene> genes = this.exampleData.getGenes();
        GenomicSet<GenomicElement> enhancer = this.exampleData.getEnhancer();
        GenomicSet<GenomicElement> boundaries = this.exampleData.getBoundaries();
        GenomicSet<GenomicElement> domains = this.exampleData.getDomains();
        PhenotypeData phenotypeData = this.exampleData.getPhenotypeData();
        HashMap<Term, HashSet<String>> term2genes = this.exampleData.getTargetTerm2targetGene();
        int regionSize = 20;
        AnnotateCNVs.boundaryOverlap(cnvs, boundaries);
        AnnotateCNVs.defineAdjacentRegionsByDistance(cnvs, regionSize);
        AnnotateCNVs.annoateOverlap(cnvs, boundaries, genes, enhancer, phenotypeData);
        AnnotateCNVs.annoateAdjacentRegions(cnvs, genes, enhancer, phenotypeData);
        InterpretCNVs.annotateEA(cnvs, term2genes);
        assertEquals("EA", cnvs.get("cnv1").getEffectMechanism("EA"));
    }

    /**
     * Test of annotateEAlowG method, of class AnnotateCNVs.
     */
    public void testAnnotateEAlowG() {
        System.out.println("annotateEAlowG");
        GenomicSet<CNV> cnvs = this.exampleData.getCnvs();
        GenomicSet<Gene> genes = this.exampleData.getGenes();
        GenomicSet<GenomicElement> enhancer = this.exampleData.getEnhancer();
        GenomicSet<GenomicElement> boundaries = this.exampleData.getBoundaries();
        GenomicSet<GenomicElement> domains = this.exampleData.getDomains();
        PhenotypeData phenotypeData = this.exampleData.getPhenotypeData();
        HashMap<Term, HashSet<String>> term2genes = this.exampleData.getTargetTerm2targetGene();
        int regionSize = 20;
        AnnotateCNVs.boundaryOverlap(cnvs, boundaries);
        AnnotateCNVs.defineAdjacentRegionsByDistance(cnvs, regionSize);
        AnnotateCNVs.annoateOverlap(cnvs, boundaries, genes, enhancer, phenotypeData);
        AnnotateCNVs.annoateAdjacentRegions(cnvs, genes, enhancer, phenotypeData);
        InterpretCNVs.annotateEAlowG(cnvs, term2genes);
        assertEquals("EAlowG", cnvs.get("cnv1").getEffectMechanism("EAlowG"));
    }

    /**
     * Test of tandemDuplicationEnhancerAdoption method, of class AnnotateCNVs.
     */
    public void testTandemDuplicationEnhancerAdoption() {
        System.out.println("tandemDuplicationEnhancerAdoption");
        GenomicSet<GenomicElement> domains = this.exampleData.getDomains();
        GenomicSet<GenomicElement> boundaries = this.exampleData.getBoundaries();
        GenomicSet<CNV> dups = this.exampleData.getDups();
        GenomicSet<Gene> genes = this.exampleData.getGenes();
        GenomicSet<GenomicElement> enhancers = this.exampleData.getEnhancer();
        PhenotypeData phenotypeData = this.exampleData.getPhenotypeData();
        AnnotateCNVs.boundaryOverlap(dups, boundaries);
        AnnotateCNVs.annoateOverlap(dups, boundaries, genes, enhancers, phenotypeData);
        AnnotateCNVs.overlapPhenogramScore(dups, phenotypeData);
        AnnotateCNVs.defineOverlappedDomainRegions(dups, domains);
        InterpretCNVs.tandemDuplicationEnhancerAdoption(dups, genes, enhancers, phenotypeData);
        System.out.println(dups.get("dup1"));
        System.out.println(dups.get("dup1").getBoundaryOverlap());
        System.out.println(dups.get("dup1").getLeftOverlappedDomainRegion());
        System.out.println(dups.get("dup1").getRightOverlappedDomainRegion());
        assertEquals("TanDupEA", dups.get("dup1").getEffectMechanism("TanDupEA"));
        System.out.println(dups.get("dup2"));
        System.out.println(dups.get("dup2").getBoundaryOverlap());
        System.out.println(dups.get("dup2").getLeftOverlappedDomainRegion());
        System.out.println(dups.get("dup2").getRightOverlappedDomainRegion());
        System.out.println(dups.get("dup2").getOverlapPhenogramScore());
        assertEquals("onlyGDE", dups.get("dup2").getEffectMechanism("TanDupEA"));
    }

    /**
     * Test of annotateTDBD method, of class AnnotateCNVs.
     */
    public void testAnnotateTDBD() throws IOException {
        System.out.println("annotateTDBD");
        this.exampleData = new ExampleData();
        GenomicSet<CNV> cnvs = this.exampleData.getCnvs();
        GenomicSet<Gene> genes = this.exampleData.getGenes();
        GenomicSet<GenomicElement> enhancer = this.exampleData.getEnhancer();
        GenomicSet<GenomicElement> boundaries = this.exampleData.getBoundaries();
        GenomicSet<GenomicElement> domains = this.exampleData.getDomains();
        PhenotypeData phenotypeData = this.exampleData.getPhenotypeData();
        HashMap<Term, HashSet<String>> terms2genes = this.exampleData.getTargetTerm2targetGene();
        AnnotateCNVs.boundaryOverlap(cnvs, boundaries);
        AnnotateCNVs.defineAdjacentRegionsByDomains(cnvs, domains);
        AnnotateCNVs.annotateOverlappedGenes(cnvs, genes);
        AnnotateCNVs.annotateAdjacentGenes(cnvs, genes);
        AnnotateCNVs.phenogramScore(cnvs, phenotypeData);
        AnnotateCNVs.annotateAdjacentEnhancers(cnvs, enhancer);
        InterpretCNVs.annotateTDBD(cnvs, terms2genes);
        assertEquals("TDBD", cnvs.get("cnv1").getEffectMechanism("TDBD"));
        assertEquals("GDE", cnvs.get("cnv2").getEffectMechanism("TDBD"));
        assertEquals("TDBD", cnvs.get("cnv3").getEffectMechanism("TDBD"));
        assertEquals("NoData", cnvs.get("cnv4").getEffectMechanism("TDBD"));
    }
    /**
     * Test of annotateTDBDjustByScore method, of class AnnotateCNVs.
     */
    //    @Test
    //    public void testAnnotateTDBDjustByScore() {
    //        System.out.println("annotateTDBDjustByScore");
    //        GenomicSet<CNV> inversions = exampleData.getCnvs();
    //        GenomicSet<Gene> genes = exampleData.getGenes();
    //        GenomicSet<GenomicElement> enhancer = exampleData.getEnhancer();
    //        GenomicSet<GenomicElement> boundaries = exampleData.getBoundaries();
    //        GenomicSet<GenomicElement> domains = exampleData.getDomains();
    //        PhenotypeData phenotypeData = exampleData.getPhenotypeData();
    //
    //
    //        AnnotateCNVs.defineAdjacentRegionsByDomains(inversions, domains);
    //        AnnotateCNVs.boundaryOverlap(inversions, boundaries);
    //        AnnotateCNVs.annotateOverlappedGenes(inversions, genes);
    //        AnnotateCNVs.annotateAdjacentGenes(inversions, genes);
    //        AnnotateCNVs.phenogramScore(inversions, phenotypeData);
    //
    //        AnnotateCNVs.annotateTDBDjustByScore(inversions, enhancer);
    //
    //        assertTrue(inversions.get("cnv1").isTDBD());
    //        assertFalse(inversions.get("cnv2").isTDBD());
    //        assertFalse(inversions.get("cnv3").isTDBD());
    //        assertFalse(inversions.get("cnv4").isTDBD());
    //    }

    /**
     * Test of inversionEnhancerAdoption method, of class AnnotateCNVs.
     */
    public void testInversionEnhancerAdoption() {
        System.out.println("inversionEnhancerAdoption");
        GenomicSet<CNV> inversions = this.exampleData.getInvs();
        System.out.println("DEBUG" + inversions);
        GenomicSet<Gene> genes = this.exampleData.getGenes();
        GenomicSet<GenomicElement> enhancer = this.exampleData.getEnhancer();
        System.out.println("DEBUG" + enhancer);
        GenomicSet<GenomicElement> boundaries = this.exampleData.getBoundaries();
        GenomicSet<GenomicElement> domains = this.exampleData.getDomains();
        PhenotypeData phenotypeData = this.exampleData.getPhenotypeData();
        AnnotateCNVs.boundaryOverlap(inversions, boundaries);
        AnnotateCNVs.defineAdjacentRegionsByDomains(inversions, domains);
        AnnotateCNVs.defineOverlappedDomainRegions(inversions, domains);
        AnnotateCNVs.annoateOverlap(inversions, boundaries, genes, enhancer, phenotypeData);
        AnnotateCNVs.annoateAdjacentRegions(inversions, genes, enhancer, phenotypeData);
        System.out.println("DEBUG inv2 left reg:" + inversions.get("inv2").getLeftAdjacentRegion());
        System.out.println("DEBUG inv2 left e:" + inversions.get("inv2").getEnhancersInLeftRegion());
        InterpretCNVs.inversionEnhancerAdoption(inversions, genes, enhancer, phenotypeData);
        System.out.print("DEBUG invs:" + inversions);
        assertEquals("EnhancerInvEA", inversions.get("inv1").getEffectMechanism("InvEA"));
        assertEquals("GeneInvEA", inversions.get("inv2").getEffectMechanism("InvEA"));
        assertEquals("noInvEA", inversions.get("inv3").getEffectMechanism("InvEA"));
    }

    /**
     * Test of annotateTDBDjustByScore method, of class AnnotateCNVs.
     */
    public void testAnnotateTDBDjustByScore() {
        System.out.println("annotateTDBDjustByScore");
        GenomicSet<CNV> cnvs = this.exampleData.getCnvs();
        GenomicSet<Gene> genes = this.exampleData.getGenes();
        GenomicSet<GenomicElement> enhancer = this.exampleData.getEnhancer();
        GenomicSet<GenomicElement> boundaries = this.exampleData.getBoundaries();
        GenomicSet<GenomicElement> domains = this.exampleData.getDomains();
        PhenotypeData phenotypeData = this.exampleData.getPhenotypeData();
        AnnotateCNVs.boundaryOverlap(cnvs, boundaries);
        AnnotateCNVs.defineAdjacentRegionsByDomains(cnvs, domains);
        AnnotateCNVs.annoateOverlap(cnvs, boundaries, genes, enhancer, phenotypeData);
        AnnotateCNVs.annoateAdjacentRegions(cnvs, genes, enhancer, phenotypeData);
        InterpretCNVs.annotateTDBDjustByScore(cnvs);
        assertEquals("TDBD", cnvs.get("cnv1").getEffectMechanism("newTDBD"));
        assertEquals("GDE", cnvs.get("cnv2").getEffectMechanism("newTDBD"));
        assertEquals("TDBD", cnvs.get("cnv3").getEffectMechanism("newTDBD"));
        assertEquals("NoData", cnvs.get("cnv4").getEffectMechanism("newTDBD"));
    }


    
}
