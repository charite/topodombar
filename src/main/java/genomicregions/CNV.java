/*
 * Copyright (c) 2014, Jonas Ibn-Salem <ibnsalem@molgen.mpg.de>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

package genomicregions;

import java.util.ArrayList;
import java.util.List;
import java.util.HashSet;
import ontologizer.go.Term;
import org.apache.commons.lang3.StringUtils; // provides a join(iterable, char) function

public class CNV extends GenomicElement {
    
    /**
     * Type of CNV ("loss" or "gain"). This field can be later used to indicate 
     * more complex structural variations.
     */
    private String type;
    /**
     * List of phenotypes as HPO term IDs.
     */
    private List<String> phenotpyes;
    
    /**
     * Phenotype terms of the CNV carrier as {@link HashSet} of {@link Term} objects.
     */
    private HashSet<Term> phenotypeTerms;

    
    /**
     * Target term or phenotype category  as single general HPO term ID.
     */    
    private String targetTerm;
    
    // annotations:
    
    /**
     * List of overlapping boundaries.
     */
    private GenomicSet<GenomicElement> boundaryOverlap;
    
    /**
     * true if CNV overlaps a boundary element.
     */
    private boolean hasBoundaryOverlap;
    
    /**
     * List of overlapping genes (any overlap)
     */
    private GenomicSet<Gene> geneOverlap;
    
    /** 
     * Phenogram score of all genes overlapped by the CNV.
     */
    private Double overlapPhenogramScore;
    
    /**
     * Constructor for CNV object.
     * Construct a {@link GenomicElement} and sets all {@link CNV} specific 
     * annotations to default values.
     * 
     * @param chr   Chromosome identifier
     * @param start Start coordinate (0-based)
     * @param end   End coordinate (0-based, half-open)
     * @param name  Name or ID of the CNV
     * 
     * @throws IllegalArgumentException 
     */
    public CNV(String chr, int start, int end, String name) throws IllegalArgumentException {
        super(chr, start, end, name);
        
        this.overlapPhenogramScore = -1.0;
        this.geneOverlap = new GenomicSet();
        this.boundaryOverlap = new GenomicSet();
        
        // set default values for annotations
        type = ".";
        phenotpyes = new ArrayList();
        targetTerm = ".";
    }
    
    /**
     * Constructor for {@link CNV} object.
     * Construct an CNV object by taking all annotations as arguments.
     * 
     * @param chr   Chromosome identifier
     * @param start Start coordinate (0-based)
     * @param end   End coordinate (0-based, half-open)
     * @param name  Name or ID of the CNV
     * @param type  CNV type (loss or gain)
     * @param phenotypes    List of HPO term IDs that represent the phenotypes used to annotate the patient carrying the CNV.
     * @param targetTerm A unspecific target term as a HPO term ID that is used to group patients in cohort or associate patients to tissues for which enhancer data his available.
     */
    public CNV(String chr, int start, int end, String name, String type, List<String> phenotypes, String targetTerm){

        // consturct an CVN object using the constructor of the {@link GenomicElement} super calss
        super(chr, start, end, name);
        
        this.overlapPhenogramScore = -1.0;
        this.geneOverlap = new GenomicSet();
        this.boundaryOverlap = new GenomicSet();
        
        // add annotations
        this.type = type;
        this.phenotpyes = phenotypes;
        this.targetTerm = targetTerm;        
        
    }
    
    /**
     * This function constructs {@link String} that represents an output line
     * for a TAB separated file like BED files.
     * 
     * @return a TAB-separated output line to write BED like files.
     */
    @Override
    public String toOutputLine(){
        // For columns with multiple elements, separate them by semiclon ';'
        String phenotypeCol = StringUtils.join(getPhenotpyes(), ';');
        String boundaryOverlapCol = getBoundaryOverlap().allNamesAsString();
        String geneOverlapCol = getGeneOverlap().allNamesAsString();
        String overlapPhenogramScoreCol = (getOverlapPhenogramScore() != -1) ? getOverlapPhenogramScore().toString() : ".";
        
        // return generic line (chr, start, end, name) and the additional specific columns:
        return super.toOutputLine()
                + "\t" 
                + StringUtils.join(new String[]{
                    getType(), 
                    phenotypeCol, getTargetTerm(),
                    boundaryOverlapCol,
                    geneOverlapCol,
                    overlapPhenogramScoreCol
                }, '\t');
    }

    /**
     * This functions returns a header line for a TAB-separated output file.
     * It contains the head labels for all columns written by the function 
     * {@link toOutputLine} and is specific for the {@link CNV} class.
     * 
     * @return header line for tab separated output file 
     */
    @Override
    public String getOutputHeaderLine(){
        
        // return generic line (chr, start, end, name) and the additional specific columns:
        return super.getOutputHeaderLine()
                + "\t" 
                + StringUtils.join(new String[]{
                    "type", 
                    "phenotypes", 
                    "targetTerm",
                    "boundaryOverlap",
                    "geneOverlap",
                    "overlapPhenogramScore"
                }, '\t');
    }

    /**
     * Type of CNV ("loss" or "gain"). This field can be later used to indicate 
     * more complex structural variations.
     * 
     * @return the type
     */
    public String getType() {
        return type;
    }

    /**
     * Type of CNV ("loss" or "gain"). This field can be later used to indicate 
     * more complex structural variations.
     * 
     * @param type the type to set
     */
    public void setType(String type) {
        this.type = type;
    }

    /**
     * List of phenotypes as HPO term IDs.
     * @return the phenotpyes
     */
    public List<String> getPhenotpyes() {
        return phenotpyes;
    }

    /**
     * List of phenotypes as HPO term IDs.
     * @param phenotpyes the phenotpyes to set
     */
    public void setPhenotpyes(List<String> phenotpyes) {
        this.phenotpyes = phenotpyes;
    }

    /**
     * Phenotype terms of the {@link CNV} carrier as {@link HashSet} of {@link Term} objects.
     * @return the phenotypeTerms
     */
    public HashSet<Term> getPhenotypeTerms() {
        return phenotypeTerms;
    }

    /**
     * Phenotype terms of the {@link CNV} carrier as {@link HashSet} of {@link Term} objects.
     * @param phenotypeTerms the phenotypeTerms to set
     */
    public void setPhenotypeTerms(HashSet<Term> phenotypeTerms) {
        this.phenotypeTerms = phenotypeTerms;
    }

    /**
     * Target term or phenotype category as single general HPO term ID
     * @return the targetTerm
     */
    public String getTargetTerm() {
        return targetTerm;
    }

    /**
     * Target term or phenotype category as single general HPO term ID
     * @param targetTerm the targetTerm to set
     */
    public void setTargetTerm(String targetTerm) {
        this.targetTerm = targetTerm;
    }

    /**
     * {@link GenomicSet} of overlapping boundaries.
     * @return the boundaryOverlap
     */
    public GenomicSet<GenomicElement> getBoundaryOverlap() {
        return boundaryOverlap;
    }

    /**
     * {@link GenomicSet} of overlapping boundaries.
     * @param boundaryOverlap the boundaryOverlap to set
     */
    public void setBoundaryOverlap(GenomicSet<GenomicElement> boundaryOverlap) {
        this.boundaryOverlap = boundaryOverlap;
    }

    /**
     * Should be {@code true} if CNV overlaps a boundary element.
     * @return the hasBoundaryOverlap
     */
    public boolean isHasBoundaryOverlap() {
        return hasBoundaryOverlap;
    }

    /**
     * Should be {@code true} if CNV overlaps a boundary element.
     * @param hasBoundaryOverlap the hasBoundaryOverlap to set
     */
    public void setHasBoundaryOverlap(boolean hasBoundaryOverlap) {
        this.hasBoundaryOverlap = hasBoundaryOverlap;
    }

    /**
     * {@link GenomicSet} of overlapping {@link Gene}s (any overlap)
     * @return the geneOverlap
     */
    public GenomicSet<Gene> getGeneOverlap() {
        return geneOverlap;
    }

    /**
     * {@link GenomicSet} of overlapping {@link Gene}s (any overlap)
     * @param geneOverlap the geneOverlap to set
     */
    public void setGeneOverlap(GenomicSet<Gene> geneOverlap) {
        this.geneOverlap = geneOverlap;
    }

    /**
     * Phenogram score of all genes overlapped by the CNV.
     * @return the overlapPhenogramScore
     */
    public Double getOverlapPhenogramScore() {
        return overlapPhenogramScore;
    }

    /**
     * Phenogram score of all genes overlapped by the CNV.
     * @param overlapPhenogramScore the overlapPhenogramScore to set
     */
    public void setOverlapPhenogramScore(Double overlapPhenogramScore) {
        this.overlapPhenogramScore = overlapPhenogramScore;
    }

    /**
     * Add a phenotype {@link Term} to the set of {@Term}s.
     * 
     * @param t 
     */
    public void addPhenotypeTerm(Term t) {
        this.phenotypeTerms.add(t);
    }

}
