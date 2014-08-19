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
   
import io.Utils;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import ontologizer.go.Term;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils; // provides a join(iterable, char) function
import phenotypeontology.PhenotypeData;

public class CNV extends GenomicElement {

    /**
     * A Comparator that compares {@link CVN} objects by there 
     * TDBD effect mechanism annotation.
     * This comparator can be used to sort CNVs in the output file.
     */
    public static final Comparator<CNV> EFFECTMECHANISM_TDBD_ORDER = new CNV.effectMechanismTdbdComparator();

    /**
     * possible effect mechanism classes used as key in the {@link effectMechanism} map.
     * @return the effectMechanismClasses
     */
    public static String[] getEffectMechanismClasses() {
        return effectMechanismClasses;
    }
    
    /**
     * Type of CNV ("loss", "gain", "inversion"). This field can be later used to indicate 
     * more complex structural variations.
     */
    private final String type;
    
    /**
     * Phenotype terms of the CNV carrier as {@link HashSet} of {@link Term} objects.
     */
    private final HashSet<Term> phenotypes;

    
    /** Target term or phenotype category as single general HPO term ID. */    
    private final Term targetTerm;
    
    
    /** List of overlapping boundaries. */
    private GenomicSet<GenomicElement> boundaryOverlap;    
   
    /** List of overlapping genes (any overlap) */
    private GenomicSet<Gene> genesInOverlap;
    
    /** Adjacent genomic region on the left (5') site of the CNV. */
    private GenomicElement leftAdjacentRegion;
    
    /** Adjacent genomic region on the right (3') site of the CNV. */
    private GenomicElement rightAdjacentRegion;
    
    /** {@link GenomicSet} of {@link Genes} in the left adjacent region of this CNV. */
    private GenomicSet<Gene> genesInLeftRegion;
    
    /** {@link GenomicSet} of {@link Genes} in the right adjacent region of this CNV. */
    private GenomicSet<Gene> genesInRightRegion;
    
    /** {@link GenomicSet} of enhancers in the left adjacent region of this CVN. */
    private GenomicSet<GenomicElement> enhancersInLeftRegion;
    
    /** {@link GenomicSet} of enhancers in the right adjacent region of this CVN. */
    private GenomicSet<GenomicElement> enhancersInRightRegion;

    /** Overlapped genomic region in the domain overlapping the 3' end of the CNV */
    private GenomicElement leftOverlappedDomainRegion;
    
    /** Overlapped genomic region in the domain overlapping the 5' end of the CNV */
    private GenomicElement rightOverlappedDomainRegion;
    
    /** Phenogram score of all genes overlapped by the CNV. */
    private Double overlapPhenogramScore;
    
    /** Phenogram score of genes in the left adjacent region. */
    private Double leftAdjacentPhenogramScore;

    /** Phenogram score of genes in the right adjacent region. */
    private Double rightAdjacentPhenogramScore;
    
    /**
     * Holds the final effect mechanism annotations. For each effect class (like 
     * TDBD, EA, or EAlowG), its holds the mechanism that best explains the patients
     * phenotype, based on the genomic context and phenotypic similartites.<br><br>
     * The key word "TDBD" indicator that this CNV is a topological domain boundary disruption (TDBD),
     * gene dosage effect (GDE), both (Mixed) or not explainable (NoData).<br>
     * "EA"  indicator that this CNV corresponds to the category enhancer adoption (EA),
     * gene dosage effect (GDE), both (Mixed) or not explainable (NoData).<br>
     * "EAlowG" indicator that this CNV corresponds to the category enhancer adoption 
     * with low score in overlap (EAlowG), gene dosage effect (GDE), both (Mixed)
     * or not explainable (NoData).
     */
    private HashMap<String, String> effectMechanism;

    /**
     * possible effect mechanism classes used as key in the {@link effectMechanism} map.
     */
    private final static String [] effectMechanismClasses 
            = new String [] {"TDBD", "newTDBD", "EA", "EAlowG", "TanDupEA", "InvEA"};
    
//    /** indicator that this CNV is a topological domain boundary disruption (TDBD),
//     * gene dosage effect (GDE), both (Mixed) or not explainable (NoData). */
//    private String effectMechanismTDBD = ".";
//    
//    /** indicator that this CNV corresponds to the category enhancer adoption (EA),
//     * gene dosage effect (GDE), both (Mixed) or not explainable (NoData). */
//    private String effectMechanismEA = ".";
//    
//    /** indicator that this CNV corresponds to the category enhancer adoption 
//     * with low score in overlap (EAlowG), gene dosage effect (GDE), both (Mixed)
//     * or not explainable (NoData). */
//    private String effectMechanismEAlowG = ".";
    
    
    
    /**
     * Constructor for CNV object.
     * Construct a {@link GenomicElement} and sets all {@link CNV} specific 
     * annotations to default values.
     * 
     * @param chr   Chromosome identifier
     * @param start Start coordinate (0-based)
     * @param end   End coordinate (0-based, half-open)
     * @param name  Name or ID of the CNV
     * @param type  CNV type (e.g. loss or gain)
     * 
     * @throws IllegalArgumentException 
     */
    public CNV(String chr, int start, int end, String name, String type) throws IllegalArgumentException {
        super(chr, start, end, name);
        
        this.overlapPhenogramScore = -1.0;
        this.genesInOverlap = new GenomicSet();
        this.boundaryOverlap = new GenomicSet();
        
        // set default values for annotations
        this.type = type;
        this.phenotypes = new HashSet<Term>();
        this.targetTerm = null;
        this.boundaryOverlap = new GenomicSet<GenomicElement>();
        this.genesInOverlap = new GenomicSet<Gene>();
        this.overlapPhenogramScore = -1.0;
        this.genesInLeftRegion = new GenomicSet<Gene>();
        this.leftAdjacentPhenogramScore = -1.0;
        this.genesInRightRegion = new GenomicSet<Gene>();
        this.rightAdjacentPhenogramScore = -1.0;
        this.enhancersInLeftRegion = new GenomicSet<GenomicElement>();
        this.enhancersInRightRegion = new GenomicSet<GenomicElement>();
        
        // set default dot "." for effect mechanism class annotations
        this.effectMechanism = new HashMap<String, String>();
        for (String mechanismClass : CNV.effectMechanismClasses){
            this.effectMechanism.put(mechanismClass, ".");
        }
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
    public CNV(String chr, int start, int end, String name, 
                String type, HashSet<Term> phenotypes, Term targetTerm){

        // consturct an CVN object using the constructor of the {@link GenomicElement} super calss
        super(chr, start, end, name);
        
        this.overlapPhenogramScore = -1.0;
        this.genesInOverlap = new GenomicSet();
        this.boundaryOverlap = new GenomicSet();
        
        // add annotations
        this.type = type;
        this.phenotypes = phenotypes;
        this.targetTerm = targetTerm;        
        
        // set default annotaton for other fealds
        this.boundaryOverlap = new GenomicSet<GenomicElement>();
        this.genesInOverlap = new GenomicSet<Gene>();
        this.overlapPhenogramScore = -1.0;
        this.genesInLeftRegion = new GenomicSet<Gene>();
        this.leftAdjacentPhenogramScore = -1.0;
        this.genesInRightRegion = new GenomicSet<Gene>();
        this.rightAdjacentPhenogramScore = -1.0;
        this.enhancersInLeftRegion = new GenomicSet<GenomicElement>();
        this.enhancersInRightRegion = new GenomicSet<GenomicElement>();        

        // set default dot "." for effect mechanism class annotations
        this.effectMechanism = new HashMap<String, String>();
        for (String mechanismClass : CNV.effectMechanismClasses){
            this.effectMechanism.put(mechanismClass, ".");
        }
    }
    
    /**
     * Create a header line as TAB-separated string with all column identifiers 
     * for the output file.
     * 
     * @return output file header as TAB separated column descriptions
     */
    public static String getHeaderLine(){
        
        String [] cnvColumns = new String[]{
                    "#chr",
                    "start",
                    "end",
                    "name",
                    "type", 
                    "phenotypes", 
                    "targetTerm",
                    "boundaryOverlap",
                    "overlapGenes",
                    "overlapScore",
                    "leftAdjacentGenes",
                    "leftAdjacentScore",
                    "rightAdjacentGenes",
                    "rightAdjacentScore",
                    "leftAdjacentEnhancers",
                    "rightAdjacentEnhancers",
                };
        // add the effect mechanism columns and return a TAB-separated string.
        return StringUtils.join(ArrayUtils.addAll(cnvColumns, CNV.getEffectMechanismClasses()), '\t');
    }
    
    /**
     * This function constructs {@link String} that represents an output line
     * for a TAB separated file like BED files.
     * 
     * @return a TAB-separated output line to write BED like files.
     */
    @Override
    public String toOutputLine(){
        
        //convert phenotpye terms to Strings
        HashSet<String> phenotypesIDs = new HashSet<String>();
        for (Term t : this.phenotypes){
            phenotypesIDs.add(t.getIDAsString()); 
        }
        
        // For columns with multiple elements, separate them by semiclon ';'
        String phenotypeCol = StringUtils.join(phenotypesIDs, ';');
        String boundaryOverlapCol = this.boundaryOverlap.allNamesAsString();
        String overlapPhenogramScoreStr = (this.overlapPhenogramScore != -1.0) ? this.overlapPhenogramScore.toString() : ".";
        String leftAdjacentPhenogramScoreStr = (this.leftAdjacentPhenogramScore != -1.0) ? this.leftAdjacentPhenogramScore.toString() : ".";
        String rightAdjacentPhenogramScoreStr = (this.rightAdjacentPhenogramScore != -1.0) ? this.rightAdjacentPhenogramScore.toString() : ".";
        
        // return generic line (chr, start, end, name) and the additional specific columns:
        return super.toOutputLine()
            + "\t" 
            + StringUtils.join(new String[]{
                this.type, 
                phenotypeCol, this.targetTerm.getIDAsString(),
                Integer.toString(this.boundaryOverlap.size()),
                this.genesInOverlap.allNamesAsString(),
                overlapPhenogramScoreStr,
                this.genesInLeftRegion.allNamesAsString(),
                leftAdjacentPhenogramScoreStr,
                this.genesInRightRegion.allNamesAsString(),
                rightAdjacentPhenogramScoreStr,
                this.enhancersInLeftRegion.allNamesAsString(),
                this.enhancersInRightRegion.allNamesAsString(), 
                this.effectMechanism.get("TDBD")
            }, '\t');
            
    }
    
    /**
     * This function constructs {@link String} that represents an output line
     * for a TAB separated file like BED files.
     * 
     * @param phenotypeData a {@link PhenotypeData} object to calculate phenoMatch scores
     * @return a TAB-separated output line to write BED like files.
     */
    public String getOutputLineWithRelevantGenes(PhenotypeData phenotypeData){
        
        
        //convert phenotpye terms to Strings
        HashSet<String> phenotypesIDs = new HashSet<String>();
        for (Term t : this.phenotypes){
            phenotypesIDs.add(t.getIDAsString()); 
        }
        
        // For columns with multiple elements, separate them by semiclon ';'
        String phenotypeCol = StringUtils.join(phenotypesIDs, ';');
        String boundaryOverlapCol = this.boundaryOverlap.allNamesAsString();
        String overlapPhenogramScoreStr = (this.overlapPhenogramScore != -1.0) ? Utils.roundToString(this.overlapPhenogramScore) : ".";
        String leftAdjacentPhenogramScoreStr = (this.leftAdjacentPhenogramScore != -1.0) ? Utils.roundToString(this.leftAdjacentPhenogramScore) : ".";
        String rightAdjacentPhenogramScoreStr = (this.rightAdjacentPhenogramScore != -1.0) ? Utils.roundToString(this.rightAdjacentPhenogramScore) : ".";
        
        // get GeneSymbol and score for all 
        ArrayList<String> overlapTargetGenes = new ArrayList<String>();
        for (Gene g : this.genesInOverlap.values()){
            double geneScore = phenotypeData.phenoMatchScore(this.phenotypes, g);
            
            // only if the gene is associated with phenotype from the patient
            // (that is if phenomatch score is > 0), add it to the output set.
            if (geneScore > 0){
                overlapTargetGenes.add( g.getSymbol() + ":" + Utils.roundToString(geneScore));
            }
        }
        
        ArrayList<String> leftTargetGenes = new ArrayList<String>();
        for (Gene g : this.genesInLeftRegion.values()){
            double geneScore = phenotypeData.phenoMatchScore(this.phenotypes, g);
            
            // only if the gene is associated with phenotype from the patient
            // (that is if phenomatch score is > 0), add it to the output set.
            if (geneScore > 0){
                leftTargetGenes.add( g.getSymbol() + ":" + Utils.roundToString(geneScore)); 
            }
        }
        
        ArrayList<String> rightTargetGenes = new ArrayList<String>();
        for (Gene g : this.genesInRightRegion.values()){
            double geneScore = phenotypeData.phenoMatchScore(this.phenotypes, g);
            
            // only if the gene is associated with phenotype from the patient
            // (that is if phenomatch score is > 0), add it to the output set.
            if (geneScore > 0){
                rightTargetGenes.add( g.getSymbol() + ":" + Utils.roundToString(geneScore)); 
            }
        }

        // put together all CNV specific annotations
        String [] cnvAnnotations = new String[]{
                this.type, 
                phenotypeCol, this.targetTerm.getIDAsString(),
                Integer.toString(this.boundaryOverlap.size()),
                overlapTargetGenes.isEmpty() ? "." : StringUtils.join(overlapTargetGenes,';'),
                overlapPhenogramScoreStr,
                leftTargetGenes.isEmpty() ? "." : StringUtils.join(leftTargetGenes, ';') ,
                leftAdjacentPhenogramScoreStr,
                rightTargetGenes.isEmpty() ? "." : StringUtils.join(rightTargetGenes, ';') ,
                rightAdjacentPhenogramScoreStr,
                this.enhancersInLeftRegion.allNamesAsString(),
                this.enhancersInRightRegion.allNamesAsString(), 
            };
        
        // put together all effect mechanism annotation columns
        int nClasses = CNV.getEffectMechanismClasses().length;
        String [] effectMechanismAnnotations = new String [nClasses];
        for (int i=0; i<nClasses; i++){
            effectMechanismAnnotations[i] = this.effectMechanism.get(CNV.getEffectMechanismClasses()[i]);
        }
                
        // return generic line (chr, start, end, name) and the additional 
        // specific columns separated by TAB character:
        return super.toOutputLine()
            + "\t" 
            + StringUtils.join(ArrayUtils.addAll(cnvAnnotations, effectMechanismAnnotations), '\t');
    }
    

    /**
     * Type of CNV ("loss", "gain" or "inversion"). This field can be later used to indicate 
     * more complex structural variations.
     * 
     * @return the type
     */
    public String getType() {
        return type;
    }

    /**
     * Phenotype terms of the {@link CNV} carrier as {@link HashSet} of {@link Term} objects.
     * @return 
     */
    public HashSet<Term> getPhenotypes() {
        return this.phenotypes;
    }

    /**
     * Target term or phenotype category as single general HPO term ID
     * @return the targetTerm
     */
    public Term getTargetTerm() {
        return this.targetTerm;
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
    public boolean hasBoundaryOverlap() {
        return ! boundaryOverlap.isEmpty();
    }

    /**
     * {@link GenomicSet} of overlapping {@link Gene}s (any overlap)
     * @return the genesInOverlap
     */
    public GenomicSet<Gene> getGenesInOverlap() {
        return genesInOverlap;
    }

    /**
     * {@link GenomicSet} of overlapping {@link Gene}s (any overlap)
     * @param genesInOverlap the genesInOverlap to set
     */
    public void setGenesInOverlap(GenomicSet<Gene> genesInOverlap) {
        this.genesInOverlap = genesInOverlap;
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
     * This method is just for creating artificial CNV objects for testing.
     * @param t 
     */
    public void addPhenotypeTerm(Term t) {
        this.phenotypes.add(t);
    }

    /**
     * Adjacent genomic region on the left (5') site of the CNV.
     * @return the leftAdjacentRegion
     */
    public GenomicElement getLeftAdjacentRegion() {
        return leftAdjacentRegion;
    }

    /**
     * Adjacent genomic region on the left (5') site of the CNV.
     * @param leftAdjacentRegion the leftAdjacentRegion to set
     */
    public void setLeftAdjacentRegion(GenomicElement leftAdjacentRegion) {
        this.leftAdjacentRegion = leftAdjacentRegion;
    }

    /**
     * Adjacent genomic region on the right (3') site of the CNV.
     * @return the rightAdjacentRegion
     */
    public GenomicElement getRightAdjacentRegion() {
        return rightAdjacentRegion;
    }

    /**
     * Adjacent genomic region on the right (3') site of the CNV.
     * @param rightAdjacentRegion the rightAdjacentRegion to set
     */
    public void setRightAdjacentRegion(GenomicElement rightAdjacentRegion) {
        this.rightAdjacentRegion = rightAdjacentRegion;
    }

    /**
     * Phenogram score of genes in the left adjacent region.
     * @return the leftAdjacentPhenogramScore
     */
    public Double getLeftAdjacentPhenogramScore() {
        return leftAdjacentPhenogramScore;
    }

    /**
     * Phenogram score of genes in the left adjacent region.
     * @param leftAdjacentPhenogramScore the leftAdjacentPhenogramScore to set
     */
    public void setLeftAdjacentPhenogramScore(Double leftAdjacentPhenogramScore) {
        this.leftAdjacentPhenogramScore = leftAdjacentPhenogramScore;
    }

    /**
     * Phenogram score of genes in the right adjacent region.
     * @return the rightAdjacentPhenogramScore
     */
    public Double getRightAdjacentPhenogramScore() {
        return rightAdjacentPhenogramScore;
    }

    /**
     * Phenogram score of genes in the right adjacent region.
     * @param rightAdjacentPhenogramScore the rightAdjacentPhenogramScore to set
     */
    public void setRightAdjacentPhenogramScore(Double rightAdjacentPhenogramScore) {
        this.rightAdjacentPhenogramScore = rightAdjacentPhenogramScore;
    }

    /**
     * {@link GenomicSet} of {@link Genes} in the left adjacent region of this CNV.
     * @return the genesInLeftRegion
     */
    public GenomicSet<Gene> getGenesInLeftRegion() {
        return genesInLeftRegion;
    }

    /**
     * {@link GenomicSet} of {@link Genes} in the left adjacent region of this CNV.
     * @param genesInLeftRegion the genesInLeftRegion to set
     */
    public void setGenesInLeftRegion(GenomicSet<Gene> genesInLeftRegion) {
        this.genesInLeftRegion = genesInLeftRegion;
    }

    /**
     * {@link GenomicSet} of {@link Genes} in the right adjacent region of this CNV.
     * @return the genesInRightRegion
     */
    public GenomicSet<Gene> getGenesInRightRegion() {
        return genesInRightRegion;
    }

    /**
     * {@link GenomicSet} of {@link Genes} in the right adjacent region of this CNV.
     * @param genesInRightRegion the genesInRightRegion to set
     */
    public void setGenesInRightRegion(GenomicSet<Gene> genesInRightRegion) {
        this.genesInRightRegion = genesInRightRegion;
    }

    /**
     * {@link GenomicSet} of enhancers in the left adjacent region of this CVN.
     * @return the enhancersInLeftRegion
     */
    public GenomicSet<GenomicElement> getEnhancersInLeftRegion() {
        return enhancersInLeftRegion;
    }

    /**
     * {@link GenomicSet} of enhancers in the left adjacent region of this CVN.
     * @param enhancersInLeftRegion the enhancersInLeftRegion to set
     */
    public void setEnhancersInLeftRegion(GenomicSet<GenomicElement> enhancersInLeftRegion) {
        this.enhancersInLeftRegion = enhancersInLeftRegion;
    }

    /**
     * {@link GenomicSet} of enhancers in the right adjacent region of this CVN.
     * @return the enhancersInRightRegion
     */
    public GenomicSet<GenomicElement> getEnhancersInRightRegion() {
        return enhancersInRightRegion;
    }

    /**
     * {@link GenomicSet} of enhancers in the right adjacent region of this CVN.
     * @param enhancersInRightRegion the enhancersInRightRegion to set
     */
    public void setEnhancersInRightRegion(GenomicSet<GenomicElement> enhancersInRightRegion) {
        this.enhancersInRightRegion = enhancersInRightRegion;
    }
    
    /**
     * Set the effect mechanism. For each effect class (like 
     * TDBD, EA, or EAlowG), this is the mechanism that best explains the patients
     * phenotype, based on the genomic context and phenotypic similartites.<br><br>
     * The mechanism class "TDBD" tries to explain this CNV is a topological 
     * domain boundary disruption (TDBD),
     * gene dosage effect (GDE), both (Mixed) or not explainable (NoData).<br>
     * The class "EA"  indicates that this CNV corresponds to the category enhancer adoption (EA),
     * gene dosage effect (GDE), both (Mixed) or not explainable (NoData).<br>
     * Finally, the class "EAlowG" indicates that this CNV corresponds to the category enhancer adoption 
     * with low score in overlap (EAlowG), gene dosage effect (GDE), both (Mixed)
     * or not explainable (NoData).
     * 
     * @param mechanimsClass effect mechanism class
     * @param mechanism the most likely effect mechansim for the given class
     */
    public void setEffectMechanism(String mechanimsClass, String mechanism){
        
        // if input mechanimsClass is one of the initialized clases 
        if ( this.effectMechanism.containsKey(mechanimsClass) ){
            // set the annotation
            this.effectMechanism.put(mechanimsClass, mechanism);
        
        // if any not initialized effect mechnism calss is used as input
        }else{
            // thorw an exception and exit
            System.err.printf("[ERROR] Try to set annotation dict with wrong effect class."
                    + "Effect class '%s' is not an supporte effect mechanism. Exit now.", mechanimsClass );
            System.exit(1);
        }
    }
    
    /**
     * Returns the most likely effect mechanism for the given class. For each effect class (like 
     * TDBD, EA, or EAlowG), this is the mechanism that best explains the patients
     * phenotype, based on the genomic context and phenotypic similartites.<br><br>
     * The mechanism class "TDBD" tries to explain this CNV is a topological 
     * domain boundary disruption (TDBD),
     * gene dosage effect (GDE), both (Mixed) or not explainable (NoData).<br>
     * The class "EA"  indicates that this CNV corresponds to the category enhancer adoption (EA),
     * gene dosage effect (GDE), both (Mixed) or not explainable (NoData).<br>
     * Finally, the class "EAlowG" indicates that this CNV corresponds to the category enhancer adoption 
     * with low score in overlap (EAlowG), gene dosage effect (GDE), both (Mixed)
     * or not explainable (NoData).
     * 
     * @param mechanismClass
     * @return the most likely effect mechanism
     */
    public String getEffectMechanism(String mechanismClass){

        return this.effectMechanism.get(mechanismClass);

    }

    /**
     * Overlapped genomic region in the domain overlapping the 3' end of the CNV
     * @return the leftOverlappedDomainRegion
     * @see AnnotateCNVs.defineOverlappedDomainRegions
     */
    public GenomicElement getLeftOverlappedDomainRegion() {
        return leftOverlappedDomainRegion;
    }

    /**
     * Overlapped genomic region in the domain overlapping the 3' end of the CNV
     * @param leftOverlappedDomainRegion the leftOverlappedDomainRegion to set
     * @see AnnotateCNVs.defineOverlappedDomainRegions
     */
    public void setLeftOverlappedDomainRegion(GenomicElement leftOverlappedDomainRegion) {
        this.leftOverlappedDomainRegion = leftOverlappedDomainRegion;
    }

    /**
     * Overlapped genomic region in the domain overlapping the 5' end of the CNV
     * @return the rightOverlappedDomainRegion
     * @see AnnotateCNVs.defineOverlappedDomainRegions
     */
    public GenomicElement getRightOverlappedDomainRegion() {
        return rightOverlappedDomainRegion;
    }

    /**
     * Overlapped genomic region in the domain overlapping the 5' end of the CNV
     * @param rightOverlappedDomainRegion the rightOverlappedDomainRegion to set
     * @see AnnotateCNVs.defineOverlappedDomainRegions
     */
    public void setRightOverlappedDomainRegion(GenomicElement rightOverlappedDomainRegion) {
        this.rightOverlappedDomainRegion = rightOverlappedDomainRegion;
    }

    /**
     * A comparison function, which imposes a total ordering on some collection 
     * of {@link CNV} objects by there TDBD effect mechanism annotation.
     * 
     * See {@see Comparator}.
     */
    private static class effectMechanismTdbdComparator implements Comparator<CNV> {

        /**
         * Holds the rank of the ordering of each effect mechanism.
         */
        private final static HashMap<String, Integer> effectRank;
        
        static{
            effectRank = new HashMap<String, Integer>();
            effectRank.put("TDBD", 0);
            effectRank.put("Mixed", 1);
            effectRank.put("GDE", 2);
            effectRank.put("NoData", 3);
        }
        
        @Override
        public int compare(CNV e1, CNV e2) {
            
            // get order of each effenct mechanism string
            Integer rankE1 = effectRank.get(e1.getEffectMechanism("TDBD"));
            Integer rankE2 = effectRank.get(e2.getEffectMechanism("TDBD"));
            
            // compare the the order ranks of the effect mechanisms  
            return Integer.signum(rankE1.compareTo(rankE2));
        }
    }
}
