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
   
import annotation.AnnotateCNVs;
import annotation.EffectAnnotation;
import annotation.GeneDosage;
import annotation.InteractionChange;
import io.Utils;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import ontologizer.go.Term;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils; // provides a join(iterable, char) function
import phenotypeontology.PhenotypeData;
import phenotypeontology.TermMatching;

/**
 * This class implements a copy number variation (CNV) object. The CNV has a 
 * defined type and location in a reference genome. Furthermore it is associated to 
 * phenotypes observed in the patient that carries the CNV. Several members of 
 * this object can be used to annotate the CNV with respect to other genomic
 * elements and potential effect mechanisms. 
 * 
 * @author jonas
 */
public class CNV extends GenomicElement {
    
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
    
    /** List of overlapping genes in overlapping TADs */
    private GenomicSet<Gene> genesInOverlapTADs;
    
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
     * A list of {@link GeneDosage} {@link EffectAnnotations} that has this CNV 
     * on phenotypically relevant genes.
     */
    private ArrayList<GeneDosage> geneDosageAnnotations;

    /**
     * A list of {@link InteractionChange} {@link EffectAnnotations} that has this CNV 
     * on phenotypically relevant genes.
     */
    private ArrayList<InteractionChange> interactionChangeAnnotations;
    
    /**
     * possible effect mechanism classes used as key in the {@link effectMechanism} map.
     */
    private final static String [] effectMechanismClasses 
            = new String [] {"TDBD", "newTDBD", "EA", "EAlowG", "TanDupEA", "InvEA"};

    /**
     * Possible CNV effect mechanism classes maped to the possible effect annotations.
     */
    private final  static HashMap<String, String []> effectMechansim2effects;
    static{
        effectMechansim2effects = new HashMap();
        effectMechansim2effects.put("TDBD", new String [] {"TDBD", "Mixed", "GDE", "NoData", "NA"});    
        effectMechansim2effects.put("newTDBD", new String [] {"TDBD", "Mixed", "GDE", "NoData", "NA"});    
        effectMechansim2effects.put("EA", new String [] {"EA", "Mixed", "GDE", "NoData", "NA"});    
        effectMechansim2effects.put("EAlowG", new String [] {"EAlowG", "Mixed", "GDE", "NoData", "NA"});    
        effectMechansim2effects.put("InvEA", new String [] {"EandGInvEA", "EnhancerInvEA", "GeneInvEA", "NoInvEA", "NA"});    
        /** 
         * Inversion enhancer removal (InvER).
         * Here the inversion removes an enhancer that is normally associated to phenotypically relevant gene
         * Or alternatively the gene gets inverted (moves to an other domain).
         */
        effectMechansim2effects.put("InvER", new String [] {"InvertedEnhancer", "InvertedGene", "NoInvER", "NA"});    
        effectMechansim2effects.put("TanDupEA", new String [] {"TanDupEA", "onlyGDE", "NoData", "NA"});    
    }
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
     * Return the the possible annotations for a given effect mechanism class.
     * 
     * @param effectMechanismClass mechanism class for which the possible annotations should be returned
     * @return the the possible annotations for a given effect mechanism class
     */
    public static String [] possibleEeffectAnnotations(String effectMechanismClass) {
        //  check if input is a effect mechanism class:
        // TODO use Enums for the classes an remove the check
        if (! Arrays.asList(CNV.effectMechanismClasses).contains(effectMechanismClass)){
            System.err.printf("[ERROR] The given String '%s' is not a valid effect mechanism class", effectMechanismClass);
            System.exit(1);
        }
        return effectMechansim2effects.get(effectMechanismClass);
    }
    
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
        
        // set empty list for {@link GeneDosage} and {@link InteractionChange} effect annotations
        this.geneDosageAnnotations = new ArrayList<GeneDosage>();
        this.interactionChangeAnnotations = new ArrayList<InteractionChange>();

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
        this.genesInOverlapTADs = new GenomicSet<Gene>();
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

        // set empty list for {@link GeneDosage} and {@link InteractionChange} effect annotations
        this.geneDosageAnnotations = new ArrayList<GeneDosage>();
        this.interactionChangeAnnotations = new ArrayList<InteractionChange>();
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
     * Create a header line as TAB-separated string with all column identifiers 
     * for the simple output file format.
     * 
     * @return simple output file header as TAB separated column descriptions
     */
    public static String getSimpleOutputHeader(){
        
        String [] cnvColumns = new String[]{
                    "#chr",
                    "start",
                    "end",
                    "name",
                    "mostLikelyEffect", 
                    "gene", 
                    "enhancer",
                    "affectedGenes",
                };
        // add the effect mechanism columns and return a TAB-separated string.
        return StringUtils.join(cnvColumns, '\t');
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
     * This function constructs  a ArrayList of {@link String} that represents 
     * output lines for each gene overlapped by this CNV.
     * 
     * @param phenotypeData a {@link PhenotypeData} object to calculate phenoMatch scores
     * @param genes set of genes to consider
     * @return a TAB-separated output line to write BED like files.
     */
    public ArrayList<String> getOverlappedGenesOutputLine(PhenotypeData phenotypeData, GenomicSet<Gene> genes){
        
        ArrayList<String> outLines = new ArrayList<String>();
        
        //convert phenotpye terms to Strings
        HashSet<String> phenotypesIDs = new HashSet<String>();
        for (Term t : this.phenotypes){
            phenotypesIDs.add(t.getIDAsString()); 
        }
        
        // For columns with multiple elements, separate them by semiclon ';'
        String phenotypeCol = StringUtils.join(phenotypesIDs, ';');

        // initialize maximal score per CNV
        Double maxScore = 0.0;
        
        // if CNV does overlap any gene add output line for each gene with score > 0
        if (!genes.isEmpty() ){
                        
            for (Gene g : genes.values()){


                String geneSymbol = g.getSymbol();
                double geneScore = phenotypeData.phenoMatchScore(this.phenotypes, g);
                
                TermMatching termMatching = phenotypeData.phenoMatchScoreWithMatching(this.phenotypes, g);
                
                // only if there is a score larger than zero output the gene
                if (geneScore > 0){
                    
                    String geneScoreStr = Utils.roundToString(geneScore);

                    String [] cnvAnnotations = new String[]{
                            phenotypeCol, 
                            geneSymbol,
                            geneScoreStr,
                            termMatching.getOutputColumns()
                        };

                    // put togeter all annotation string separated by TAB
                    String outLineGene = super.toOutputLine()
                        + "\t" 
                        + StringUtils.join(cnvAnnotations, '\t');            

                    // append to output lines
                    outLines.add(outLineGene);
                    
                    // upate maxScore
                    if (geneScore > maxScore){
                        maxScore = geneScore;
                    }
                }
            }
        }
        
        // if CNV does not overlap any gene or overlapped genes have score 0
        if( maxScore == 0.0){
            String outLineCNV = super.toOutputLine()
                    + "\t" 
                    + StringUtils.join(new String[]{
                        phenotypeCol, 
                        ".",
                        Utils.roundToString(0.0),
                        ".",
                        ".",
                        ".",
                        Utils.roundToString(0.0)
                    }, '\t');
            // append to output lines
            outLines.add(outLineCNV);
        }
        
        return outLines;
        
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
     * {@link GenomicSet} of {@link Gene}s overlapping TADs that overlap the CNV.
     * @return 
     */
    public GenomicSet<Gene> getGenesInOverlapTADs() {
        return genesInOverlapTADs;
    }

    /**
     * {@link GenomicSet} of {@link Gene}s overlapping TADs that overlap the CNV.
     * @param genesInOverlapTADs 
     */
    public void setGenesInOverlapTADs(GenomicSet<Gene> genesInOverlapTADs) {
        this.genesInOverlapTADs = genesInOverlapTADs;
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
     * add an {@link GeneDosage} {@link EffectAnnotation} to this cnv.
     * @param gde 
     */
    public void addGeneDosageAnnotaion(GeneDosage gde){
        this.geneDosageAnnotations.add(gde);
    }

    /**
     * add an {@link InteractionChange} {@link EffectAnnotation} to this cnv.
     * @param ice 
     */
    public void addInteractionChangeAnnotaion(InteractionChange ice){
        this.interactionChangeAnnotations.add(ice);
    }

    /**
     * calculates the most likely effect mechanism form the lists of {@link GeneDosage}
     * and {@link InteractionChange} {@link EffectAnnotation}s and outputs a line
     * for the simple CNV output format.
     * This line is TAB separated and includes the following fields:
     * <ul>
     * <li>chr -> chromosom of CNV</li>
     * <li>start -> start position of CNV</li>
     * <li>end -> end position of CNV</li>
     * <li>name -> the name of the CNV</li>
     * <li>effect -> an indicator of the effect mechanism</li> 
     * <li>gene -> the gene symbols of the involved gene</li>
     * <li>enhancer -> one enhancer elements involved (only one is reported, in case of several possible)</li>
     * <li>possibleEffects -> the number of possible {@link EffectAnnotation} that also might explain the phenotype </li>
     * </ul>
     * @return 
     */
    public String getMostLikelyEffectOutputLine(){
        
        String outLine = super.toOutputLine();
        String effect = null;
        String geneSymbol = ".";
        String enhancerStr = ".";
        Integer possibleEffects = this.geneDosageAnnotations.size() + this.interactionChangeAnnotations.size();

        // find best {@link GeneDosage} annotation
        Double gdeScore = 0.0;
        GeneDosage bestGeneDosage = null;
        
        for (GeneDosage gde: this.geneDosageAnnotations){
            
            if (gde.getScore() > gdeScore){
                bestGeneDosage = gde;
                gdeScore = gde.getScore();
            }
        }
        
        // find best {@link InteractionChange} annotation
        Double interactionScore = 0.0;
        InteractionChange bestInteraction = null;
        
        for (InteractionChange ice : this.interactionChangeAnnotations){
            
            if (ice.getScore() > interactionScore){
                bestInteraction = ice;
                interactionScore = ice.getScore();
            }
        }
        
        // if no explaination was found at all
        if (possibleEffects == 0){
            effect = "noData";
            
        }else{
            // check if GDE or interaction is the better explaination
            if (gdeScore >= interactionScore){

                effect = "GDE_" + bestGeneDosage.getDosageChange();
                geneSymbol = bestGeneDosage.getGene().getSymbol();
                enhancerStr = ".";

                // if interaction score is greater
            }else{

                effect = "interaction_" + bestInteraction.getInteractionChange();
                geneSymbol = bestInteraction.getGene().getSymbol();
                enhancerStr = bestInteraction.getEnhancers().values().iterator().next().getName();
            }
        }
        
        // add fields to output line
        outLine += "\t" + 
            StringUtils.join(
                new String[]{effect, geneSymbol, enhancerStr, possibleEffects.toString(),
                // TODO remove this column after debug
                this.effectMechanism.get("newTDBD")}, 
                '\t');
        
        return outLine;
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
    
    public void debugPrint(){
        System.out.println("DEBUG CNV: " + this.toString());
        System.out.println("DEBUG CNV type " + this.type);
        
        System.out.println("DEBUG CNV phenotype: " + this.phenotypes);
        System.out.println("DEBUG CNV targetTerm: " + this.targetTerm);
        
        System.out.println("DEBUG CNV boundaryOverlap: " + this.boundaryOverlap);
        System.out.println("DEBUG CNV genesInOverlap: " + this.genesInOverlap);
        System.out.println("DEBUG CNV overlapPhenogramScore: " + this.overlapPhenogramScore);
        
        System.out.println("DEBUG CNV leftAdjacentRegion: " + this.leftAdjacentRegion);
        System.out.println("DEBUG CNV enhancersInLeftRegion: " + this.enhancersInLeftRegion);
        System.out.println("DEBUG CNV genesInLeftRegion: " + this.genesInLeftRegion);
        System.out.println("DEBUG CNV leftAdjacentPhenogramScore: " + this.leftAdjacentPhenogramScore);
        
        System.out.println("DEBUG CNV rightAdjacentRegion: " + this.rightAdjacentRegion);
        System.out.println("DEBUG CNV enhancersInRightRegion: " + this.enhancersInRightRegion);
        System.out.println("DEBUG CNV genesInRightRegion: " + this.genesInRightRegion);
        System.out.println("DEBUG CNV rightAdjacentPhenogramScore: " + this.rightAdjacentPhenogramScore);
        
        System.out.println("DEBUG CNV leftOverlappedDomainRegion: " + this.leftOverlappedDomainRegion);
        System.out.println("DEBUG CNV rightOverlappedDomainRegion: " + this.rightOverlappedDomainRegion);
        
        System.out.println("DEBUG CNV effectMechanism: " + this.effectMechanism);
    }
}
