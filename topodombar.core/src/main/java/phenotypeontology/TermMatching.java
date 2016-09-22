/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phenotypeontology;

import io.Utils;
import java.util.ArrayList;
import ontologizer.go.Term;
import org.apache.commons.lang3.StringUtils;


/**
 * Implements a matching between terms of a patient and terms associated to a gene together with non-zero phonomatch scores.
 * 
 * @author ibnsalem
 */
public class TermMatching {

    /**
     * Constructor.
     */
    public TermMatching() {
        this.geneTerms = new ArrayList<>();
        this.patientTerms = new ArrayList<>();
        this.score = new ArrayList<>();
    }
    
    private ArrayList<Term> geneTerms;
    private ArrayList<Term> patientTerms;
    private ArrayList<Double> score;
    
    
    
    /**
     * adds a pair of terms with score.
     * @param gt
     * @param pt
     * @param s 
     */
    public void addTermPair(Term gt, Term pt, Double s){
        
        this.geneTerms.add(gt);
        this.patientTerms.add(pt);
        this.score.add(s);
        
    }
    
    /**
     * Retunrs a string of three columns: the gene pehnotpyes, the patient phenotpye, and the scores between them.
     * Each separated by ';'. Columns are separated by '\t'.
     * @return string
     */
    public String getOutputColumns(){

        //convert phenotpye terms to Strings
        ArrayList<String> gPstr = new ArrayList<>();
        for (Term t : this.geneTerms){
            gPstr.add(t.getIDAsString()); 
        }        

        //convert phenotpye terms to Strings
        ArrayList<String> pPstr = new ArrayList<>();
        for (Term t : this.patientTerms){
            pPstr.add(t.getIDAsString()); 
        }        
        
        ArrayList<String> sStr = new ArrayList<>();
        for (Double s : this.score){
            sStr.add(Utils.roundToString(s));
        }
        
        String gPstrCol = StringUtils.join(gPstr, ';');
        String tPstrCol = StringUtils.join(pPstr, ';');
        String sStrCol = StringUtils.join(sStr, ';');
        
        
        return gPstrCol + tPstrCol + sStrCol;
    }
}
