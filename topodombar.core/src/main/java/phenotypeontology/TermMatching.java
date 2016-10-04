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
 * Implements a matching between terms of a patient and terms associated to a 
 * gene together with non-zero phonomatch scores.
 * 
 * @author ibnsalem
 */
public class TermMatching {

    private final ArrayList<Term> patientTerms;
    private final ArrayList<Term> geneTerms;
    private final ArrayList<Double> score;

    /**
     * Constructor.
     */
    public TermMatching() {
        this.patientTerms = new ArrayList<>();
        this.geneTerms = new ArrayList<>();
        this.score = new ArrayList<>();
    }
    
    
    /**
     * adds a pair of terms with score.
     */
    public void addTermPair(TermPair tp){
        
        this.patientTerms.add(tp.getPp());
        this.geneTerms.add(tp.getGp());
        this.score.add(tp.getS());
        
    }
    
    public TermPair getMax(){
        
        int maxIdx = argmax(score);
        
        TermPair maxPair = new TermPair(
                this.patientTerms.get(maxIdx), 
                this.geneTerms.get(maxIdx), 
                this.score.get(maxIdx));
        
        return(maxPair);
    }
    
    /**
     * Retunrs a string of three columns: the gene pehnotpyes, the patient phenotpye, and the scores between them.
     * Each separated by ';'. Columns are separated by '\t'.
     * @return string
     */
    public String getOutputColumns(){


        //convert phenotpye terms to Strings
        ArrayList<String> pPstr = new ArrayList<>();
        for (Term t : this.patientTerms){
            pPstr.add(t.getIDAsString()); 
        }        

        //convert phenotpye terms to Strings
        ArrayList<String> gPstr = new ArrayList<>();
        for (Term t : this.geneTerms){
            gPstr.add(t.getIDAsString()); 
        }        

        ArrayList<String> sStr = new ArrayList<>();
        for (Double s : this.score){
            sStr.add(Utils.roundToString(s));
        }
        
        String pPstrCol = StringUtils.join(pPstr, ';');
        String gPstrCol = StringUtils.join(gPstr, ';');
        String sStrCol = StringUtils.join(sStr, ';');
        
        
        return pPstrCol + "\t" + gPstrCol + "\t" + sStrCol;
    }

    /**
     * Returns the (first) index with the maximum currMax
     * @param score
     * @return 
     */
    private int argmax(ArrayList<Double> score) {
        
        int idx=0;
        Double currMax = score.get(0);
        
        for (int i=1; i < score.size(); i++){
            if (score.get(i) > currMax){
                idx = i;
            }
            
        }
        
        return(idx);
    }
}
