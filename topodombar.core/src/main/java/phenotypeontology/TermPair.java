/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phenotypeontology;

import java.util.Comparator;
import ontologizer.go.Term;

/**
 * Implements a pair of two single terms with a score
 * @author jonas
 */
public class TermPair {
    
    private final Term pp;
    private final Term gp;
    private final Double s;
    private final Term lca; // lowest common ancestor 
    
    /**
     * Constructor
     * @param pp
     * @param gp
     * @param s
     */
    public TermPair(Term pp, Term gp, Double s, Term lca){
        
        this.pp = pp;
        this.gp = gp;
        this.s = s;
        this.lca = lca;
    }

    public Term getPp() {
        return pp;
    }

    public Term getGp() {
        return gp;
    }

    public Double getS() {
        return s;
    }

    public Term getLca() {
        return lca;
    }
    
    public String toString(){
        return pp.toString() + ":" + gp.toString() + ":" + lca.toString() + ":" + s.toString();
    }
    
    /**
     * A comparison function, which imposes a total ordering on some collection 
     * of {@link TermPair} objects by there phenotype similarity score.
     * 
     * See {@see Comparator}.
     */
    private static class termPairComparator implements Comparator<TermPair> {


        @Override
        public int compare(TermPair tp1, TermPair tp2) {
            
            // compare the the order ranks of the effect mechanisms  
            return Integer.signum(tp1.s.compareTo(tp2.s));
        }
    }

    /**
     * A Comparator that compares {@link CVN} objects by there 
     * TDBD effect mechanism annotation.
     * This comparator can be used to sort CNVs in the output file.
     */
    public static final Comparator<TermPair> TERM_PAIR_SCORE_ORDER = new TermPair.termPairComparator();
    
    
}
