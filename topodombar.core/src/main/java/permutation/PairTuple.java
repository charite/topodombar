/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package permutation;

/**
 * A tuple with two elements x and y of (potentially) different types.
 * 
 * @author jonas
 * @param <X> type of first element
 * @param <Y> type of second element
 */
public class PairTuple<X, Y> {
    
    public final X x;
    public final Y y;
    
    public PairTuple(X x, Y y) { 
        this.x = x; 
        this.y = y; 
    }    
    
    
}
