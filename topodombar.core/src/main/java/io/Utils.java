/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package io;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;


/**
 *
 * @author jonas
 */
public class Utils {
    
    
    private static DecimalFormatSymbols dfs = new DecimalFormatSymbols();
    private static DecimalFormat df = new DecimalFormat("0.00");
    
    
    
    public static String roundToString(double d){
        
        // force the decimal separator to be a dot (even in german environemtns)
        dfs.setDecimalSeparator('.');
        df.setDecimalFormatSymbols(dfs);
        
        return df.format(d);
    }

    public static String roundToString(String s){
        
        // force the decimal separator to be a dot (even in german environemtns)
        dfs.setDecimalSeparator('.');
        df.setDecimalFormatSymbols(dfs);
        
        return df.format(s);
    }
}
