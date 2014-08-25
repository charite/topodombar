/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package io;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;


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
    
    /**
     * convertes an {@link Integer} array into an {@code int} array.
     * @param integerArray
     * @return 
     */
    public static int [] toIntArray(Integer [] integerArray){
        
        int [] intArray = new int [integerArray.length];
        for (int idx=0; idx < integerArray.length; idx++) {
            intArray[idx] = integerArray[idx];
        }
        
        return intArray;
    }


    /**
     * convertes an Integer array into an double array.
     * @param integerArray
     * @return 
     */
    public static double [] toDoubleArray(Integer [] integerArray){
        
        double [] doubleArray = new double [integerArray.length];
        for (int idx=0; idx < integerArray.length; idx++) {
            
            doubleArray[idx] = (double) integerArray[idx];

        }
        
        return doubleArray;
    }

    /**
     * Calculates the fraction of values that are greater or equal to a given
     * value x.
     * 
     * @param values
     * @param x
     * @return 
     */
    public static double fractionGraterEqualX (double [] values, double x){
        
        int n = values.length;
        int counter = 0;
        for (double v : values){
            // increatse counter if values are greater or equal
            if ( v >= x){
                counter++;
            }
        }
        
        // calc frac of values 
        return counter / (double) n; 
        
    }
}

