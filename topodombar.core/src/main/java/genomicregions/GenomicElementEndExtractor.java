/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package genomicregions;

import de.charite.compbio.jannovar.impl.intervals.IntervalEndExtractor;

/**
 * Implements the getBegin and getEnd functions for use in IntervalArray.
 * @author jonas
 */

public class GenomicElementEndExtractor<T extends GenomicElement> implements IntervalEndExtractor<T> {

        @Override
        public int getBegin(T e) {
                return e.getStart();
        }

        @Override
        public int getEnd(T e) {
                return e.getEnd();
        }

    }

