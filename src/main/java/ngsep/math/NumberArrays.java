/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.math;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class NumberArrays {

	public static byte [] toByteArray (Collection<Byte> values) {
		byte[] ret = new byte[values.size()];
	    int i = 0;
	    for (Byte e : values) ret[i++] = e.byteValue();
	    return ret;
	}
	public static int [] toIntArray (Collection<Integer> values) {
		int[] ret = new int[values.size()];
	    int i = 0;
	    for (Integer e : values) ret[i++] = e.intValue();
	    return ret;
	}
	public static double [] toDoubleArray (Collection<Double> values) {
		double[] ret = new double[values.size()];
	    int i = 0;
	    for (Double d : values) ret[i++] = d.doubleValue();
	    return ret;
	}
	public static List<Integer> toIntegerList(int [] array) {
		List<Integer> ret = new ArrayList<Integer>();
		for(int i=0;i<array.length;i++) ret.add(array[i]);
		return ret;
	}
	public static List<Double> toDoubleList(double [] array) {
		List<Double> ret = new ArrayList<Double>();
		for(int i=0;i<array.length;i++) ret.add(array[i]);
		return ret;
	}
	public static void initializeIntMatrix (int [][] matrix) {
		initializeIntMatrix(matrix, 0);
	}
	public static void initializeIntMatrix (int [][] matrix, int value) {
		for(int i=0;i<matrix.length;i++) {
			for(int j=0;j<matrix[i].length;j++) {
				matrix[i][j] = value;
			}
		}
	}
	public static void initializeDoubleMatrix (double [][] matrix) {
		initializeDoubleMatrix(matrix,0.0);
	}
	public static void initializeDoubleMatrix (double [][] matrix, double value) {
		for(int i=0;i<matrix.length;i++) {
			for(int j=0;j<matrix[i].length;j++) {
				matrix[i][j] = value;
			}
		}
	}
	public static void accumulate(double[][] cumulative, double[][] next) {
		for(int i=0;i<cumulative.length;i++) {
			for(int j=0;j<cumulative[i].length;j++) {
				cumulative[i][j]+=next[i][j];
			}
		}
		
	}
	public static int getIndexMaximum (int [] numbers) {
		return getIndexMaximum(numbers,-1);
	}
	public static int getIndexMaximum (int [] numbers, int ignoreIndex) {
		int idxMax = -1;
		for(int i=0;i<numbers.length;i++) {
			if(i!=ignoreIndex && (idxMax==-1 || numbers[idxMax]<numbers[i])) {
				idxMax = i;
			}
		}
		return idxMax;
	}
	public static int getSum (int [] numbers) {
		int sum = 0;
		for(int i=0;i<numbers.length;i++) sum+=numbers[i];
		return sum;
	}
	
	public static int getIndexMaximum (double [] numbers) {
		return getIndexMaximum(numbers,-1);
	}
	public static int getIndexMaximum (double [] numbers, int ignoreIndex) {
		int idxMax = -1;
		for(int i=0;i<numbers.length;i++) {
			if(i!=ignoreIndex && (idxMax==-1 || numbers[idxMax]<numbers[i])) {
				idxMax = i;
			}
		}
		return idxMax;
	}
	public static double getSum (double [] numbers) {
		double sum = 0;
		for(int i=0;i<numbers.length;i++) sum+=numbers[i];
		return sum;
	}
	public static int getIndexMaximum (byte [] numbers) {
		return getIndexMaximum(numbers,-1);
	}
	public static int getIndexMaximum (byte [] numbers, int ignoreIndex) {
		int idxMax = -1;
		for(int i=0;i<numbers.length;i++) {
			if(i!=ignoreIndex && (idxMax==-1 || numbers[idxMax]<numbers[i])) {
				idxMax = i;
			}
		}
		return idxMax;
	}
}
