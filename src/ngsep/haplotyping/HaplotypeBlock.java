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
package ngsep.haplotyping;

import java.util.ArrayList;
import java.util.List;

import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.GenomicVariant;

public class HaplotypeBlock {

	/**
	 * Represents the matrix of fragments and variants.
	 */
	private List<HaplotypeFragment> matrix;

	/**
	 * Represents the list of variants.
	 */
	private List<? extends GenomicVariant> variants;
	
	/**
	 * Represents a haplotype.
	 */
	private byte haplotype[];
	
	/**
	 * Constructor that initializes the attributes of a HaplotypeBlock with the given parameters.
	 * @param variants.
	 */
	public HaplotypeBlock(List<? extends GenomicVariant> variants) 
	{
		this.variants = variants;
		matrix = new ArrayList <HaplotypeFragment>();
		haplotype = null;

	}
	/**
	 * Add fragment to the matrix
	 * @param firstColumn where valid allele calls are found
	 * @param alleleCalls Calls starting from the given column
	 */
	public void addFragment(int firstColumn, byte[] alleleCalls) {
		// TODO Auto-generated method stub
		
	}
		
	/**
	 * Returns the call in a given position in the matrix of fragments.
	 * @param i
	 * @param j
	 * @return allele.
	 */
	public byte getAllele(int i, int j)
	{
		HaplotypeFragment row = matrix.get(i);
		byte allele = row.getCall(j);
		return allele;
	}
	
	/**
	 * Returns the haplotype.
	 * @return haplotype.
	 */
	public byte [] getHaplotype()
	{
		return haplotype;
	}
	
	/**
	 * Returns the variant in the given position in the list of variants.
	 * <b> pre: </b> The list of variants has been initialized.
	 * @param i
	 * @return genomic variant.
	 */
	public GenomicVariant getVariant(int i)
	{
		GenomicVariant temp = variants.get(i);
		return  temp;
	}
	
	/**
	 * Returns Hamming distance between two fragments
	 * <b> pre: </b> The matrix of fragments has been initialized.
	 * @param row1. Row1 < Row2
	 * @param row2.
	 * tener en cuenta los maximos
	 * @return Hamming distance between two fragments.
	 */
	public int getHammingDistance(int row1, int row2) 
	{
		int distance = 0;
		for(int i = getFirstColumn(row2) ; i < getLastColumn(row1); i++)
		{
			byte allele1 = getAllele(row1, i);
			byte allele2 = getAllele(row2, i);
			if( allele1 != allele2 && (allele1 != CalledGenomicVariant.ALLELE_UNDECIDED && allele2!= CalledGenomicVariant.ALLELE_UNDECIDED))
			{
				distance ++;
			}
		}
			
		return distance;
	}
	
	/**comparacon cruzada de primeros contra ultimos.
	 * <b> pre: </b> The matrix of fragments has been initialized.
	 * @param row1
	 * @param row2
	 * @return True when the two fragments overlap.
	 */
	public boolean overlap(int row1, int row2)
	{
		boolean overlap = false;
		int initPosRow1 = getFirstColumn(row1);
		int lastPosRow1 = getLastColumn(row1);
		int initPosRow2 = getFirstColumn(row2);
		int lastPosRow2 = getLastColumn(row2);
		int a = Math.min(initPosRow1, initPosRow2);
		int b = Math.min(lastPosRow1, lastPosRow2);
		int c = Math.max(initPosRow1, initPosRow2);
		int d = Math.max(lastPosRow1, lastPosRow2);
		
		if( a < d && b > c )
		{
			overlap = true;
		}
		return overlap;
	}
	
	/**
	 * Returns the first column with a valid call in a given row.
	 * @param row
	 * @return Last column.
	 */
	public int getFirstColumn(int row)
	{
		HaplotypeFragment pos = matrix.get(row);
		int firstColumn = pos.getFirstColumn();
		return firstColumn;
	}
	
	/**
	 * Returns the last column with a valid call in a given row.
	 * @param row
	 * @return Last column.
	 */
	public int getLastColumn(int row)
	{
		HaplotypeFragment pos = matrix.get(row);
		int lastColumn = pos.getLastColumn();
		return lastColumn;
	}
	/**
	 * Returns the number of fragments in the block.
	 * @return Number of fragments.
	 */
	public int getNumFragments()
	{
		int numberFragments = matrix.size();
		return numberFragments;
	}
	
	/**
	* Calculates the score of two fragments according to their hamming distance.
	* If the call is the same in both fragments it adds -1, if it is different it adds +1, if either is ALLELE_UNDECIDED it adds nothing.
	* <b> pre: </b> The matrix of fragments has been initialized.
 	* @param row1. 
 	* @param row2.
 	* @return hamming score.
 	*/
	public int getHamming2(int row1 , int row2)
	{
		int score = 0;
		for(int i = getFirstColumn(row2) ; i < getLastColumn(row1); i++)
		{
			byte allele1 = getAllele(row1, i);
			byte allele2 = getAllele(row2, i);
			if( allele1 != allele2 && (allele1 != CalledGenomicVariant.ALLELE_UNDECIDED && allele2!= CalledGenomicVariant.ALLELE_UNDECIDED))
			{
				score ++;
			} else if (allele1 == allele2 && (allele1 != CalledGenomicVariant.ALLELE_UNDECIDED && allele2!= CalledGenomicVariant.ALLELE_UNDECIDED))
			{
				score --;
			}
		}
			
		return score;
	}
	
	/**
	 * Returns the number of variants
	 * @return number of variants.
	 */
	public int getNumVariants()
	{
		return variants.size();
	}
	
	/**
	 * Set haplotype. 
	 * @param haplotype
	 */
	public void setHaplotype(byte [] haplotype)
	{
		this.haplotype = haplotype;
	}

	
}
