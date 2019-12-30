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
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import ngsep.variants.CalledGenomicVariant;
import ngsep.variants.CalledSNV;
import ngsep.variants.GenomicVariant;

public class HaplotypeBlock 
{

	/**
	 * Represents the matrix of fragments and variants.
	 */
	private List<HaplotypeFragment> matrix;

	/**
	 * Represents the list of variants.
	 */
	private List<CalledGenomicVariant> calls;
	
	/**
	 * Represents a haplotype.
	 */
	private byte haplotype[];
	
	/**
	 * Indicates if the matrix is already sorted
	 */
	private boolean sorted = true;

	
	/**
	 * Constructor that initializes the attributes of a HaplotypeBlock with the given parameters.
	 * @param variants.
	 */
	public HaplotypeBlock(List<CalledGenomicVariant> calls) 
	{
		this.calls = calls;
		matrix = new ArrayList <HaplotypeFragment>();
		haplotype = null;

	}
	
	/**
	 * Add fragment to the matrix
	 * @param firstColumn where valid allele calls are found
	 * @param alleleCalls Calls starting from the given column
	 */
	public void addFragment(int firstColumn, byte[] alleleCalls) {
		HaplotypeFragment fragment = new HaplotypeFragment(firstColumn, alleleCalls);
		matrix.add(fragment);
		sorted = false;
	}
		
	/**
	 * Returns the call in a given position in the matrix of fragments.
	 * @param i row of the matrix
	 * @param j column of the matrix
	 * @return byte Allele call at position i,j
	 */
	public byte getAllele(int i, int j)
	{
		sort();
		HaplotypeFragment row = matrix.get(i);
		byte allele = row.getCall(j);
		return allele;
	}
	
	/**
	 * Returns the haplotype phasing the variants.
	 * @return byte [] haplotype configuration.
	 */
	public byte [] getHaplotype()
	{
		return haplotype;
	}
	
	/**
	 * Returns the variant in the given position in the list of variants.
	 * <b> pre: </b> The list of variants has been initialized.
	 * @param column of the matrix
	 * @return GenomicVariant associated with the column
	 */
	public GenomicVariant getVariant(int column)
	{
		return calls.get(column);
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
		sort();
		int score = 0;
		int lastColRow1 = getLastColumn(row1);
		for(int i = getFirstColumn(row2) ; i <=lastColRow1 ; i++) {
			byte allele1 = getAllele(row1, i);
			byte allele2 = getAllele(row2, i);
			score+=getHammingScore(allele1, allele2, false);
		}	
		return score;
	}
	
	/**
	 * Calculates the score of two fragments according to their hamming distance.
	 * If the call is the same in both fragments it adds -1, if it is different it adds +1, if either is ALLELE_UNDECIDED it adds nothing.
	 * <b> pre: </b> The matrix of fragments has been initialized.
 	 * @param row1. 
 	 * @param row2.
 	 * @return hamming score.
 	 */
	public int getHamming2(int row1, int row2)
	{
		sort();
		int score = 0;
		int lastColRow1 = getLastColumn(row1);
		for(int i = getFirstColumn(row2) ; i <=lastColRow1 ; i++) 
		{
			byte allele1 = getAllele(row1, i);
			byte allele2 = getAllele(row2, i);
			score+=getHammingScore(allele1, allele2, true);
		}
		return score;
	}
	
	/**
	 * Calculates the hamming2 score of a haplotype against a fragment
	 * If the call is the same in both fragments it adds -1, if it is different it adds +1, if either is ALLELE_UNDECIDED it adds nothing.
	 * <b> pre: </b> The matrix of fragments has been initialized.
 	 * @param haplotype with length equal to the number of variants 
 	 * @param row of the matrix to calculate the score
 	 * @return int Modified hamming distance score as defined above
 	 */
	public int getHamming2(byte [] haplotype, int row)
	{
		sort();
		int score = 0;
		int lastColRow = getLastColumn(row);
		for(int j = getFirstColumn(row) ; j <=lastColRow ; j++) {
			byte allele1 = haplotype[j];
			byte allele2 = getAllele(row, j);
			score+=getHammingScore(allele1, allele2, true);
		}
		return score;
	}
	
	/**
	 * Calculates the hamming score of two alleles
	 * @param allele1
	 * @param allele2
	 * @param type2
	 * @return
	 */
	private int getHammingScore (byte allele1, byte allele2, boolean type2) 
	{
		if(allele1 != CalledGenomicVariant.ALLELE_UNDECIDED && allele2!= CalledGenomicVariant.ALLELE_UNDECIDED) {
			if( allele1 != allele2)
			{
				return 1;
			} else if(type2)
			{
				return -1;
			}
		}
		return 0;
	}
	
	/**
	 * Checks if two fragments overlap
	 * <b> pre: </b> The matrix of fragments has been initialized.
	 * @param row1 First row to compare
	 * @param row2 Second row to compare
	 * @return True when the two fragments overlap.
	 */
	public boolean overlap(int row1, int row2)
	{
		sort();
		return getFirstColumn(row1) <= getLastColumn(row2) && getFirstColumn(row2) <= getLastColumn(row1);
	}
	
	/**
	 * Returns the first column with a valid call in a given row.
	 * @param row
	 * @return Last column.
	 */
	public int getFirstColumn(int row)
	{
		sort();
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
		sort();
		HaplotypeFragment fragment = matrix.get(row);
		return fragment.getLastColumn();
	}
	
	/**
	 * Returns the number of fragments in the block.
	 * @return Number of fragments.
	 */
	public int getNumFragments()
	{
		return matrix.size();
	}
	
	/**
	 * Return the HF in a n position of the matrix
	 * @param n fragment position to get
	 * @return Haplotype fragment.
	 */
	public HaplotypeFragment getHaplotypeFragment(int n)
	{
		return matrix.get(n);
	}

	/**
	 * Obtains the calls that are in a column of the Haplotype Block 
	 * @param j
	 * @return
	 */
	public ArrayList<Byte> getColumn(int j)
	{
		ArrayList<Byte> column = new ArrayList<Byte>();
		for(int i =0; i<matrix.size();i++)
		{
		column.add(getAllele(i, j));
		}
		
		return column;
		
	}
	
	/**
	 * Returns the number of variants
	 * @return number of variants.
	 */
	public int getNumVariants()
	{
		return calls.size();
	}
	
	/**
	 * Return the number of non-undecided calls within a specific fragment
	 * @param row where the fragment is located
	 * @return int Number of non undecided calls
	 */
	public int getFragmentCalls(int row) 
	{
		int firstJ = getFirstColumn(row);
		int lastJ = getLastColumn(row);
		int count = 0;
		for(int j=firstJ;j<=lastJ;j++)
		{
			if(getAllele(row, j)!=CalledGenomicVariant.ALLELE_UNDECIDED) 
			{
				count++;
			}
		}
		return count;
	}
	
	/**
	 * Changes the haplotype corresponding to the given block. 
	 * @param haplotype new haplotype
	 */
	public void setHaplotype(byte [] haplotype)
	{
		this.haplotype = haplotype;
	}
	
	/**
	 * Sorts the matrix by first position of the fragment
	 */
	private void sort() 
	{
		if(sorted) return;
		Collections.sort(matrix, new Comparator<HaplotypeFragment>() 
		{

			@Override
			public int compare(HaplotypeFragment f1, HaplotypeFragment f2) 
			{
				return f1.getFirstColumn()-f2.getFirstColumn();
			}
		});
		sorted = true;
		
	}
	
	/**
	 * Phase the calls within the block using the given haplotype
	 */
	public void phaseCallsWithHaplotype() 
	{
		for(int i=0;i<haplotype.length;i++)
		{
			CalledGenomicVariant call = calls.get(i);
			if(call instanceof CalledSNV) ((CalledSNV)call).setPhasingCN2(haplotype[i]==CalledGenomicVariant.ALLELE_ALTERNATIVE);
		}
	}
	
	/**
	 * Deletes the fragment j of the haplotype block
	 * @param j
	 */
	public void deleteFragment(int j) 
	{
		matrix.remove(j);
		sorted = false;
		
	}
	
	/**
	 * Return the number of calls in the Haplotype Block
	 * @return
	 */
	public int getCallsLenght()
	{
		return calls.size();
	}
}

