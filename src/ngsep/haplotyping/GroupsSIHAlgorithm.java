/*******************************************************************************
 * SingleIndividualHaplotyper - Efficient heuristic algorithms for the SIH problem
 * Copyright 2011 Jorge Duitama, Christian Chavarro, Daniel Bautista
 *
 * This file is part of SingleIndividualHaplotyper.
 *
 *     SingleIndividualHaplotyper is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     SingleIndividualHaplotyper is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with SingleIndividualHaplotyper.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
package ngsep.haplotyping;

import java.util.ArrayList;
import java.util.Arrays;

import ngsep.variants.CalledGenomicVariant;

public class GroupsSIHAlgorithm implements SIHAlgorithm 
{
	private boolean [] cut;
	private byte [] haplotype;

	/**
	 * Build the haplotypes using the interface implementation of SIHAlgorithm
	 */
	@Override
	public void buildHaplotype(HaplotypeBlock block)
	{
		int numFragments = block.getNumFragments();

		Long init= System.currentTimeMillis();
		int overlapCountMax=0;
		cut = new boolean [block.getNumFragments()];
		initCut(block);
		haplotype=CutHaplotypeTranslator.getHaplotype(block, cut, CutHaplotypeTranslator.CONSENSUS_COMBINED);
		//HashMap<Double, ArrayList<HaplotypeFragment>> fragmentsArray;

		ArrayList<Integer> fragmentsArray= new ArrayList<Integer>();
		int fragmentMaxId=0;
		ArrayList<Integer> fragmentsArrayActual = new ArrayList<Integer>();
		int overlapActual=0;
		for(int i=0;i<numFragments;i++) 
		{	//HashMap<Double, ArrayList<HaplotypeFragment>> fragmentsArrayActual= new HashMap<Double, ArrayList<HaplotypeFragment>>();

			for(int j=i+1;j<numFragments;j++) 
			{

				if(!block.overlap(i, j)) 
				{
					break;
				}
				overlapActual++;
				double score = getScore(block, i,j,true);
				if(score<=500.0)
				{

					fragmentsArrayActual.add(j);

				}



			}
			if(overlapCountMax<overlapActual)
			{
				fragmentMaxId=i;
				overlapCountMax=overlapActual;
				fragmentsArray= fragmentsArrayActual;
			}


		}
		calculateHaplotype(fragmentMaxId, fragmentsArray, block);
		block.setHaplotype(haplotype);
		Long fin= System.currentTimeMillis();
		System.out.println("tiempo ejecucion " + (fin-init));
	}

	/**
	 * Calculates the haplotype of the haplotype block
	 * @param fragmentMaxId fragment of the HB with the maximun overlap
	 * @param fragmentsArray fragments that overlaps with the maximun overlap fragment
	 * @param block HB that contains the fragments calls
	 */
	private void calculateHaplotype(int fragmentMaxId, ArrayList<Integer> fragmentsArray, HaplotypeBlock block) 
	{
		// TODO Auto-generated method stub
		HaplotypeFragment f= block.getHaplotypeFragment(fragmentMaxId);
		int columns = f.getCalls().length;



		for (int j=0;j<columns;j++)
		{
			byte hap =block.getAllele(fragmentMaxId, j);
			int count0=0;
			int count1=0;

			for(int i= 0; i<fragmentsArray.size();i++)
			{
				if(block.getAllele(j, i)==CalledGenomicVariant.ALLELE_UNDECIDED)
				{

				}
				else if(block.getAllele(j, i)==CalledGenomicVariant.ALLELE_REFERENCE)
				{
					count0++;
				}
				else
				{
					count1++;
				}


			}
			if(count0>count1)
			{
				hap= CalledGenomicVariant.ALLELE_REFERENCE;
			}
			else
			{
				hap =CalledGenomicVariant.ALLELE_ALTERNATIVE;
			}
			if(haplotype.length>j)
			{
				if(haplotype[j]!=hap) 
				{
					updateHaplotype(haplotype, block, j, true);
				}
			}

		}



	}

	/**
	 * Inits a random haplotype cut in the haplotype block given as parameter
	 * @param b HB to init the cut
	 */
	private void initCut(HaplotypeBlock b) 
	{
		byte [] hap = new byte [b.getNumVariants()];

		for(int i=0;i<hap.length;i++)
		{
			hap[i] = CalledGenomicVariant.ALLELE_UNDECIDED;
		}
		int n = b.getNumFragments();
		boolean [] assigned = new boolean [n];
		Arrays.fill(assigned, false);
		//Get the fragment with larger number of calls
		int maxC=0;
		int maxI=0;
		for(int i=0;i<n;i++) 
		{
			int calls = b.getFragmentCalls(i);
			if(calls>maxC) 
			{
				maxC=calls;
				maxI=i;
			}
		}
		assigned[maxI] = true;
		cut[maxI] = false;
		updateHaplotype(hap, b, maxI, false);
		//Assign the other fragments
		for(int h=0;h<n-1;h++) 
		{
			int maxAbsScore =0;
			int maxScore =0;
			//Find the best next fragment
			for(int i=0;i<n;i++) 
			{
				if(!assigned[i]) 
				{
					int score = b.getHamming2(hap,i);
					int absScore = Math.abs(score);
					if(absScore > maxAbsScore)
					{
						maxAbsScore = absScore;
						maxI = i;
						maxScore = score;
					}
				}
			}
			assigned[maxI] = true;
			cut[maxI] = (maxScore==maxAbsScore);
			updateHaplotype(hap, b, maxI, cut[maxI]);
		}
	}


	/**
	 * Calculate the hamming score for two rows
	 * @param block
	 * @param row1
	 * @param row2
	 * @param useQualityScores
	 * @return
	 */
	public double getScore(HaplotypeBlock block, int row1, int row2, boolean useQualityScores)
	{

		int score = block.getHamming2(row1, row2);
		return score;
	}

	/**
	 * Updates the haplotype in a row
	 * @param hap haplotype to modify
	 * @param b  HB to modify the haplotype
	 * @param row to modify
	 * @param reverse true if the reverse is the new haplotype in the row in order to modify
	 */
	private void updateHaplotype(byte [] hap, HaplotypeBlock b, int row, boolean reverse)
	{
		int firstJ = b.getFirstColumn(row);
		int lastJ = b.getLastColumn(row);
		for(int j=firstJ;j<=lastJ;j++) 
		{
			byte allele = b.getAllele(row, j);
			if(hap[j]==CalledGenomicVariant.ALLELE_UNDECIDED && allele!=CalledGenomicVariant.ALLELE_UNDECIDED) 
			{
				//TODO: Fix for multiallelic
				if(!reverse) 
				{
					hap[j] = allele;
				} else if (allele == CalledGenomicVariant.ALLELE_REFERENCE) 
				{
					hap[j] = CalledGenomicVariant.ALLELE_ALTERNATIVE;
				} else
				{
					hap[j] = CalledGenomicVariant.ALLELE_REFERENCE;
				}

			}
		}

	}


	//	/**
	//	 * Updates the haplotypes cut in a position
	//	 * @param b
	//	 * @param haplotype
	//	 */
	//	private void updateCut(HaplotypeBlock b, byte [] haplotype) 
	//	{
	//		int n = b.getNumFragments();
	//		for(int i=0;i<n;i++) {
	//			int score = b.getHamming2(haplotype,i);
	//			if(score != 0) {
	//				cut[i] = score > 0;
	//			}
	//			
	//		}
	//	}

}