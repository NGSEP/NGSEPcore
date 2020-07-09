/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2019 Jorge Duitama, Christian Chavarro, Murray D. Patterson
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
import java.util.Arrays;

import ngsep.variants.CalledGenomicVariant;

public class HapChatSIHAlgorithm  implements SIHAlgorithm
{
	//The real PACBIO error rate is less than the constant, but this implementation 
	//has the same assumptions that the original paper.
	public static final Double ERROR_RATE_PACBIO=0.15;
	public static final Integer K_NUMBER_OF_CORRECTIONS=14; 
	private HaplotypeBlock block;
	private byte [] haplotype;
	private boolean [] cut;

	/**
	 * Build the haplotypes using the interface implementation of SIHAlgorithm
	 * This implementation is based in the HapChat Paper publication
	 */
	@Override
	public void buildHaplotype(HaplotypeBlock block) 
	{
		this.block=block;
		cut = new boolean [block.getNumFragments()];
		initCut(block);
		haplotype=CutHaplotypeTranslator.getHaplotype(block, cut, CutHaplotypeTranslator.CONSENSUS_COMBINED);
		preProcessHB();
		columCorrections();
		block.setHaplotype(haplotype);
	}

	/**
	 * Calculate the columns corrections of the calls in the HB
	 */
	private void columCorrections() 
	{
		int j=2;
		while(j<block.getCallsLenght())
		{
			//System.err.println("Entra a corregir columnas" + j);
		
				computeDistanceMin(j-1,j);
		
			j++;
		}
	}
	/**
	 * This method computes the distance for kMEC between the columns i and j
	 * @param i
	 * @param j
	 */
	private void computeDistanceMin(int i, int j) 
	{
		ArrayList<ArrayList<Byte>> kCorrectionsJ=new ArrayList<ArrayList<Byte>>();
		ArrayList<ArrayList<Byte>> kCorrectionsI=new ArrayList<ArrayList<Byte>>();
		ArrayList<Integer> distancesI= new ArrayList<Integer>(); 
		ArrayList<Integer> distancesJ= new ArrayList<Integer>(); 
		if(i>0&&j>0) 
		{
			ArrayList<Byte> iColumn = block.getColumn(i);
			kCorrectionsI=computekCorrections(iColumn);
			for (int z=0;z<kCorrectionsI.size();z++)
			{
				distancesI.add(getHammingScore(iColumn, kCorrectionsI.get(z)));
			}

			ArrayList<Byte> jColumn = block.getColumn(j);
			kCorrectionsJ=computekCorrections(jColumn);
			for (int z=0;z<kCorrectionsI.size();z++)
			{
				distancesJ.add(getHammingScore(jColumn, kCorrectionsJ.get(z)));
			}
			
			int minDistance=Integer.MAX_VALUE;
			ArrayList<Byte> jColumnCorrection =new ArrayList<Byte>();
			ArrayList<Byte> iColumnCorrection =new ArrayList<Byte>();
			for (int k = 0; k < kCorrectionsI.size(); k++) 
			{
				for (int l = 0; l < kCorrectionsI.size(); l++) 
				{
					int distance = getHammingScore(kCorrectionsI.get(k),kCorrectionsI.get(l));
					if(distance<minDistance)
					{
						minDistance=distance;
						jColumnCorrection =kCorrectionsI.get(l);
						iColumnCorrection =kCorrectionsI.get(k);
					}
				}
			}
			calculateHaplotype(jColumnCorrection, j);
			calculateHaplotype(iColumnCorrection, i);
		}


	}

	/**
	 * Computes all the k corrections of the column	
	 * @param iColumn list of the calls values in the different reads of the HB, this values represent the same column
	 * @return the K-Corrections of the column iColumn
	 */
	private ArrayList<ArrayList<Byte>> computekCorrections(ArrayList<Byte> iColumn) 
	{
		ArrayList<ArrayList<Byte>> kCorrections=new ArrayList<ArrayList<Byte>>();
		int i=iColumn.size();
		for(int j=0;j<i;j++)
		{
			kCorrections.add(kCorrection(iColumn));
		}
		return kCorrections;
	}

	/**
	 * This method recreates the preprocesing step of the HapChat Algorithm
	 */
	private void preProcessHB() 
	{
		for(int i=0;i<block.getNumFragments()-1;i++)
		{
			for (int j = 0; j < block.getNumFragments()-1; j++)
			{
				//First the method verifies if two fragments overlaps in the HB
				if(block.overlap(i, j))
				{
					calculateLikehood(i,j);
				}

			}
		}
	}

	/**
	 * This method calculated the likehood probability of the i and j  fragments. 
	 * In order two know if these tow fragments belongs to the same haplotype
	 * @param i 
	 * @param j
	 */
	private void calculateLikehood(int i, int j) 
	{
		int m=0;
		int x=0;
		HaplotypeFragment first = block.getHaplotypeFragment(i);
		HaplotypeFragment second = block.getHaplotypeFragment(j);
		int k;
		if(first.getCalls().length>second.getCalls().length)
		{
			k=first.getCalls().length;
		}
		else
		{
			k=second.getCalls().length;
		}
		for (int h=0; h < k; h++) 
		{
			if(block.getAllele(i, h)==block.getAllele(j, h))
			{
				m++;
			}
			else
			{
				x++;
			}
		}

		double ps=Math.pow(1-ERROR_RATE_PACBIO, 2*m)*Math.pow(ERROR_RATE_PACBIO, x)*Math.pow((1-ERROR_RATE_PACBIO/3), x);
		double pd =Math.pow(ERROR_RATE_PACBIO, m)*Math.pow((1-ERROR_RATE_PACBIO/3), m)*Math.pow(1-ERROR_RATE_PACBIO, x);
		double proablitiy= ps/pd;
		//Merges the two HF and delete the shortest for the HB
		if(proablitiy>1000000) 
		{
			if(first.getCalls().length>second.getCalls().length)
			{
				block.deleteFragment(i);
			}
			else
			{
				block.deleteFragment(j);
			}

		}

	}

	/**
	 * Calculates the hamming distance of the columns Mj and Cj
	 * @param Mj
	 * @param Cj
	 * @return the hamming distance of the columns
	 */
	private int getHammingScore( ArrayList<Byte> Mj, ArrayList<Byte> Cj)
	{
		int hamming=0;
		int tamanioMJ= Mj.size();
		int tamanioCJ= Cj.size();
		if(tamanioCJ==tamanioMJ)
		{
			for (int i=0; i<tamanioCJ;i++)
			{
				if(Mj.get(i)!=Cj.get(i))
				{
					hamming++;
				}
			}
		}
		return hamming;

	}

	/**
	 * Calculates a modification of the values of the calls in the colums M.
	 * The details of this method are defined in the original HapChat Paper
	 * @param m column to modified
	 * @return a modification of the values in the original m array
	 */
	private ArrayList<Byte> kCorrection(ArrayList<Byte> m)
	{
		ArrayList<Byte> mc=m;
		int numberCorrectionsActual=0;
		for(int i=0; i<m.size()&&numberCorrectionsActual<K_NUMBER_OF_CORRECTIONS;i++)
		{
			int correct= (int) ((Math.random()*100)%4);
			byte actual=m.get(i);
			if(correct==0&& actual!=CalledGenomicVariant.ALLELE_UNDECIDED)
			{
				if(actual==CalledGenomicVariant.ALLELE_REFERENCE)
				{
					mc.set(i, CalledGenomicVariant.ALLELE_ALTERNATIVE);
				}
				else
				{
					mc.set(i, CalledGenomicVariant.ALLELE_REFERENCE);
				}
				numberCorrectionsActual++;
			}
		}
		return mc;
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

	/**
	 * Calculate the haplotype based on a haplotype calculate on a given column
	 * @param haplotype calculated until now
	 * @param j column to calculate the haplotype
	 */
	private void calculateHaplotype(ArrayList<Byte> haplo, int j) 
	{

		int count0=0;
		int count1=0;

		for(int i= 0; i<haplo.size();i++)
		{
			if(haplo.get(i)==CalledGenomicVariant.ALLELE_UNDECIDED)
			{

			}
			else if(haplo.get(i)==CalledGenomicVariant.ALLELE_REFERENCE)
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
			haplotype[j]= CalledGenomicVariant.ALLELE_REFERENCE;
		}
		else
		{
			haplotype[j] =CalledGenomicVariant.ALLELE_ALTERNATIVE;
		}
	}



}



