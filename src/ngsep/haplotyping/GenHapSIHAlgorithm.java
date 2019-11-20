/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2019 Jorge Duitama, Christian Chavarro, Daniel Bautista
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
import java.util.Collections;
import java.util.List;

import ngsep.variants.CalledGenomicVariant;

public class GenHapSIHAlgorithm implements SIHAlgorithm 
{	
	private boolean [] cut;
	private byte [] haplotype;
	private byte [] haplotype0;
	private byte [] haplotype1;
	private ArrayList<boolean[]> cutsPoblations = new ArrayList<boolean[]>();
	private ArrayList<Integer> fitnessArray= new ArrayList<Integer>();
	private HaplotypeBlock block;


	@Override
	public void buildHaplotype(HaplotypeBlock block) 
	{
		this.block=block;


		cut = new boolean [block.getNumFragments()];
		initCut(block);
		haplotype0=CutHaplotypeTranslator.getHaplotype(block, cut, CutHaplotypeTranslator.CONSENSUS_COMBINED);

		initCut(block);
		haplotype1=CutHaplotypeTranslator.getHaplotype(block, cut, CutHaplotypeTranslator.CONSENSUS_COMBINED);

		//Inicia 100 cortes para analizar las posibilidades, donde 100 se define en el paper como 100 luego de las pruebas

		for (int i=0; i<100;i++)
		{ 
			boolean [] cutA = new boolean [block.getNumFragments()];
			initCut(block, cutA);
			cutsPoblations.add(cutA);
		}

		calculateCuts(block);




		calculateHaplotype();
		haplotype=CutHaplotypeTranslator.getHaplotype(block, cut, CutHaplotypeTranslator.CONSENSUS_COMBINED);
		block.setHaplotype(haplotype);
	}

	private void calculateCuts(HaplotypeBlock block) 
	{
		int fitnessMaxAnterior=0;
		int fitnessMax=0;
		if(fitnessArray.isEmpty()==false) 
		{
			fitnessMax=calculateFitnessMax();	
		}
		int countStop=0;
		if(fitnessMaxAnterior==fitnessMax)
		{
			countStop=0;
		}
		else
		{
			countStop++;
		}
	
		for (int i=0; i<100&&countStop!=25;i++)
		{ 
	
			for (int j = 0; j<block.getNumFragments();j++)
			{
				HaplotypeFragment actualFragment= block.getHaplotypeFragment(j);
				if(actualFragment.getCall(j)==CalledGenomicVariant.ALLELE_REFERENCE)
				{
					cutsPoblations.get(i).clone()[j]=true;
				}
				else if (actualFragment.getCall(j)==CalledGenomicVariant.ALLELE_ALTERNATIVE)
				{
					cutsPoblations.get(i).clone()[j]=false;
				}
			}
			calculateHaplotype(haplotype0);
			calculateHaplotype(haplotype1);
			calculateFitness();
			recalculateCuts();
		}
	
	
	}

	



	private void recalculateCuts()
	{
		ArrayList<boolean[]> newCutsPoblations = cutsPoblations;
		int pDeletes=0;
		int maxFitness = maxIndex(fitnessArray);
		for(;(int)pDeletes<newCutsPoblations.size()*0.1;pDeletes++)
		{
			int toDelete= (int) (Math.random()*100);

			if(toDelete!=maxFitness && toDelete<newCutsPoblations.size())
			{
				newCutsPoblations.remove(toDelete);
			}

		}

		if(newCutsPoblations.size()<=(int)cutsPoblations.size()*0.9)
		{
			int toDelete2= (int) Math.abs(newCutsPoblations.size()-(int)cutsPoblations.size()*0.9);
			for(;toDelete2>0;toDelete2--)
			{
				int toDelete= (int) (Math.random()*100);

				if(toDelete!=maxFitness && toDelete<newCutsPoblations.size())
				{
					newCutsPoblations.remove(toDelete);
				}

			}
		}
		for(int i=0; i<10;i++)
		{
			newCutsPoblations.add(mutateOrCross());
		}
		cutsPoblations=newCutsPoblations;
	}

	private boolean[] mutateOrCross() 
	{


		if((Math.random()*20)%2==1)
		{
			return mutate();
		}else
		{
			return crossover();
		}
	}
	private boolean[] crossover() 
	{
		int toCross1=(int) ((int)cutsPoblations.size()*Math.random());
	
		boolean[] cutCrossed= cutsPoblations.get(toCross1);
		int toCross2=(int) ((int)cutsPoblations.size()*Math.random());
		boolean[] cutCross2= cutsPoblations.get(toCross2);
		int avgLengh = (int)(cutCrossed.length+cutCross2.length)/2;
		int toCroos= (int) (avgLengh+Math.random()*avgLengh);
		for (int i = toCroos; toCroos < cutCross2.length; i++) 
		{
			cutCrossed[i]=cutCross2[i];
		}
		return cutCrossed;
	}

	public boolean[] mutate() 
	{

		int toMutate=(int) ((int)cutsPoblations.size()/Math.random());
		boolean[] cutMutate= cutsPoblations.get(toMutate);
		for(int i=0;i<cutMutate.length;i++)
		{
			if(Math.random()<0.05)
			{
				cutMutate[i]=!cutMutate[i];
			}
		}

		return cutMutate;

	}




	public  int maxIndex(List<Integer> list)
	{
		Integer i=0, maxIndex=-1, max=null;
		for (Integer x : list) {
			if ((x!=null) && ((max==null) || (x>max))) 
			{
				max = x;
				maxIndex = i;
			}
			i++;
		}
		return maxIndex;
	}

	private void calculateHaplotype(byte [] haplotype)
	{

		for (int j=0;j<block.getNumVariants();j++)
		{
			int count0=0;
			int count1=0;

			for (int i=0;i<block.getNumFragments();i++)
			{


				if(block.getAllele(i, j)==CalledGenomicVariant.ALLELE_REFERENCE)
				{
					//no hace nada
					count0++;
				}
				else if(block.getAllele(i, j)==CalledGenomicVariant.ALLELE_ALTERNATIVE)
				{
					count1++;
				}

			}
			byte hap;
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
					haplotype[j]=hap;
				}
			}

		}
	}

	private void calculateFitness() 
	{ 
		for (int j = 0; j<cutsPoblations.size();j++)
		{
			boolean [] cutActual = cutsPoblations.get(j);
			int fitnessValue=0;
			for (int i=0;i<cutActual.length;i++)
			{
				if(cutActual[i])
				{
					fitnessValue+=block.getHamming2(haplotype0, i);
				}
				else
				{
					fitnessValue+=block.getHamming2(haplotype1, i);

				}
			}
			fitnessArray.add(fitnessValue);	
		}

	}

	private int calculateFitnessMax() 
	{
		return Collections.max(fitnessArray);
	}

	private void initCut(HaplotypeBlock b, boolean[] cutA)
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
		cutA[maxI] = false;
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
			cutA[maxI] = (maxScore==maxAbsScore);
			updateHaplotype(hap, b, maxI, cutA[maxI]);
		}

	}

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

	private void calculateHaplotype() 
	{
	haplotype=haplotype0;

		for (int j=0;j<haplotype0.length;j++)
		{
			if(haplotype0[j]==haplotype1[j]&& haplotype0[j]!=CalledGenomicVariant.ALLELE_UNDECIDED
					&& haplotype1[j]!=CalledGenomicVariant.ALLELE_UNDECIDED) 
			{
				setHaplotypeSol(haplotype0[j],j);
			}
			else if(haplotype0[j]!=CalledGenomicVariant.ALLELE_UNDECIDED)
			{
				setHaplotypeSol(haplotype0[j],j);
			}
			else
			{
				setHaplotypeSol(haplotype1[j],j);
			}


		}



	}

	private void setHaplotypeSol(byte b, int j) 
	{
		
		if(haplotype.length>j)
		{
			if(haplotype[j]!=b) 
			{
				updateHaplotype(haplotype, block, j, true);
			}
		}

	}

}
