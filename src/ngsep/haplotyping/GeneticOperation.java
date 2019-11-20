/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2019 Jorge Duitama
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


public class GeneticOperation 
{


	public int crossover(Chromosome parent1, Chromosome parent2, ArrayList<Gene> genes, double cross_pointN) 
	{
		int cross_point= (int)cross_pointN;

		ArrayList<Gene> genesP1 = parent1.getGenes();
		ArrayList<Gene> genesP2 = parent2.getGenes();

		int numberGenes = parent1.getGenes().size();
		int half = (int) ((numberGenes*1.0) / 2.0);

		if(cross_point != -1.0)
		{
			if(cross_point >= half)
			{
				for(int i=0; i < (cross_point-half); i++)
					genes.add(genesP2.get(i));

				for(int i=cross_point-half; i < cross_point; i++)
					genes.add(genesP1.get(i));

				for(int i=cross_point; i < numberGenes; i++)
					genes.add(genesP2.get(i));			
			}
			else
			{
				for(int i=0; i < cross_point; i++)
					genes.add(genesP1.get(i));

				for(int i=cross_point; i < (cross_point+half); i++)
					genes.add(genesP2.get(i));

				for(int i=cross_point+half; i < numberGenes; i++)
					genes.add(genesP1.get(i));
			}
			return cross_point;
		}	
		else
		{	//uniform_int_distribution<int> distributionInt(i+1, genes.size()-1);
			//int endFlip = (int) (i+1 + Math.random()*(genes.size()-1));
			//uniform_int_distribution<int> distributionInt(0, numberGenes-1);
			int randNum = (int)(Math.random()*(numberGenes-1));

			if(randNum >= half)
			{
				for(int i=0; i < (randNum-half); i++)
					genes.add(genesP1.get(i));

				for(int i=randNum-half; i < randNum; i++)
					genes.add(genesP2.get(i));

				for(int i=randNum; i < numberGenes; i++)
					genes.add(genesP1.get(i));
			}
			else
			{			
				for(int i=0; i < randNum; i++)
					genes.add(genesP2.get(i));

				for(int i=randNum; i < (randNum+half); i++)
					genes.add(genesP1.get(i));

				for(int i=randNum+half; i < numberGenes; i++)
					genes.add(genesP2.get(i));
			}
			return randNum;
		}
	}

	public void mutate(ArrayList<Gene> genes, double rate, boolean flip) 
	{
		double random = Math.random();
		for(int i=0; i < genes.size(); i++)
		{
			if(random < rate)
			{
				if(flip)
				{
					if(i < (genes.size()-1))
					{
						//uniform_int_distribution<int> distributionInt(i+1, genes.size()-1);
						int endFlip = (int) (i+1 + Math.random()*(genes.size()-1));
						for(int j=i; j < endFlip; j++)
						{
							genes.get(j).mutatePosition();
						}
						break;
					}
					else
					{
						genes.get(i).mutatePosition();
					}
				}
				else
				{
					genes.get(i).mutatePosition();
				}
			}
		}
	}

}
