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

public class Chromosome 
{
	private double fitness;
	private int cross_point;

	private ArrayList<Integer> indexC1 = new ArrayList<Integer>();
	private ArrayList<Integer> indexC2 = new ArrayList<Integer>();

	private ArrayList<Integer> h1 = new ArrayList<Integer>();
	private ArrayList<Integer> h2 = new ArrayList<Integer>();

	public ArrayList<Gene> genes;

	public Chromosome()
	{

	}

	public Chromosome(ArrayList<ArrayList<Integer>> M, ArrayList<ArrayList<Integer>> M_weight, ArrayList<ArrayList<Integer>> listIndex)
	{
		for(int i = 0; i < M.size(); i++)
		{
			Gene gene = new Gene();
			//The source code in c++ uses the method push_back from the Vector structure, this in Arraylist
			// is represented by the operation .add()
			genes.add(gene);
		}
		calculateFitness(M, M_weight, listIndex);
	}

	public Chromosome (ArrayList<ArrayList<Integer>> M, ArrayList<ArrayList<Integer>> M_weight, ArrayList<ArrayList<Integer>> listIndex,
			double mut_rate, Chromosome parent1, Chromosome parent2, double cross_pointIn, boolean flip)
	{
		GeneticOperation op = new GeneticOperation();
		cross_point = op.crossover(parent1, parent2, genes, cross_pointIn);
		op.mutate(genes, mut_rate, flip);
		calculateFitness(M, M_weight, listIndex);
	}



	public Chromosome(ArrayList<Chromosome> pop) {
		// TODO Auto-generated constructor stub
	}

	public void calculateFitness(ArrayList<ArrayList<Integer>> M, ArrayList<ArrayList<Integer>> M_weight,
			ArrayList<ArrayList<Integer>> listIndex) 
	{
		int m = M.size();
		int n = M.get(0).size();

		int num0 = 0;
		int num1 = 0;

		for(int i=0; i < m; i++)
		{
			if(genes.get(i).getValue() == 0)
				num0++;
			else
				num1++;
		}

		ArrayList<ArrayList<Integer>> C1 = new ArrayList<ArrayList<Integer>>(); 
		ArrayList<ArrayList<Integer>> C2 = new ArrayList<ArrayList<Integer>>(); 
		int count1 = 0;
		int	count2 = 0;
		for(int i=0; i < m; i++)
		{
			if(genes.get(i).getValue() == 0)
			{
				for(int j=0; j < n; j++)
				{
					//C1[count1][j] = M[i][j];
					C1.get(count1).set(j, M.get(i).get(j));
				}
				count1++;
				indexC1.add(i);
			}
			else
			{
				for(int j=0; j < n; j++)
				{
					//	C2[count2][j] = M[i][j];
					C2.get(count2).set(j, M.get(i).get(j));
				}
				count2++;
				indexC2.add(i);
			}
		}
		fitness = errorFunction(n, num0, num1, C1, C2, indexC1, indexC2, M_weight, listIndex);

	}

	private ArrayList<Integer> calculate_h(ArrayList<Integer> N0, ArrayList<Integer> N1)
	{
		int n = N0.size();
		ArrayList<Integer> h = new ArrayList<Integer>();

		for(int j=0; j < n; j++)
		{

			if( (N0.get(j) == 0) & (N1.get(j) == 0) )
			{
				h.set(j, -1);
			}
			else
			{
				if(N0.get(j) > N1.get(j))
				{
					h.set(j, 0);
				}
				else
				{
					h.set(j, 1);
				}
			}
		}

		return h;
	}

	private int distanceFragment(ArrayList<Integer> f1, ArrayList<Integer> f2)
	{
		int sumErr = 0;

		for(int i=0; i < f1.size(); i++)
		{
			sumErr += singleDistance(f1.get(i),f2.get(i));
		}
		return sumErr;
	}
	private int singleDistance(int x, int y)
	{
		{
			if( (x == -1) & (y == -1) )
			{
				return 0;
			}
			else
			{
				if( (x != -1) & (y != -1) & (x != y) )
				{
					return 1;
				}
				else
				{
					return 0;
				}
			}
		}

	}

	private void calculate_n0_n1(ArrayList<ArrayList<Integer>> C, ArrayList<ArrayList<Integer>> M_weight,
			ArrayList<ArrayList<Integer>> listIndex, ArrayList<Integer> indexC, ArrayList<Integer> N0, ArrayList<Integer> N1)
	{
		int m = C.size();
		int idx = 0;

		for(int i=0; i < m; i++)
		{
			for(int j=0; j < listIndex.get(indexC.get(i)).get(1); j++)
				// for(int j=listIndex[indexC[i]][0]; j < listIndex[indexC[i]][1]; j++)
			{
				if(C.get(i).get(j)== 0)
				{
					idx = indexC.get(i);
					N0.set(j, N0.get(j)+M_weight.get(idx).get(j));
				}
				else
				{
					idx = indexC.get(i);
					//N1[j] += M_weight[idx][j];
					N1.set(j, N0.get(j)+M_weight.get(idx).get(j));
					// if(C[i][j] == 1)
					// {
					// 	idx = indexC[i];
					// 	N1[j] += M_weight[idx][j];
					// }	
				}
			}
		}
	}

	private	double errorFunction(int n, int num0, int num1, ArrayList<ArrayList<Integer> > C1, ArrayList<ArrayList<Integer> > C2, 
			ArrayList<Integer> indexC1, ArrayList<Integer> indexC2, ArrayList<ArrayList<Integer> > M_weight, ArrayList<ArrayList<Integer> > listIndex)
	{
		double error1 = 0;
		double error2 = 0;


		ArrayList<Integer> h1 = new ArrayList<Integer>();
		ArrayList<Integer> h2 = new ArrayList<Integer>();
		ArrayList<Integer> N0_1 = new ArrayList<Integer>();
		ArrayList<Integer> N1_1 = new ArrayList<Integer>();
		ArrayList<Integer> N0_2 = new ArrayList<Integer>();
		ArrayList<Integer> N1_2 = new ArrayList<Integer>();

		if(num0==1)
		{
			calculate_n0_n1(C1, M_weight, listIndex, indexC1, N0_1, N1_1);
			h1 = calculate_h(N0_1, N1_1);
		}

		if(num1==1)
		{
			calculate_n0_n1(C2, M_weight, listIndex, indexC2, N0_2, N1_2);
			h2 = calculate_h(N0_2, N1_2);
		}

		if( (num0 > 0) & (num1 > 0) )
		{
			for(int j=0; j < h1.size(); j++)
			{
				if( (h1.get(j) == h2.get(j)) & (h1.get(j) == -1))
				{
					continue;
				}
				else
				{
					if(h1.get(j) == h2.get(j))
					{
						if(h1.get(j) == 0)
						{
							if(N0_1.get(j) > N0_2.get(j))
								h2.set(j, 1);
							else
								h1.set(j, 1);
						}
						else
						{
							if(N1_1.get(j) > N1_2.get(j))
								h2.set(j, 0);
							else
								h1.set(j, 0);
						}
					}
					else
					{
						// no ambigue positions
						if(h1.get(j) == -1 & h2.get(j)==0)
							h1.set(j, 1);
						else if(h1.get(j) == -1 & h2.get(j)==1)
							h1.set(j, 0);
						else if(h2.get(j) == -1 & h1.get(j)==0)
							h2.set(j, 1);
						else if(h2.get(j) == -1 & h1.get(j)==1)
							h2.set(j, 0);

					}
				}
			}
		}

		if(num0==1)
		{
			for(int j=0; j < C1.size(); j++)
			{
				error1 += distanceFragment(C1.get(j), h1);
			}
		}
		if(num1==1)
		{
			for(int j=0; j < C2.size(); j++)
			{
				error2 += distanceFragment(C2.get(j), h2);
			}
		}

		return error1+error2;
	}

	public double getFitness()
	{ 
		return fitness;
	}

	public int getCrossPoint()
	{ 
		return cross_point;
	}

	public ArrayList<Integer> getIndexC1() 
	{
		return indexC1;
	}

	public ArrayList<Integer> getIndexC2() 
	{
		return indexC2;
	}

	public ArrayList<Integer> getH1() 
	{
		return h1;
	}

	public ArrayList<Integer> getH2()
	{
		return h2;
	}

	public ArrayList<Gene> getGenes() 
	{
		return genes;
	}

}
