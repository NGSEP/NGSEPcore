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
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

public class GeneticAlgorithm
{

	private boolean verbose;
	private boolean saving;
	private int block;
	private String OUTPUT_NAME_FIT;
	private String OUTPUT_NAME_INFO;
	private String OUTPUT_NAME_HAP;

	private int POP_SIZE;
	private double MUT_RATE;
	private double CROSS_RATE;
	private int GENERATIONS;
	private int CHILDREN_PER_GEN;
	private int elitism;
	private int numInd;
	private String selection;
	private ArrayList<ArrayList<Integer> > M = new ArrayList<ArrayList<Integer>>();
	private ArrayList<ArrayList<Integer> > M_weight = new ArrayList<ArrayList<Integer>>();
	private ArrayList<ArrayList<Integer> > listIndex = new ArrayList<ArrayList<Integer>>();

	private String settings;
	private HaplotypeBlock haplotypeBlock;


	public void evolve(ArrayList<Chromosome> pop, int num_opt, int elitism, String method, int numberInd, Chromosome best)
	{

	}


	private void initialize(ArrayList<HaplotypeFragment> pop, int num_opt, HaplotypeFragment best)
	{
		for(int i=0; i< pop.size();i++)
		{
		
		}
	}



	private void evolve()
	{

	}
	private void getPartitions(int m, int n, ArrayList<Gene> genes, ArrayList<ArrayList<Integer> > M, ArrayList<ArrayList<Integer> > C1, ArrayList<ArrayList<Integer> > C2)
	{
		int num0 = 0;
		int num1 = 0;

		for(int i=0; i < m; i++)
		{
			if(genes.get(i).getValue() == 0)
				num0++;
			else
				num1++;
		}

		//		C1 = vector<vector<int> >(num0, vector<int>(n));
		//		C2 = vector<vector<int> >(num1, vector<int>(n));
		C1 = new ArrayList<ArrayList<Integer>>(); 
		C2 = new ArrayList<ArrayList<Integer>>(); 

		int count1 = 0;
		int	count2 = 0;
		for(int i=0; i < m; i++)
		{
			if(genes.get(i).getValue() == 0)
			{
				for(int j=0; j < n; j++)
				{
					C1.get(count1).set(i, M.get(i).get(j));
				}
				count1++;
			}
			else
			{
				for(int j=0; j < n; j++)
				{
					C2.get(count2).set(j, M.get(i).get(j));
				}
				count2++;
			}
		}
	}
	private void readParameters()
	{
		GENERATIONS = 100;
		POP_SIZE    = 100;
		CROSS_RATE  = 0.9;
		MUT_RATE    = 0.05;
		elitism     = 1;
		selection   = "tournament";
		numInd      = (int)(POP_SIZE*0.1);
		
	}
	
	public GeneticAlgorithm(HaplotypeBlock block)
	{
		startGA(block);
	}
	public void startGA(HaplotypeBlock block)
	{
		haplotypeBlock=block;
		readParameters();
		
		ArrayList<Chromosome> pop = new ArrayList<Chromosome>();
		
	}
	
	
}
