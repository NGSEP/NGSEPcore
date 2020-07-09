/*******************************************************************************
 * SingleIndividualHaplotyper - Efficient heuristic algorithms for the SIH problem
 * Copyright 2011 Jorge Duitama
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
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class FragmentsCutBuilder {
	private List<Vertex> graph;
	private List<Edge> allEdges;
	private boolean [] cut;
	/**
	 * 
	 * @param block
	 */
	public FragmentsCutBuilder(HaplotypeBlock block) {
		this(block,false);
	}
	/**
	 * @param block
	 */
	public FragmentsCutBuilder(HaplotypeBlock block, boolean useQualityScores) {

		super();
		int numFragments = block.getNumFragments();
		graph = new ArrayList<Vertex>(numFragments);
		allEdges = new ArrayList<Edge>(numFragments);
		cut = new boolean [numFragments];
		for(int i=0;i<numFragments;i++) {
			graph.add(new Vertex(i));
		}
		for(int i=0;i<numFragments;i++) {
			for(int j=i+1;j<numFragments;j++) {
				if(!block.overlap(i, j)) {
					break;
				}
				double score = getScore(block, i,j,useQualityScores);
				if(score!=0) {
					Edge e = new Edge(i, j, score);
					graph.get(i).addEdge(e);
					graph.get(j).addEdge(e);
					allEdges.add(e);
				}
			}
		}
		Collections.sort(allEdges, new EdgesComparator());
	}
	public double [] calcAgreementScoreNodes() {
		double [] answer = new double [graph.size()];
		int i=0;
		for(Vertex v:graph) {
			for(Edge e: v.getEdges()) {
				answer[i]+=Math.abs(e.getWeight());
			}
			i++;
		}
		return answer;
	}
	public double getScore(HaplotypeBlock block, int row1, int row2, boolean useQualityScores) {
		//if(!useQualityScores) {
		int score = block.getHamming2(row1, row2);
		return score;
		/*int overlap = f1.getOverlappingCount(f2);
			int disagree = f1.getHammingDistance(f2);
			int agree = overlap - disagree;
			return disagree*disagree - agree*agree;
			double score2 = Math.pow(score, 2);
			if(score >= 0)
				return (int) Math.round(score2);
			else {
				return (int) -Math.round(score2);
			}*/
		//}
		//else return f1.getDistanceWithQuals(f2);
	}
	/**
	 * @return the cut
	 */
	public boolean[] getCut() {
		return cut;
	}

	public void calculateMaxCut() {
		//randomizeCut();
		boolean [] bestCut = new boolean [cut.length];
		double maxScore = 0;
		//double iters = allEdges.size()+1;
		double iters = Math.sqrt(allEdges.size())+1;
		for(int i=0;i<allEdges.size() && i<iters;i++) {
			Edge e=allEdges.get(i);
			if(e.getWeight()>0) {
				//System.out.println("Starting with edge: "+e.getPos1() +" - "+e.getPos2());
				initCut(e);
				boolean improvement = true;
				while(improvement) {
					heuristic1();
					improvement = heuristic2();
				}

				double score = calculateScore (cut);
				//System.out.println("Cut score: "+score);
				if(maxScore < score) {
					maxScore = score;
					copy(bestCut,cut);
				}
			}

		}
		copy(cut,bestCut);
	}
	private double calculateScore(boolean[] cut) {
		double answer = 0;
		for(Edge e:allEdges) {
			if(cut[e.getPos1()]!=cut[e.getPos2()]) {
				answer+=e.getWeight();
			}
		}
		return answer;
	}

	private void copy(boolean[] cutDest, boolean[] cutSource) {
		for(int i=0;i<cutDest.length;i++) {
			cutDest[i] = cutSource[i];
		}
	}

	/*private void randomizeCut() {
		Random r = new Random();
		for(int i=0;i<cut.length;i++) {
			cut[i] = r.nextBoolean();
			assigned[i]=true;
		}
	}*/
	private void initCut(Edge e) {
		boolean [] assigned = new boolean[graph.size()];
		Arrays.fill(assigned, false);
		assigned[e.getPos1()] = true;
		cut[e.getPos1()] = false;
		assigned[e.getPos2()] = true;
		cut[e.getPos2()] = true;
		int nAssigned = 2;
		while(nAssigned<cut.length) {
			Vertex vMax = null;
			double max = 0;
			boolean group = false;
			for(Vertex v:graph) {
				if(!assigned[v.getPos()]) {
					double diff = getAssignmentDiff(v,assigned);
					double absDiff = Math.abs(diff);
					if(vMax == null || absDiff>max) {
						max = absDiff;
						vMax = v;
						group = diff<0;
					}
				}
			}
			assigned[vMax.getPos()] = true;
			cut[vMax.getPos()] = group;
			nAssigned++;
		}
	}
	private void initCutRandom(Edge e) 
	{
		boolean [] assigned = new boolean[graph.size()];
		Arrays.fill(assigned, false);
		assigned[e.getPos1()] = true;
		cut[e.getPos1()] = false;
		assigned[e.getPos2()] = true;
		cut[e.getPos2()] = true;
		int nAssigned = 2;
		while(nAssigned<cut.length) {
			Vertex vMax = null;
			double max = 0;
			boolean group = false;
			for(Vertex v:graph) {
				if(!assigned[v.getPos()]) {
					double diff = getAssignmentDiff(v,assigned);
					double absDiff = Math.abs(diff);
					if(vMax == null || absDiff>max) {
						max = absDiff;
						vMax = v;
						group = diff<0;
					}
				}
			}
			assigned[vMax.getPos()] = true;
			cut[vMax.getPos()] = group;
			nAssigned++;
		}
	}
	private double getAssignmentDiff(Vertex v, boolean [] assigned) {
		int pos1 = v.getPos(); 
		List<Edge> edges = v.getEdges();
		double answer = 0;
		for(Edge e:edges) {
			int pos2 = e.getPos1();
			if(pos1 == pos2) {
				pos2 = e.getPos2();
			}
			if(assigned[pos2]){
				if(cut[pos2]) {
					answer += e.getWeight();
				} else {
					answer -= e.getWeight();
				}
			}
		}
		return answer;
	}

	private void heuristic1() {
		Vertex maxVertex ;
		do {
			maxVertex = null;
			double improvement = 0;
			for(Vertex v:graph) {
				double diff = getFlipDifference(v, cut);
				if(diff > improvement) {
					maxVertex = v;
					improvement = diff;
				}
			}
			if(maxVertex!=null) {
				flipVertex(maxVertex.getPos());
			}
		} while (maxVertex !=null);
	}
	private boolean heuristic2() {
		boolean totalImp = false;
		Edge maxEdge;
		do {
			maxEdge = null;

			double improvement = 0;
			for(Edge e:allEdges) {
				int pos1 = e.getPos1();
				int pos2 = e.getPos2();

				Vertex v1 = graph.get(pos1);
				Vertex v2 = graph.get(pos2);
				double diff1 = getFlipDifference(v1, cut);
				double diff2 = getFlipDifference(v2, cut);
				boolean same = cut[pos1]==cut[pos2];
				if(same) {
					diff1-=e.getWeight();
					diff2-=e.getWeight();
				} else {
					diff1+=e.getWeight();
					diff2+=e.getWeight();
				}
				if(diff1 + diff2 > improvement) {
					maxEdge = e;
					improvement = diff1+diff2;
				}
			}
			if(maxEdge!=null) {
				flipVertex(maxEdge.getPos1());
				flipVertex(maxEdge.getPos2());
				totalImp = true;
			}
		} while (maxEdge !=null);
		return totalImp;
	}
	private void flipVertex(int pos) {
		cut[pos] = !cut[pos];
	}
	private double getFlipDifference (Vertex v, boolean [] cut) 
	{
		int pos1 = v.getPos();
		boolean g1 = cut[pos1]; 
		List<Edge> edges = v.getEdges();
		double answer = 0;
		for(Edge e:edges) {
			int pos2 = e.getPos1();
			if(pos1 == pos2) {
				pos2 = e.getPos2();
			}
			boolean g2 = cut[pos2];
			if(g1 == g2) {
				answer += e.getWeight();
			} else {
				answer -= e.getWeight();
			}
		}
		return answer;
	}

	/**
	 * Implementation of the strategy for the new RefHap Algorithm
	 */
	public void calculateMaxCutStrategy2() 
	{
		boolean [] bestCut = new boolean [cut.length];
		double maxScore = 0;
		//double iters = Math.sqrt(allEdges.size())+1;
		for(int i=0;i<allEdges.size()/2;i++) 
		{
			Edge e=allEdges.get(i);
			if(e.getWeight()>0) 
			{
				initCutRandom(e);
				boolean improvement = true;
				int pr=1;
				while(improvement && pr<300) 
				{
					improvement = heuristic2();
				} pr++;

				double score = calculateScore (cut);
				if(maxScore < score) 
				{
					maxScore = score;
					copy(bestCut,cut);
				}
			}

		}
		copy(cut,bestCut);

	}
	/**
	 * Creates a partition picking random edges. 
	 * This partition may not be complete
	 * Some vertices can not be assigned
	 */
	private void heuristic3(double probNotExchange)
	{
	Vertex maxVertex ;
		do {
			maxVertex = null;
			double improvement = 0;
			for(Vertex v:graph)
			{
				double diff = getFlipDifference(v, cut);
				if(diff*(1-probNotExchange) > improvement) 
				{
					
					maxVertex = v;
					improvement = diff;
				}
			}
			if(maxVertex!=null) 
			{
				flipVertex(maxVertex.getPos());
			}
		} while (maxVertex !=null);
	}
	
	private boolean heuristic4(double probNotExchange) 
	{
		boolean totalImp = false;
		Edge maxEdge;
		do {
			maxEdge = null;

			double improvement = 0;
			for(Edge e:allEdges) {
				int pos1 = e.getPos1();
				int pos2 = e.getPos2();

				Vertex v1 = graph.get(pos1);
				Vertex v2 = graph.get(pos2);
				double diff1 = getFlipDifference(v1, cut);
				double diff2 = getFlipDifference(v2, cut);
				boolean same = cut[pos1]==cut[pos2];
				if(same) {
					diff1-=e.getWeight();
					diff2-=e.getWeight();
				} else {
					diff1+=e.getWeight();
					diff2+=e.getWeight();
				}
				if((diff1 + diff2)*(1-probNotExchange) > improvement) {
					maxEdge = e;
					improvement = diff1+diff2;
				}
			}
			if(maxEdge!=null) {
				flipVertex(maxEdge.getPos1());
				flipVertex(maxEdge.getPos2());
				totalImp = true;
			}
		} while (maxEdge !=null);
		return totalImp;
	}
	/*
	 * The MaxCut Strategy 3 reduces the running of the edges exchange and the vertex exchange strategies
	 */
	public void calculateMaxCutStrategy3() 
	{
		boolean [] bestCut = new boolean [cut.length];
		double maxScore = 0;
		//double iters = allEdges.size()+1;
		double iters = Math.sqrt(allEdges.size())+1;
		for(int i=0;i<allEdges.size() && i<iters;i++) {
			Edge e=allEdges.get(i);
			if(e.getWeight()>0) {
				//System.out.println("Starting with edge: "+e.getPos1() +" - "+e.getPos2());
				initCut(e);
				boolean improvement = true;
				while(improvement)
				{
					double prob = Math.random()*(0.05);
					heuristic3(prob);
					improvement = heuristic4(prob);
				}

				double score = calculateScore (cut);
				//System.out.println("Cut score: "+score);
				if(maxScore < score) {
					maxScore = score;
					copy(bestCut,cut);
				}
			}
		
	}

}
class Vertex 
{
	private int pos;
	private List<Edge> edges = new ArrayList<Edge>();
	public Vertex(int pos) 
	{
		super();
		this.pos = pos;
	}

	public void deleteEdges(Vertex v) 
	{
		// TODO Auto-generated method stub
		
		for(Edge e:edges)
		{
			if(e.getPos1()==pos || e.getPos2()==pos)
			{
				edges.remove(e);
			}
		}
		
	}

	/**
	 * @return the pos
	 */
	public int getPos() {
		return pos;
	}

	/**
	 * @return the edges
	 */
	public List<Edge> getEdges() {
		return edges;
	}

	public void addEdge (Edge e) {
		edges.add(e);
	}
}
class Edge {
	private int pos1;
	private int pos2;
	private double weight;
	public Edge(int pos1, int pos2, double weight) {
		super();
		this.pos1 = pos1;
		this.pos2 = pos2;
		this.weight = weight;
	}
	/**
	 * @return the weight
	 */
	public double getWeight() {
		return weight;
	}
	/**
	 * @param weight the weight to set
	 */
	public void setWeight(double weight) {
		this.weight = weight;
	}
	/**
	 * @return the pos1
	 */
	public int getPos1() {
		return pos1;
	}
	/**
	 * @return the pos2
	 */
	public int getPos2() {
		return pos2;
	}

}
class EdgesComparator implements Comparator<Edge> {

	@Override
	public int compare(Edge e1, Edge e2) {
		return (int)(e2.getWeight() -e1.getWeight());
	}
}
}