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
package ngsep.assembly;

import java.util.Map;

import JSci.maths.statistics.NormalDistribution;
import ngsep.math.Distribution;
import ngsep.math.PhredScoreHelper;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class AssemblySequencesRelationshipScoresCalculator {
	private boolean useIndels = false;
	private int debugIdx = -1;
	private Map<Integer,Double> averageIKBPs;
	private double alphaIKBP = 0.001;
	//private double alphaIKBP = 0;
	
	
	
	public boolean isUseIndels() {
		return useIndels;
	}
	public void setUseIndels(boolean useIndels) {
		this.useIndels = useIndels;
	}
	
	
	public Map<Integer, Double> getAverageIKBPs() {
		return averageIKBPs;
	}
	public void setAverageIKBPs(Map<Integer, Double> averageIKBPs) {
		this.averageIKBPs = averageIKBPs;
	}
	public int calculateScore(AssemblySequencesRelationship relationship, NormalDistribution[] edgesDists, Map<Integer,Distribution> byLengthSumIKBPDists) {
		double evProp = relationship.getEvidenceProportion();
		double overlapProportion=relationship.getOverlap();
		double w1 = 1;
		double w2 = 1;
		
		if(relationship instanceof AssemblyEmbedded) {
			overlapProportion = 1;
		} else {
			AssemblyEdge edge = (AssemblyEdge) relationship;
			int l1 = edge.getVertex1().getRead().getLength();
			int l2 = edge.getVertex2().getRead().getLength();
			overlapProportion/=Math.max(l1, l2);
		}
		double maxIKBP = getMaxAverageIKBP (relationship);
		//return edge.getCoverageSharedKmers();
		//return edge.getRawKmerHits();
		double score = (relationship.getOverlap()*w1+relationship.getWeightedCoverageSharedKmers()*w2)*evProp;
		//if(logRelationship(relationship)) System.out.println("Relationship: "+relationship+" Evidence proportion: "+evProp+" score: "+score);
		//double score = relationship.getOverlap()*evProp+relationship.getWeightedCoverageSharedKmers()*Math.sqrt(overlapProportion);
		//double score = (relationship.getOverlap()+relationship.getWeightedCoverageSharedKmers())*evProp*evProp;
		if(useIndels) {
			NormalDistribution globalIKBP = edgesDists[5];
			double avg = Math.max(globalIKBP.getMean(), maxIKBP);
			NormalDistribution normalDistIkbp = new NormalDistribution(avg,Math.max(avg,globalIKBP.getVariance()));
			double pValueIKBP = 1-normalDistIkbp.cumulative(relationship.getIndelsPerKbp());
			int key = relationship.getLengthSum()/1000;
			Distribution byLength = byLengthSumIKBPDists.get(key);
			if(byLength!=null && byLength.getCount()>20) {
				avg = Math.max(byLength.getAverage(), maxIKBP);
				double normalVariance = Math.max(avg, byLength.getVariance());
				normalVariance = Math.min(normalVariance, avg*avg);
				NormalDistribution normalDist = new NormalDistribution(avg,normalVariance);
				pValueIKBP = 1-normalDist.cumulative(relationship.getIndelsPerKbp());
				//count = byLength.getCount();
				if(logRelationship(relationship)) System.out.println("Relationship: "+relationship+" key: "+key+" IKBP avg: "+byLength.getAverage()+" variance: "+normalDist.getVariance()+" pvalIkbp: "+pValueIKBP);
			} else if(logRelationship(relationship)) System.out.println("Relationship: "+relationship+" key: "+key+" pvalIkbp: "+pValueIKBP);
			if(pValueIKBP>alphaIKBP) pValueIKBP = 0.5;
			//pValueIKBP*=2;
			score*=pValueIKBP;
		}
		return (int)Math.round(score);
	}
	private double getMaxAverageIKBP(AssemblySequencesRelationship relationship) {
		Double avgIKBP1;
		Double avgIKBP2;
		if(relationship instanceof AssemblyEmbedded) {
			AssemblyEmbedded embedded = (AssemblyEmbedded)relationship;
			avgIKBP1 = averageIKBPs.get(embedded.getHostId());
			avgIKBP2 = averageIKBPs.get(embedded.getSequenceId());
		} else {
			AssemblyEdge edge = (AssemblyEdge) relationship;
			avgIKBP1 = averageIKBPs.get(edge.getVertex1().getSequenceIndex());
			avgIKBP2 = averageIKBPs.get(edge.getVertex2().getSequenceIndex());
		}
		if(avgIKBP1==null) avgIKBP1=1.0;
		if(avgIKBP2==null) avgIKBP2=1.0;
		return Math.max(avgIKBP1, avgIKBP2);
	}
	public int calculateCost(AssemblySequencesRelationship relationship, NormalDistribution[] edgesDists, Map<Integer,Distribution> byLengthSumIKBPDists) {
		NormalDistribution overlapD = edgesDists[0];
		NormalDistribution cskD = edgesDists[1];
		NormalDistribution wcskD = edgesDists[2];
		NormalDistribution wcskpD = edgesDists[3];
		NormalDistribution evPropD = edgesDists[4];
		NormalDistribution indelsKbpD = edgesDists[5];
		double maxIKBP = getMaxAverageIKBP(relationship);
		double avg = Math.max(indelsKbpD.getMean(), maxIKBP);
		NormalDistribution normalDistIkbp = new NormalDistribution(avg,Math.max(avg,indelsKbpD.getVariance()));
		//double prop = (double)edge.getCoverageSharedKmers()/edge.getOverlap();
		//int cost = (int)Math.round(1000*(1.5-prop));
		//return cost;
		//NormalDistribution niTP = new NormalDistribution(200,40000);
		int overlap = relationship.getOverlap();
		double cumulativeOverlap = overlapD.cumulative(overlap);
		//if(pValueOTP>0.5) pValueOTP = 1- pValueOTP;
		//int cost1 = PhredScoreHelper.calculatePhredScore(cumulativeOverlap);
		double cost1 = 100.0*(1-cumulativeOverlap);
		double cumulativeCSK = cskD.cumulative(relationship.getCoverageSharedKmers());
		double cost2 = 100.0*(1-cumulativeCSK);
		double cumulativeWCSK = wcskD.cumulative(relationship.getWeightedCoverageSharedKmers());
		double cost3 = 100.0*(1-cumulativeWCSK);
		//double pValueWCPTP = wcskpD.cumulative((double)relationship.getWeightedCoverageSharedKmers()/(overlap+1));
		int cost4 = PhredScoreHelper.calculatePhredScore(Math.min(0.05, cumulativeWCSK));
		double pValueEvProp = evPropD.cumulative(relationship.getEvidenceProportion());
		if(pValueEvProp>0.5) pValueEvProp = 0.5;
		int cost5 = PhredScoreHelper.calculatePhredScore(pValueEvProp);
		//double cost5 = 100.0*(1.0-relationship.getEvidenceProportion());
		double pValueIKBP = 1-normalDistIkbp.cumulative(relationship.getIndelsPerKbp());
		
		if(pValueIKBP>alphaIKBP) pValueIKBP = 0.5;
		int cost6 = PhredScoreHelper.calculatePhredScore(pValueIKBP);
		Distribution byLength = byLengthSumIKBPDists.get(relationship.getLengthSum()/1000);
		int cost7 = cost6;
		double pValueIKBP2 = pValueIKBP;
		if(byLength!=null && byLength.getCount()>20) {
			avg = Math.max(byLength.getAverage(), maxIKBP);
			double normalVariance = Math.max(avg, byLength.getVariance());
			normalVariance = Math.min(normalVariance, avg*avg);
			NormalDistribution normalDist = new NormalDistribution(avg,normalVariance);
			pValueIKBP2 = 1-normalDist.cumulative(relationship.getIndelsPerKbp());
			if( logRelationship(relationship)) System.out.println("CalculateCost. maxIKBP: "+maxIKBP+" final avg: "+avg+" variance group: "+byLength.getVariance()+" finalVar: "+normalVariance+" pval2: "+pValueIKBP2); 
			if(pValueIKBP2>alphaIKBP) pValueIKBP2 = 0.5;
			cost7 = PhredScoreHelper.calculatePhredScore(pValueIKBP2);
			//if(relationship.getIndelsPerKbp()>avg) cost7+=relationship.getIndelsPerKbp();
		}
		
		double costD = 0;
		costD+=cost1;
		//cost += cost2;
		//costD += cost3;
		costD += cost4;
		if(useIndels) costD += cost7;
		//if(useIndels) costD/=(pValueIKBP2+0.5);
		costD += cost5;
		//costD/=relationship.getEvidenceProportion();
		
		
		
		int cost = (int)(100.0*costD);
		
		//cost+= (int) (1000000*(1-pValueOTP)*(1-pValueCTP));
		//cost+= (int) (1000*(1-pValueOTP)*(1-pValueWCTP));

		if( logRelationship(relationship)) System.out.println("CalculateCost. Values "+cumulativeOverlap+" "+cumulativeWCSK+" "+pValueEvProp+" "+pValueIKBP+" "+pValueIKBP2+" costs: "+cost1+" "+cost3+" "+cost4+" "+cost5+" "+cost6+" "+cost7+" cost: " +cost+ " Rel: "+relationship);
		
		return cost;
	}
	private boolean logRelationship(AssemblySequencesRelationship relationship) {
		if (relationship instanceof AssemblyEdge) {
			AssemblyEdge edge = (AssemblyEdge)relationship;
			return edge.getVertex1().getSequenceIndex()==debugIdx || edge.getVertex2().getSequenceIndex()==debugIdx;
		} else if (relationship instanceof AssemblyEmbedded) {
			AssemblyEmbedded embedded = (AssemblyEmbedded)relationship;
			return embedded.getSequenceId() == debugIdx;
		}
		return false;
	}
}
