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
	private int debugIdx = -1;
	private Map<Integer,Double> averageIKBPVertices;
	private Map<Integer,Double> averageIKBPEmbedded;
	private double alphaIKBP = 0.01;
	//private double alphaIKBP = 0;
	
	
	public Map<Integer, Double> getAverageIKBPVertices() {
		return averageIKBPVertices;
	}
	public void setAverageIKBPVertices(Map<Integer, Double> averageIKBPVertices) {
		this.averageIKBPVertices = averageIKBPVertices;
	}
	public Map<Integer, Double> getAverageIKBPEmbedded() {
		return averageIKBPEmbedded;
	}
	public void setAverageIKBPEmbedded(Map<Integer, Double> averageIKBPEmbedded) {
		this.averageIKBPEmbedded = averageIKBPEmbedded;
	}
	public int calculateScore(AssemblySequencesRelationship relationship, NormalDistribution[] edgesDists) {
		NormalDistribution overlapD = edgesDists[0];
		double evProp = relationship.getEvidenceProportion();
		double overlapRel = (double)relationship.getOverlap()/overlapD.getMean();
		if(overlapRel>1) overlapRel = 1;
		double maxIKBP = getMaxAverageIKBP (relationship);
		//return edge.getCoverageSharedKmers();
		//return edge.getRawKmerHits();
		//double score = relationship.getWeightedCoverageSharedKmers();
		//double score = relationship.getWeightedCoverageSharedKmers()*evProp;
		//double score = (relationship.getOverlap()*w1+relationship.getWeightedCoverageSharedKmers()*w2)*evProp;
		//double score = (0.001*relationship.getOverlap())*relationship.getWeightedCoverageSharedKmers();
		//double score = overlapRel*relationship.getWeightedCoverageSharedKmers();
		//double score = overlapRel*relationship.getWeightedCoverageSharedKmers()*evProp/Math.sqrt(1+relationship.getIndelsPerKbp());
		double score = Math.sqrt(relationship.getOverlap())*(0.01*relationship.getWeightedCoverageSharedKmers())*Math.sqrt(evProp);
		//double score = (0.001*relationship.getOverlap())*relationship.getWeightedCoverageSharedKmers()*evProp;
		//if(logRelationship(relationship)) System.out.println("Relationship: "+relationship+" Evidence proportion: "+evProp+" score: "+score);
		//double score = relationship.getOverlap()*evProp+relationship.getWeightedCoverageSharedKmers()*Math.sqrt(overlapProportion);
		//double score = (relationship.getOverlap()+relationship.getWeightedCoverageSharedKmers())*evProp*evProp;
		
		/*NormalDistribution globalIKBP = edgesDists[5];
		double globalMean = globalIKBP.getMean();
		double avg = Math.max(globalMean, maxIKBP);
		//Average is not controlled in the chimera detection for embedded relationships
		//if(relationship instanceof AssemblyEmbedded) avg = Math.min(globalMean*2, avg);
		NormalDistribution normalDistIkbp = new NormalDistribution(avg,Math.max(avg,globalIKBP.getVariance()));
		double pValueIKBP = 1-normalDistIkbp.cumulative(relationship.getIndelsPerKbp());
		if(pValueIKBP>0.0001) pValueIKBP = 0.5;
		//pValueIKBP*=2;
		score*=pValueIKBP;*/
		
		return (int)Math.round(score);
	}
	private double getMaxAverageIKBP(AssemblySequencesRelationship relationship) {
		Double avgIKBP1;
		Double avgIKBP2;
		if(relationship instanceof AssemblyEmbedded) {
			AssemblyEmbedded embedded = (AssemblyEmbedded)relationship;
			avgIKBP1 = averageIKBPEmbedded.get(embedded.getHostId());
			avgIKBP2 = averageIKBPEmbedded.get(embedded.getSequenceId());
		} else {
			AssemblyEdge edge = (AssemblyEdge) relationship;
			avgIKBP1 = averageIKBPVertices.get(edge.getVertex1().getUniqueNumber());
			avgIKBP2 = averageIKBPVertices.get(edge.getVertex2().getUniqueNumber());
		}
		if(avgIKBP1==null) avgIKBP1=1.0;
		if(avgIKBP2==null) avgIKBP2=1.0;
		return Math.max(avgIKBP1, avgIKBP2);
	}
	public int calculateCost(AssemblySequencesRelationship relationship, NormalDistribution[] edgesDists) {
		NormalDistribution overlapD = edgesDists[0];
		NormalDistribution cskD = edgesDists[1];
		NormalDistribution wcskD = edgesDists[2];
		NormalDistribution overlapSD = edgesDists[3];
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
		int cost1 = PhredScoreHelper.calculatePhredScore(Math.min(0.5, cumulativeOverlap));
		//double cost1 = 20.0*(1-cumulativeOverlap);
		
		double cumulativeCSK = cskD.cumulative(relationship.getCoverageSharedKmers());
		double cost2 = 100.0*(1-cumulativeCSK);
		double cumulativeWCSK = wcskD.cumulative(relationship.getWeightedCoverageSharedKmers());
		//double cost3 = 100.0*(1-cumulativeWCSK);
		
		int cost3 = PhredScoreHelper.calculatePhredScore(Math.min(0.5, cumulativeWCSK));
		
		//TODO: Save overlapSD for embedded 
		double pValueOverlapSD = 0.5;
		if(relationship instanceof AssemblyEdge) pValueOverlapSD = 1-overlapSD.cumulative(((AssemblyEdge)relationship).getOverlapStandardDeviation());
		int cost4 = PhredScoreHelper.calculatePhredScore(Math.min(0.05, pValueOverlapSD));
		
		double cumulativeEvProp = evPropD.cumulative(relationship.getEvidenceProportion());
		int cost5 = PhredScoreHelper.calculatePhredScore(Math.min(0.05, cumulativeEvProp));
		//double cost5 = 100.0*(1.0-relationship.getEvidenceProportion());
		double pValueIKBP = 1-normalDistIkbp.cumulative(relationship.getIndelsPerKbp());
		
		//if(pValueIKBP>alphaIKBP) pValueIKBP = 0.5;
		//int cost6 = PhredScoreHelper.calculatePhredScore(pValueIKBP);
		double cost6 = Math.max(indelsKbpD.getMean(), relationship.getIndelsPerKbp());
		
		double costD = 0;
		costD+=cost1;
		//cost += cost2;
		costD += cost3;
		//costD += cost4;
		costD += cost5;
		costD += cost6;
		//costD/=relationship.getEvidenceProportion();
		
		int cost = (int)(1000.0*costD);
		
		cost+= (int) (100.0*(1-cumulativeOverlap));
		//cost+= (int) (1000*(1-pValueOTP)*(1-pValueWCTP));

		if( logRelationship(relationship)) System.out.println("CalculateCost. Values "+cumulativeOverlap+" "+cumulativeWCSK+" "+cumulativeEvProp+" "+pValueIKBP+" costs: "+cost1+" "+cost3+" "+cost4+" "+cost5+" "+cost6+" cost: " +cost+ " Rel: "+relationship);
		
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
