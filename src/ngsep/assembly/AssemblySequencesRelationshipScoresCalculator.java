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
	
	
	public boolean isUseIndels() {
		return useIndels;
	}
	public void setUseIndels(boolean useIndels) {
		this.useIndels = useIndels;
	}
	public int calculateScore(AssemblySequencesRelationship relationship, Distribution[] edgesStats) {
		double evProp = relationship.getEvidenceProportion();
		//return edge.getCoverageSharedKmers();
		//return edge.getRawKmerHits();
		double score = (relationship.getOverlap()+relationship.getWeightedCoverageSharedKmers())*evProp;
		//double score = (relationship.getOverlap()+relationship.getWeightedCoverageSharedKmers())*evProp*evProp;
		if(useIndels) {
			Distribution indelsKbp = edgesStats[5];
			NormalDistribution ikbp = new NormalDistribution(indelsKbp.getAverage(),2*indelsKbp.getVariance());
			double pValueIKBP = 1-ikbp.cumulative(relationship.getIndelsPerKbp());
			if(pValueIKBP>0.5) pValueIKBP = 0.5;
			score*=pValueIKBP;
		}
		return (int)Math.round(score);
	}
	public int calculateCost(AssemblySequencesRelationship relationship, Distribution[] edgesStats) {
		Distribution overlapTP = edgesStats[0];
		Distribution covTP = edgesStats[1];
		Distribution wCovTP = edgesStats[2];
		Distribution wCovPropTP = edgesStats[3];
		Distribution evPropTP = edgesStats[4];
		Distribution indelsKbpTP = edgesStats[5];
		//double prop = (double)edge.getCoverageSharedKmers()/edge.getOverlap();
		//int cost = (int)Math.round(1000*(1.5-prop));
		//return cost;
		NormalDistribution noTP = new NormalDistribution(overlapTP.getAverage(),overlapTP.getVariance()+1);
		NormalDistribution ncTP = new NormalDistribution(covTP.getAverage(),covTP.getVariance()+1);
		NormalDistribution nwcTP = new NormalDistribution(wCovTP.getAverage(),wCovTP.getVariance()+1);
		NormalDistribution nwcpTP = new NormalDistribution(wCovPropTP.getAverage(),wCovPropTP.getVariance()+0.0001);
		NormalDistribution evpTP = new NormalDistribution(evPropTP.getAverage(),evPropTP.getVariance()+0.0001);
		NormalDistribution ikbpTP = new NormalDistribution(indelsKbpTP.getAverage(),indelsKbpTP.getVariance()+1);
		//NormalDistribution niTP = new NormalDistribution(200,40000);
		int overlap;
		if(relationship instanceof AssemblyEdge) overlap = ((AssemblyEdge)relationship).getOverlap();
		else overlap = ((AssemblyEmbedded)relationship).getSequenceEvidenceEnd() - ((AssemblyEmbedded)relationship).getSequenceEvidenceStart();
		double pValueOTP = noTP.cumulative(overlap);
		//if(pValueOTP>0.5) pValueOTP = 1- pValueOTP;
		int cost1 = PhredScoreHelper.calculatePhredScore(pValueOTP);
		double pValueCTP = ncTP.cumulative(relationship.getCoverageSharedKmers());
		int cost2 = PhredScoreHelper.calculatePhredScore(pValueCTP);
		double pValueWCTP = nwcTP.cumulative(relationship.getWeightedCoverageSharedKmers());
		int cost3 = PhredScoreHelper.calculatePhredScore(pValueWCTP);
		double pValueWCPTP = nwcpTP.cumulative((double)relationship.getWeightedCoverageSharedKmers()/(overlap+1));
		int cost4 = PhredScoreHelper.calculatePhredScore(pValueWCPTP);
		double pValueEvProp = evpTP.cumulative(relationship.getEvidenceProportion());
		int cost5 = PhredScoreHelper.calculatePhredScore(pValueEvProp);
		double pValueIKBP = 1-ikbpTP.cumulative(relationship.getIndelsPerKbp());
		int cost6 = PhredScoreHelper.calculatePhredScore(pValueIKBP);
		double costD = 0;
		costD+=cost1;
		//cost += cost2;
		costD += cost3;
		costD += cost5;
		if(useIndels) costD += cost6;
		
		costD*=1000;
		int cost = (int)Math.min(1000000000, costD);
		
		//cost+= (int) (1000000*(1-pValueOTP)*(1-pValueCTP));
		cost+= (int) (1000*(1-pValueOTP)*(1-pValueWCTP));

		if( logRelationship(relationship)) System.out.println("CalculateCost. Pvalues "+pValueOTP+" "+pValueCTP+" "+pValueWCTP+" "+pValueWCPTP+" "+pValueEvProp+" "+pValueIKBP+" costs: "+cost1+" "+cost2+" "+cost3+" "+cost4+" "+cost5+" "+cost6+" cost: " +cost+ "Rel: "+relationship);
		
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
