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

import java.util.ArrayList;
import java.util.List;

import JSci.maths.statistics.NormalDistribution;
import ngsep.math.Distribution;

/**
 * 
 * @author Jorge Duitama
 *
 */
public class AssemblySequencesRelationshipFilter {
	private int debugIdx = -1;
	//Initial edge filtering
	public void filterEdgesAndEmbedded(AssemblyGraph graph, double minScoreProportionEdges) {
		Distribution lengthsDistribution = new Distribution(0, graph.getSequenceLength(0), 1);
		int n = graph.getNumSequences();
		for(int i=0;i<n;i++) lengthsDistribution.processDatapoint(graph.getSequenceLength(i));
		int medianLength = graph.getMedianLength();
		System.out.println("Median read length: "+medianLength);
		NormalDistribution distLengths = new NormalDistribution(medianLength,medianLength);
		int numEmbedded = 0;
		for (int seqId = n-1; seqId >=0; seqId--) {
			AssemblyVertex vS = graph.getVertex(seqId, true);
			AssemblyVertex vE = graph.getVertex(seqId, false);
			if(vS==null || vE==null) continue;
			
			int [] bestValues = filterEdges(graph, seqId, medianLength, minScoreProportionEdges);
			int maxScore = Math.max(bestValues[0], bestValues[1]);
			if(bestValues[0]==0 || bestValues[1]==0) {
				//if(maxScore>0) System.out.println("Zero score on one side for sequence: "+seqId+" "+graph.getSequence(seqId).getName()+ ". Scores: "+bestValues[0]+" "+bestValues[1]+" Removing vertices");
				//graph.removeVertices(seqId);
				maxScore = 0;
			}
			//if(filterEmbeddedByCost(graph, seqId, medianLength, Math.min(bestValues[2], bestValues[3]))) {
			if(filterEmbeddedByScore(graph, seqId, medianLength, distLengths, maxScore)) {
				numEmbedded++;
				//List<AssemblyEdge> edges = graph.getEdgesBySequenceId(seqId);
				//for(AssemblyEdge edge:edges) edge.setCost(10*edge.getCost());
			}
			if(seqId == debugIdx && graph.getVertex(seqId, true)!=null) System.out.println("EdgesAndEmbeddedFiltering. Edges after processing sequence: "+graph.getEdges(vS).size()+" "+graph.getEdges(vE).size()+" self edge: "+graph.getSameSequenceEdge(seqId));
		}
		System.out.println("Filtered edges and embedded. Final number of embedded sequences: "+numEmbedded);
		for (int seqId = n-1; seqId >=0; seqId--) {
			AssemblyVertex vS = graph.getVertex(seqId, true);
			AssemblyVertex vE = graph.getVertex(seqId, false);
			if(vS==null || vE==null) continue;
			int nS = graph.getEdges(vS).size()-1;
			int nE = graph.getEdges(vE).size()-1;
			int nM = Math.min(nS, nE); 
			if(seqId == debugIdx) System.out.println("Final edges for sequence: "+seqId+" "+graph.getSequence(seqId).getName()+" START "+graph.getEdges(vS)+" END "+graph.getEdges(vE));
			if(nM==0 /*|| (nM==1 && Math.max(nS, nE)>20) */) {
				//System.out.println("Disbalanced number of edges for sequence: "+seqId+" "+graph.getSequence(seqId).getName()+ ". Values: "+nS+" "+nE+" Removing vertices");
				//graph.removeVertices(seqId);
			}
		}
		//graph.pruneEmbeddedSequences();
		//System.out.println("Prunned embedded sequences.");
		//if(debugIdx>=0) System.out.println("EdgesAndEmbeddedFiltering. Final number of edges: "+graph.getEdges(graph.getVertex(debugIdx, true)).size()+" "+graph.getEdges(graph.getVertex(debugIdx, false)).size());
		//filterEdgesCloseRelationships();
		//System.out.println("Searched vertex: "+graph.getVertex(debugIdx, false));
	}
	
	private int [] filterEdges (AssemblyGraph graph, int sequenceId,int medianLength, double minScoreProportionEdges) {
		AssemblyVertex vS = graph.getVertex(sequenceId, true);
		AssemblyVertex vE = graph.getVertex(sequenceId, false);
		List<AssemblyEdge> edgesS = new ArrayList<AssemblyEdge>();
		if(vS!=null) edgesS.addAll(graph.getEdges(vS));
		
		
		int maxScoreS = 0;
		int minCostS = 100000000;
		for(AssemblyEdge edge: edgesS) {
			if(edge.isSameSequenceEdge()) continue;
			//int connectingLength = getSequenceLength(edge.getConnectingVertex(vS).getSequenceIndex());
			/*if(connectingLength<1.2*sequenceLength && edge.getOverlap() >0.8*sequenceLength)*/
			maxScoreS = Math.max(maxScoreS, edge.getScore());
			minCostS = Math.min(minCostS, edge.getCost());
			 
		}
		if(sequenceId == debugIdx) System.out.println("EdgesFiltering. Initial edges start "+edgesS.size()+" Max score start: "+maxScoreS+" minCost: "+minCostS);
		List<AssemblyEdge> edgesE = new ArrayList<AssemblyEdge>();
		if(vE!=null) edgesE.addAll(graph.getEdges(vE));
		int maxScoreE = 0;
		int minCostE = 100000000;
		for(AssemblyEdge edge: edgesE) {
			if(edge.isSameSequenceEdge()) continue;
			//int connectingLength = getSequenceLength(edge.getConnectingVertex(vE).getSequenceIndex());
			/*if(connectingLength<1.5*sequenceLength && edge.getOverlap() >0.8*sequenceLength)*/ 
			maxScoreE = Math.max(maxScoreE, edge.getScore());
			minCostE = Math.min(minCostE, edge.getCost());
		}
		if(sequenceId == debugIdx) System.out.println("EdgesFiltering. Initial edges end "+edgesE.size()+" Max score end: "+maxScoreE+" min cost: "+minCostE);
		//double limitS = minScoreProportionEdges*Math.max(maxScoreS, maxScoreE);
		double limitS = minScoreProportionEdges*maxScoreS;
		for(AssemblyEdge edge: edgesS) {
			if(edge.isSameSequenceEdge()) continue;
			double score = edge.getScore();
			if(sequenceId == debugIdx) System.out.println("EdgesFiltering. Next edge start "+edge+" limit: "+limitS);
			if(score < maxScoreS && score < limitS) {
				if(sequenceId == debugIdx) System.out.println("EdgesFiltering. Removing edge: "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber());
				graph.removeEdge(edge);
			}
		}
		//double limitE = minScoreProportionEdges*Math.max(maxScoreS, maxScoreE);
		double limitE = minScoreProportionEdges*maxScoreE;
		for(AssemblyEdge edge: edgesE) {
			if(edge.isSameSequenceEdge()) continue;
			double score = edge.getScore();
			if(sequenceId == debugIdx) System.out.println("EdgesFiltering. Next edge end "+edge+" limit: "+limitE);
			if(score < maxScoreE && score < limitE) {
				if(sequenceId == debugIdx) System.out.println("EdgesFiltering. Removing edge: "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber());
				graph.removeEdge(edge);
			}
		}
		if(sequenceId == debugIdx) System.out.println("EdgesFiltering. Final number of edges: "+graph.getEdges(vS).size()+" "+graph.getEdges(vE).size());
		int [] answer = {maxScoreS,maxScoreE,minCostS,minCostE};
		return answer;
	}

	private boolean filterEmbeddedByCost(AssemblyGraph graph, int sequenceId,int medianLength, int minCostEdges) {
		int sequenceLength = graph.getSequenceLength(sequenceId);
		double medianRelationship = 1.0*sequenceLength/(double)medianLength;
		double costInflationEmbedded = Math.max(1, 1.0/medianRelationship);
		//if(minScoreProportionEmbedded<0.5) minScoreProportionEmbedded = 0.5;
		double costLimit = minCostEdges*costInflationEmbedded;
		//double costLimit = minCostEdges;
		if(sequenceId == debugIdx) System.out.println("min cost edges: "+minCostEdges+" median rel: "+medianRelationship+" cost inflation: "+costInflationEmbedded+" cost limit: "+costLimit);
		List<AssemblyEmbedded> embeddedList= new ArrayList<AssemblyEmbedded>();
		embeddedList.addAll(graph.getEmbeddedBySequenceId(sequenceId));
		if(embeddedList.size()==0) return false;
		AssemblyEmbedded embeddedMin = null;
		double minCostEmbedded = -1;
		for(AssemblyEmbedded embedded:embeddedList) {
			int cost = embedded.getCost();
			if(sequenceId == debugIdx) System.out.println("Assembly graph. Next embedded "+embedded+" cost "+cost+" evProp: "+embedded.getEvidenceProportion()+" Indels: "+embedded.getNumIndels()+" IKBP: "+embedded.getIndelsPerKbp());
			
			if(embeddedMin==null || minCostEmbedded>cost) {
				minCostEmbedded = cost;
				embeddedMin = embedded;
			}
		}
		
		if(minCostEmbedded>costLimit) {
			//Replace embedded relationships with edges to make the sequence not embedded
			for(AssemblyEmbedded embedded:embeddedList) {
				graph.removeEmbedded(embedded);
				if(sequenceId == debugIdx) System.out.println("Adding edge replacing embedded "+embedded.getHostId()+" limits: "+embedded.getHostStart()+" "+embedded.getHostEnd()+" host length: "+graph.getSequenceLength(embedded.getHostId())+"score: "+embedded.getScore());
				addEdgeFromEmbedded(graph, embedded);
			}
			return false;
		} else {
			if(sequenceId == debugIdx) System.out.println("Sequence is embedded ");
			for(AssemblyEmbedded embedded:embeddedList) {
				if(embedded!=embeddedMin) {
					graph.removeEmbedded(embedded);
					//if(sequenceId == debugIdx) System.out.println("Assembly graph. Removed embedded host: "+embedded.getHostId()+" Embedded relations: "+embeddedMapBySequence.get(sequenceId)+" is embedded: "+isEmbedded(sequenceId));
				}
			}
			return true;
		}
	}
	private boolean filterEmbeddedByScore(AssemblyGraph graph, int sequenceId,int medianLength, NormalDistribution distLengths, int maxScoreEdges) {
		int sequenceLength = graph.getSequenceLength(sequenceId);
		double medianRelationship = 1.0*sequenceLength/(double)medianLength;
		//double minScoreProportionEmbedded = 0.8;
		//double cumulative = lengthsDistribution.getCumulativeCount(sequenceLength)/lengthsDistribution.getCount();
		double minScoreProportionEmbedded = Math.min(0.8, 0.5*medianRelationship);
		//double minScoreProportionEmbedded = Math.min(0.8, distLengths.cumulative(sequenceLength));
		//double minScoreProportionEmbedded = 0.8*cumulative;
		//This can be improved if the graph is completely calculated
		if(minScoreProportionEmbedded<0.5) minScoreProportionEmbedded = 0.5;
		
		double scoreLimit = minScoreProportionEmbedded*maxScoreEdges;
		
		List<AssemblyEmbedded> embeddedList= new ArrayList<AssemblyEmbedded>();
		embeddedList.addAll(graph.getEmbeddedBySequenceId(sequenceId));
		if(sequenceId == debugIdx) System.out.println("max score edges: "+maxScoreEdges+" minscoreprop: "+minScoreProportionEmbedded+" score limit: "+scoreLimit+" num embedded: "+embeddedList.size());
		if(embeddedList.size()==0) return false;
		AssemblyEmbedded embeddedMax = null;
		double maxScoreEmbedded = -1;
		for(AssemblyEmbedded embedded:embeddedList) {
			int score = embedded.getScore();
			if(sequenceId == debugIdx) System.out.println("Assembly graph. Next embedded "+embedded);
			
			if(embeddedMax==null || maxScoreEmbedded<score) {
				maxScoreEmbedded = score;
				embeddedMax = embedded;
			}
		}
		
		if(maxScoreEmbedded<scoreLimit) {
			//Replace embedded relationships with edges to make the sequence not embedded
			for(AssemblyEmbedded embedded:embeddedList) {
				graph.removeEmbedded(embedded);
				//if(sequenceId == debugIdx) System.out.println("Adding edge replacing embedded "+embedded);
				//addEdgeFromEmbedded(graph, embedded);
			}
			return false;
		} else {
			if(sequenceId == debugIdx) System.out.println("Sequence is embedded ");
			for(AssemblyEmbedded embedded:embeddedList) {
				if(embedded!=embeddedMax) {
					graph.removeEmbedded(embedded);
					//if(sequenceId == debugIdx) System.out.println("Assembly graph. Removed embedded host: "+embedded.getHostId()+" Embedded relations: "+embeddedMapBySequence.get(sequenceId)+" is embedded: "+isEmbedded(sequenceId));
				}
			}
			return true;
		}
	}
	

	private void addEdgeFromEmbedded(AssemblyGraph graph, AssemblyEmbedded embedded) {
		int distanceStart = embedded.getHostStart();
		int distanceEnd = graph.getSequenceLength(embedded.getHostId())-embedded.getHostEnd();
		AssemblyVertex vertexHost=null;
		AssemblyVertex vertexEmbedded=null;
		if(distanceStart<0.5*distanceEnd) {
			vertexHost = graph.getVertex(embedded.getHostId(), true);
			vertexEmbedded = graph.getVertex(embedded.getSequenceId(), embedded.isReverse());
		} else if (distanceEnd<0.5*distanceStart) {
			vertexHost = graph.getVertex(embedded.getHostId(), false);
			vertexEmbedded = graph.getVertex(embedded.getSequenceId(), !embedded.isReverse());
		}
		if(vertexHost==null || vertexEmbedded==null) return;
		int overlap = Math.max(embedded.getCoverageSharedKmers(), embedded.getSequenceEvidenceEnd()-embedded.getSequenceEvidenceStart());
		AssemblyEdge edge = new AssemblyEdge(vertexHost, vertexEmbedded, overlap);
		edge.setOverlapStandardDeviation(embedded.getHostStartStandardDeviation());
		edge.setWeightedCoverageSharedKmers(embedded.getWeightedCoverageSharedKmers());
		edge.setCoverageSharedKmers(embedded.getCoverageSharedKmers());
		edge.setNumSharedKmers(embedded.getNumSharedKmers());
		edge.setRawKmerHits(embedded.getRawKmerHits());
		edge.setRawKmerHitsSubjectStartSD(embedded.getRawKmerHitsSubjectStartSD());
		edge.setVertex1EvidenceStart(embedded.getHostEvidenceStart());
		edge.setVertex1EvidenceEnd(embedded.getHostEvidenceEnd());
		edge.setVertex2EvidenceStart(embedded.getSequenceEvidenceStart());
		edge.setVertex2EvidenceEnd(embedded.getSequenceEvidenceEnd());
		edge.setNumIndels(embedded.getNumIndels());
		edge.setNumMismatches(embedded.getNumMismatches());
		edge.setScore(embedded.getScore());
		edge.setCost(embedded.getCost());
		graph.addEdge(edge);
	}
}
