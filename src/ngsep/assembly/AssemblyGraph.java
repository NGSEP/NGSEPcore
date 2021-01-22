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
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import ngsep.math.Distribution;
import ngsep.sequences.QualifiedSequence;

/**
 * @author Jorge Duitama
 * @author Juan Camilo Bojaca
 * @author David Guevara
 */
public class AssemblyGraph {
	
	/**
	 * Sequences to build the graph. The index of each sequence is the unique identifier
	 */
	private List<QualifiedSequence> sequences;
	
	/**
	 * Sum of lengths from sequence 0 to i
	 */
	private long [] cumulativeReadLength;
	/**
	 * Start vertices indexed by sequence index
	 */
	private Map<Integer,AssemblyVertex> verticesStart;
	/**
	 * End vertices indexed by sequence index
	 */
	private Map<Integer,AssemblyVertex> verticesEnd;
	/**
	 * All vertices indexed by unique number
	 */
	private Map<Integer,AssemblyVertex> verticesByUnique;
	
	private Map<Integer,List<AssemblyEdge>> edgesMap = new HashMap<>();
	
	// Map with sequence ids as keys and embedded sequences as objects
	private Map<Integer, List<AssemblyEmbedded>> embeddedMapByHost = new HashMap<>();
	
	
	private Map<Integer, List<AssemblyEmbedded>> embeddedMapBySequence = new HashMap<>();

	private List<List<AssemblyEdge>> paths = new ArrayList<List<AssemblyEdge>>();
	
	private int numEdges = 0;
	
	/**
	 * Optional attribute to store the sequence names. Useful for the gold standard graph
	 */
	private List<String> readNames;

	public AssemblyGraph(List<QualifiedSequence> sequences) {
		int n = sequences.size();
		this.sequences = Collections.unmodifiableList(sequences);
		verticesStart = new HashMap<>(n);
		verticesEnd = new HashMap<>(n);
		verticesByUnique = new HashMap<>(n);
		edgesMap = new HashMap<>(n);
		cumulativeReadLength = new long [n];
		for (int i=0;i<sequences.size();i++) {
			QualifiedSequence seq = sequences.get(i);
			int length = seq.getLength();
			cumulativeReadLength[i]=length;
			if(i>0) cumulativeReadLength[i]+=cumulativeReadLength[i-1];
			AssemblyVertex vS = new AssemblyVertex(seq, true, i);
			verticesStart.put(i,vS);
			verticesByUnique.put(vS.getUniqueNumber(), vS);
			edgesMap.put(vS.getUniqueNumber(), new ArrayList<>());
			AssemblyVertex vE = new AssemblyVertex(seq, false, i);
			verticesEnd.put(i,vE);
			verticesByUnique.put(vE.getUniqueNumber(), vE);
			edgesMap.put(vE.getUniqueNumber(), new ArrayList<>());
			AssemblyEdge edge = new AssemblyEdge(vS, vE, length);
			edge.setAverageOverlap(length);
			edge.setMedianOverlap(length);
			edge.setFromLimitsOverlap(length);
			edge.setCoverageSharedKmers(length);
			edge.setWeightedCoverageSharedKmers(length);
			edge.setNumSharedKmers(length);
			edge.setOverlapStandardDeviation(0);
			edge.setVertex1EvidenceStart(0);
			edge.setVertex1EvidenceEnd(length-1);
			edge.setVertex2EvidenceStart(0);
			edge.setVertex2EvidenceEnd(length-1);
			addEdge(edge);
		}
	}
	
	//Modifiers
	
	public void addEdge(AssemblyEdge edge) {
		edgesMap.get(edge.getVertex1().getUniqueNumber()).add(edge);
		edgesMap.get(edge.getVertex2().getUniqueNumber()).add(edge);
		numEdges++;
	}
	
	public void removeEdge (AssemblyEdge edge) {
		List<AssemblyEdge> edges1 = edgesMap.get(edge.getVertex1().getUniqueNumber()); 
		if(edges1!=null) edges1.remove(edge);
		List<AssemblyEdge> edges2 = edgesMap.get(edge.getVertex2().getUniqueNumber());
		if(edges2!=null) edges2.remove(edge);
		numEdges--;
	}
	
	public void removeVertices(int sequenceId) {
		removeEdges(sequenceId);
		AssemblyVertex v1 = getVertex(sequenceId, true);
		AssemblyVertex v2 = getVertex(sequenceId, false);
		edgesMap.remove(v1.getUniqueNumber());
		edgesMap.remove(v2.getUniqueNumber());
		verticesStart.remove(sequenceId);
		verticesEnd.remove(sequenceId);
		verticesByUnique.remove(v1.getUniqueNumber());
		verticesByUnique.remove(v2.getUniqueNumber());
	}
	
	private void removeEdges(int sequenceId) {
		AssemblyVertex v1 = getVertex(sequenceId, true);
		List<AssemblyEdge> toRemove = new ArrayList<AssemblyEdge>();
		if(v1!=null) {
			List<AssemblyEdge> edges1 = edgesMap.get(v1.getUniqueNumber()); 
			if(edges1!=null) toRemove.addAll(edges1);
		}
		AssemblyVertex v2 = getVertex(sequenceId, false);
		if(v2!=null) {
			List<AssemblyEdge> edges2 = edgesMap.get(v2.getUniqueNumber());
			if(edges2!=null) {
				for(AssemblyEdge edge:edges2) {
					if(!edge.isSameSequenceEdge()) toRemove.add(edge);
				}
			}
		}
		for(AssemblyEdge edge:toRemove) {
			removeEdge(edge);
		}
	}
	
	public void addEmbedded(AssemblyEmbedded embeddedObject) {
		List<AssemblyEmbedded> list = embeddedMapByHost.computeIfAbsent(embeddedObject.getHostId(), key -> new LinkedList<>());
		list.add(embeddedObject);
		List<AssemblyEmbedded> list2 = embeddedMapBySequence.computeIfAbsent(embeddedObject.getSequenceId(), key -> new LinkedList<>());
		list2.add(embeddedObject);	
	}
	
	public void removeEmbedded (AssemblyEmbedded embeddedObject) {
		int hostId = embeddedObject.getHostId();
		embeddedMapByHost.get(hostId).remove(embeddedObject);
		if(embeddedMapByHost.get(hostId).size()==0) embeddedMapByHost.remove(hostId);
		int seqId = embeddedObject.getSequenceId();
		embeddedMapBySequence.get(seqId).remove(embeddedObject);
		if(embeddedMapBySequence.get(seqId).size()==0) embeddedMapBySequence.remove(seqId);
	}
	
	private void removeEmbeddedRelations(int sequenceId) {
		List<AssemblyEmbedded> embeddedList = new ArrayList<AssemblyEmbedded>();
		List<AssemblyEmbedded> emb = embeddedMapByHost.get(sequenceId);
		if(emb!=null) embeddedList.addAll(emb);
		emb = embeddedMapBySequence.get(sequenceId);
		if(emb!=null) embeddedList.addAll(emb);
		for(AssemblyEmbedded embedded: embeddedList) {
			removeEmbedded(embedded);
		}
	}
	
	public void pruneEmbeddedSequences() {
		for(int i:embeddedMapBySequence.keySet()) {
			if(verticesStart.get(i)!=null) {
				removeVertices(i);
			}
		}
	}

	/**
	 * @return the sequences
	 */
	public List<QualifiedSequence> getSequences() {
		return sequences;
	}
	public QualifiedSequence getSequence(int sequenceIdx) {
		return sequences.get(sequenceIdx);
	}
	public int getSequenceLength(int sequenceIdx) {
		return sequences.get(sequenceIdx).getLength();
	}
	public long getCumulativeLength(int sequenceIdx) {
		return cumulativeReadLength[sequenceIdx];
	}
	public int getNumSequences () {
		return sequences.size();
	}
	
	public List<String> getReadNames() {
		return readNames;
	}
	
	public int getMedianLength() {
		int n = getNumSequences();
		List<Integer> lengths = new ArrayList<Integer>(n);
		for(QualifiedSequence seq:sequences) lengths.add(seq.getLength());
		Collections.sort(lengths);
		return lengths.get(n/2);
	}

	

	public AssemblyVertex getVertex(int indexSequence, boolean start) {
		if(start) return verticesStart.get(indexSequence);
		return verticesEnd.get(indexSequence);
	}
	
	public AssemblyVertex getVertexByUniqueId(int uniqueId) {
		return verticesByUnique.get(uniqueId);
	}

	/**
	 * Return the list of embedded sequences for the given read
	 * 
	 * @param index of the read having embedded sequences
	 * @return list of embedded sequences
	 */
	public List<AssemblyEmbedded> getEmbeddedByHostId(int hostIndex) {
		List<AssemblyEmbedded> answer = embeddedMapByHost.get(hostIndex);
		if(answer == null) return new ArrayList<AssemblyEmbedded>();
		return answer;
	}
	
	/**
	 * 
	 * @param seqIndex
	 * @return List<AssemblyEmbedded> Sequences where this is embedded
	 */
	public List<AssemblyEmbedded> getEmbeddedBySequenceId(int seqIndex) {
		List<AssemblyEmbedded> answer = embeddedMapBySequence.get(seqIndex);
		if(answer == null) return new ArrayList<AssemblyEmbedded>();
		return answer;
	}
	
	public boolean isEmbedded(int sequenceId) {
		return embeddedMapBySequence.get(sequenceId)!=null;
	}
	public int getEmbeddedCount () {
		return embeddedMapBySequence.size();
	}

	public synchronized void addPath(List<AssemblyEdge> path) {
		paths.add(path);
		for(AssemblyEdge edge:path) edge.setLayoutEdge(true);
	}
	
	public List<AssemblyVertex> getVertices() {
		List<AssemblyVertex> vertices = new ArrayList<>();
		vertices.addAll(verticesByUnique.values());
		return vertices;
	}
	
	

	public int getNumEdges() {
		return numEdges;
	}

	/**
	 * @return the edges
	 */
	public List<AssemblyEdge> getEdges() {
		List<AssemblyEdge> edges = new ArrayList<>();
		for(AssemblyVertex v:verticesByUnique.values()) {
			List<AssemblyEdge> edgesVertex = edgesMap.get(v.getUniqueNumber());
			if(edgesVertex==null) continue;
			for(AssemblyEdge edge:edgesVertex) {
				//Avoid adding twice the same edge
				if(edge.getVertex1()==v) edges.add(edge);
			}
		}
		return edges;
	}
	
	public List<AssemblyEdge> getEdges(AssemblyVertex vertex) {
		return edgesMap.get(vertex.getUniqueNumber());
	}
	/**
	 * Returns the edge connecting the vertices of the given sequence id
	 * @param sequenceId
	 * @return AssemblyEdge
	 */
	public AssemblyEdge getSameSequenceEdge(int sequenceId) {
		AssemblyVertex vertex = verticesStart.get(sequenceId);
		if(vertex == null) return null;
		return getSameSequenceEdge(vertex);
	}
	/**
	 * Returns the edge connecting the given vertex with the corresponding vertex in the same sequence
	 * @param vertex
	 * @return AssemblyEdge
	 */
	public AssemblyEdge getSameSequenceEdge(AssemblyVertex vertex) {
		List<AssemblyEdge> edges = edgesMap.get(vertex.getUniqueNumber());
		for(AssemblyEdge edge:edges) {
			if(edge.getVertex1()==vertex && edge.getVertex2().getRead()==vertex.getRead()) {
				return edge;
			}
			if(edge.getVertex2()==vertex && edge.getVertex1().getRead()==vertex.getRead()) {
				return edge;
			}
		}
		throw new RuntimeException("Same sequence edge not found for vertex: "+vertex);
	}
	/**
	 * Searches for an edge between the given vertices
	 * @param v1
	 * @param v2
	 * @return
	 */
	public AssemblyEdge getEdge(AssemblyVertex v1, AssemblyVertex v2) {
		List<AssemblyEdge> edgesV1 = edgesMap.get(v1.getUniqueNumber());
		if(edgesV1 == null) return null;
		for(AssemblyEdge edge:edgesV1) {
			if(edge.getConnectingVertex(v1)==v2) return edge;
		}
		return null;
	}

	/**
	 * @return the paths
	 */
	public List<List<AssemblyEdge>> getPaths() {
		return paths;
	}

	public void updateVertexDegrees () {
		for (int v:edgesMap.keySet()) {
			AssemblyVertex vertex = verticesByUnique.get(v);
			vertex.setDegreeUnfilteredGraph(edgesMap.get(v).size());
		}
	}
	/**
	 * Calculates the distribution of vertex degrees
	 * @return Distribution of degrees of vertices
	 */
	public Distribution getVertexDegreeDistribution() {
		Distribution answer = new Distribution(0, 100, 1);
		for(List<AssemblyEdge> edges:edgesMap.values()) {
			answer.processDatapoint(edges.size());
		}
		return answer;
	}
	
	public List<AssemblyEmbedded> getAllEmbedded(int sequenceIndex) {
		Map<Integer,AssemblyEmbedded> embeddedSequencesMap = new HashMap<Integer,AssemblyEmbedded>();
		LinkedList<Integer> agenda = new LinkedList<Integer>();
		agenda.add(sequenceIndex);
		while (agenda.size()>0) {
			
			int nextSequenceIdx = agenda.removeFirst();
			List<AssemblyEmbedded> embeddedList = getEmbeddedByHostId(nextSequenceIdx);
			for(AssemblyEmbedded embedded:embeddedList) {
				int seqId = embedded.getSequenceId();
				if(embeddedSequencesMap.containsKey(seqId)) {
					System.err.println("Found two embedded relationships for sequence "+seqId+" parents: "+embedded.getHostId()+" and "+embeddedSequencesMap.get(seqId).getHostId());
					continue;
				}
				AssemblyEmbedded parentObject = embeddedSequencesMap.get(embedded.getHostId());
				if(parentObject==null) {
					embeddedSequencesMap.put(seqId, embedded);
				} else {
					int rootStartParent = parentObject.getHostStart();
					int rootStartSequence = rootStartParent+embedded.getHostStart();
					int rootEndSequence = rootStartParent+embedded.getHostEnd();
					
					boolean reverse = parentObject.isReverse()!=embedded.isReverse();
					embeddedSequencesMap.put(seqId, new AssemblyEmbedded(seqId, embedded.getRead(), reverse, sequenceIndex, rootStartSequence, rootEndSequence));
				}
				
				agenda.add(embedded.getSequenceId());
			}
					
		}
		List<AssemblyEmbedded> answer = new ArrayList<AssemblyEmbedded>();
		answer.addAll(embeddedSequencesMap.values());
		Collections.sort(answer, (a1,a2)-> a1.getHostStart()-a2.getHostStart());
		
		return answer;
	}
	
	public int [] estimateNStatisticsFromPaths () {
		List<Integer> lengths = new ArrayList<Integer>(paths.size());
		for(List<AssemblyEdge> path:paths) {
			int n = path.size();
			if(n==0) continue;
			
			//Ignore small segment at the end
			int lastLength = path.get(0).getOverlap();
			int total = lastLength;
			for (int i=1;i<n;i++) {
				AssemblyEdge edge = path.get(i);
				if(edge.isSameSequenceEdge()) {
					lastLength = edge.getVertex1().getRead().getLength();
					continue;
				}
				int l1 = edge.getVertex1().getRead().getLength();
				int l2 = edge.getVertex2().getRead().getLength();
				if(l1 == lastLength) total += l2 - edge.getOverlap();
				else total += l1 - edge.getOverlap();
			}
			lengths.add(total);
		}
		return NStatisticsCalculator.calculateNStatistics (lengths); 
	}

	public void removeVerticesChimericReads () {
		for(int i=0;i<sequences.size();i++) {
			if(isChimeric(i)) {
				removeVertices(i);
				removeEmbeddedRelations(i);
			}
		}
	}
	
	private boolean isChimeric(int sequenceId) {
		int idxDebug = -1;
		int seqLength = getSequenceLength(sequenceId);
		
		List<AssemblyEmbedded> embeddedList = new ArrayList<AssemblyEmbedded>();
		List<AssemblyEmbedded> emb = embeddedMapByHost.get(sequenceId);
		if(emb==null) emb = new ArrayList<AssemblyEmbedded>();
		if(sequenceId==idxDebug) System.out.println("Finding chimeras. Embedded sequences "+emb.size());
		embeddedList.addAll(emb);
		Collections.sort(embeddedList,(e1,e2)->e1.getHostEvidenceStart()-e2.getHostEvidenceStart());
		List<Integer> hostEvidenceEndsLeft = new ArrayList<Integer>();
		int hostPredictedEndLeft = 0;
		int hostPredictedStartRight = seqLength;
		List<Integer> hostEvidenceStartsRight = new ArrayList<Integer>();
		
		for(AssemblyEmbedded embedded:embeddedList) {
			
			int nextEvidenceStart = embedded.getHostEvidenceStart();
			int nextEvidenceEnd = embedded.getHostEvidenceEnd();
			int unknownLeft = nextEvidenceStart - embedded.getHostStart();
			int unknownRight = embedded.getHostEnd() - nextEvidenceEnd;
			if(sequenceId==idxDebug) System.out.println("Finding chimeras. Embedded "+embedded.getSequenceId()+" reverse"+embedded.isReverse()+" limits: "+nextEvidenceStart+" "+nextEvidenceEnd+" unknown: "+unknownLeft+" "+unknownRight+" count: "+embedded.getNumSharedKmers()+" CSK: "+embedded.getCoverageSharedKmers());
			
			
			if(unknownRight>1000 && unknownLeft<1000) {
				hostEvidenceEndsLeft.add(nextEvidenceEnd);
				hostPredictedEndLeft = Math.max(hostPredictedEndLeft, embedded.getHostEnd());
			}
			if(unknownLeft>1000 && unknownRight<1000) {
				hostEvidenceStartsRight.add(nextEvidenceStart);
				hostPredictedStartRight = Math.min(hostPredictedStartRight, embedded.getHostStart());
			}
		}
		AssemblyVertex vS = verticesStart.get(sequenceId);
		List<AssemblyEdge> edgesS = getEdges(vS);
		for(AssemblyEdge edge: edgesS) {
			if(edge.isSameSequenceEdge()) continue;
			int nextEvidenceEnd = (edge.getVertex1()==vS)?edge.getVertex1EvidenceEnd():edge.getVertex2EvidenceEnd();
			int unknownLeft = (edge.getVertex1()==vS)?edge.getVertex1EvidenceStart():edge.getVertex2EvidenceStart();
			int unknownRight = edge.getOverlap() - nextEvidenceEnd;
			if(sequenceId==idxDebug) System.out.println("Finding chimeras. Edge start "+edge+" evidence end: "+nextEvidenceEnd+" unknown: "+unknownLeft+" "+unknownRight+" count: "+edge.getNumSharedKmers()+" CSK: "+edge.getCoverageSharedKmers());
			if(unknownRight>1000 && unknownLeft<1000) {
				hostEvidenceEndsLeft.add(nextEvidenceEnd);
				hostPredictedEndLeft = Math.max(hostPredictedEndLeft, edge.getOverlap());
			}
		}
		AssemblyVertex vE = verticesEnd.get(sequenceId);
		List<AssemblyEdge> edgesE = getEdges(vE);
		for(AssemblyEdge edge: edgesE) {
			if(edge.isSameSequenceEdge()) continue;
			int nextEvidenceStart = (edge.getVertex1()==vE)?edge.getVertex1EvidenceStart():edge.getVertex2EvidenceStart();
			int unknownLeft = edge.getOverlap() - (seqLength-nextEvidenceStart);
			int unknownRight = seqLength - ((edge.getVertex1()==vE)?edge.getVertex1EvidenceEnd():edge.getVertex2EvidenceEnd());
			if(sequenceId==idxDebug) System.out.println("Finding chimeras. Edge end "+edge+" evidence start: "+nextEvidenceStart+" unknown: "+unknownLeft+" "+unknownRight+" count: "+edge.getNumSharedKmers()+" CSK: "+edge.getCoverageSharedKmers());
			if(unknownLeft>1000 && unknownRight<1000) {
				hostEvidenceStartsRight.add(nextEvidenceStart);
				hostPredictedStartRight = Math.min(hostPredictedStartRight, seqLength-edge.getOverlap());
			}
		}
		
		if(sequenceId==idxDebug) System.out.println("Finding chimeras. Sequence "+sequenceId+". length "+seqLength+" num unknown: "+hostEvidenceEndsLeft.size()+" "+hostEvidenceStartsRight.size());
		if(hostEvidenceStartsRight.size()<2 || hostEvidenceEndsLeft.size()<2) return false;
		Collections.sort(hostEvidenceEndsLeft,(n1,n2)->n1-n2);
		int hostEvidenceEndLeft = hostEvidenceEndsLeft.get(hostEvidenceEndsLeft.size()/2);
		Collections.sort(hostEvidenceStartsRight,(n1,n2)->n1-n2);
		int hostEvidenceStartRight = hostEvidenceStartsRight.get(hostEvidenceStartsRight.size()/2);
		
		int numCrossing = 0;
		for(AssemblyEmbedded embedded:embeddedList) {
			int unknownLeft = embedded.getHostEvidenceStart() - embedded.getHostStart();
			int unknownRight = embedded.getHostEnd() - embedded.getHostEvidenceEnd();
			if(unknownLeft<200 && unknownRight<200 && hostEvidenceEndLeft-embedded.getHostEvidenceStart()>100 && embedded.getHostEvidenceEnd()-hostEvidenceStartRight>100) {
				if(sequenceId==idxDebug) System.out.println("Embedded crossing. Sequence: "+embedded.getSequenceId());
				numCrossing++;
			}
		}
		for(AssemblyEdge edge: edgesS) {
			if(edge.isSameSequenceEdge()) continue;
			int nextEvidenceEnd = (edge.getVertex1()==vS)?edge.getVertex1EvidenceEnd():edge.getVertex2EvidenceEnd();
			int unknownLeft = (edge.getVertex1()==vS)?edge.getVertex1EvidenceStart():edge.getVertex2EvidenceStart();
			int unknownRight = edge.getOverlap() - nextEvidenceEnd;
			if(unknownLeft<200 && unknownRight<200 && nextEvidenceEnd-hostEvidenceStartRight>100) {
				if(sequenceId==idxDebug) System.out.println("Edge start crossing. Edge: "+edge);
				numCrossing++;
			}
		}
		for(AssemblyEdge edge: edgesE) {
			if(edge.isSameSequenceEdge()) continue;
			int nextEvidenceStart = (edge.getVertex1()==vE)?edge.getVertex1EvidenceStart():edge.getVertex2EvidenceStart();
			int unknownLeft = edge.getOverlap() - (seqLength-nextEvidenceStart);
			int unknownRight = seqLength - ((edge.getVertex1()==vE)?edge.getVertex1EvidenceEnd():edge.getVertex2EvidenceEnd());
			if(unknownLeft<200 && unknownRight<200 && hostEvidenceEndLeft - nextEvidenceStart > 100) {
				if(sequenceId==idxDebug) System.out.println("Edge end crossing. Edge: "+edge);
				numCrossing++;
			}
		}
		if(sequenceId==idxDebug) System.out.println("Finding chimeras. Sequence "+sequenceId+". length "+seqLength+" median evidence: "+hostEvidenceEndLeft+" "+hostEvidenceStartRight+" predicted: "+hostPredictedEndLeft+" "+hostPredictedStartRight+" num crossing: "+numCrossing);
		
		if( numCrossing<2 && hostEvidenceEndLeft - hostEvidenceStartRight < 50 && hostEvidenceStartRight - hostEvidenceEndLeft < 200 && hostEvidenceEndLeft-hostPredictedStartRight>100 && hostPredictedEndLeft-hostEvidenceStartRight>100) {
			System.out.println("Possible chimera identified for sequence "+sequenceId+". length "+seqLength+" num unknown: "+hostEvidenceEndsLeft.size()+" "+hostEvidenceStartsRight.size()+" evidence end : "+hostEvidenceEndLeft+" "+hostEvidenceStartRight+" predicted: "+hostPredictedEndLeft+" "+hostPredictedStartRight);
			return true;
		}
		
		return false;
	}
	
	public void filterEdgesAndEmbedded() {
		removeVerticesChimericReads();
		Distribution lengthsDistribution = new Distribution(0, getSequenceLength(0), 1);
		for(QualifiedSequence seq:sequences) lengthsDistribution.processDatapoint(seq.getLength());
		/*Distribution evidenceProportionEmbedded = new Distribution(0, 1, 0.01);
		Distribution cskProportionSelfEmbedded = new Distribution(0, 1, 0.01);
		Distribution wcskProportionSelfEmbedded = new Distribution(0, 1, 0.01);
		for(int seqId:embeddedMapBySequence.keySet()) {
			int length = getSequenceLength(seqId);
			AssemblyEdge edge = getSameSequenceEdge(seqId);
			if(edge == null) continue;
			int selfCSK = edge.getCoverageSharedKmers();
			int selfWCSK = edge.getWeightedCoverageSharedKmers();
			List<AssemblyEmbedded> relations = embeddedMapBySequence.get(seqId);
			for(AssemblyEmbedded embedded:relations) {
				double evidenceLength = embedded.getHostEvidenceEnd()-embedded.getHostEvidenceStart();
				evidenceProportionEmbedded.processDatapoint(evidenceLength/length);
				cskProportionSelfEmbedded.processDatapoint((double)embedded.getCoverageSharedKmers()/selfCSK);
				wcskProportionSelfEmbedded.processDatapoint((double)embedded.getCoverageSharedKmers()/selfWCSK);
			}
		}
		System.out.println("Proportion of evidence vs read length for embedded relationships");
		evidenceProportionEmbedded.printDistribution(System.out);
		System.out.println("Proportion of CSK vs self CSK for embedded relationships");
		cskProportionSelfEmbedded.printDistribution(System.out);
		System.out.println("Proportion of WCSK vs self WCSK for embedded relationships");
		wcskProportionSelfEmbedded.printDistribution(System.out);*/
		int medianLength = getMedianLength();
		System.out.println("Median read length: "+medianLength);
		int numEmbedded = 0;
		for (int seqId = sequences.size()-1; seqId >=0; seqId--) {
			if(filterEdgesAndEmbedded(seqId,medianLength, lengthsDistribution)) numEmbedded++;
		}
		System.out.println("Filtered edges and embedded. Final number of embedded sequences: "+numEmbedded);
		pruneEmbeddedSequences();
		System.out.println("Prunned embedded sequences");
		//filterEdgesCloseRelationships();
	}

	public boolean filterEdgesAndEmbedded(int sequenceId,int medianLength, Distribution lengthsDistribution) {
		int debugIdx = -1;
		int sequenceLength = getSequenceLength(sequenceId);
		AssemblyVertex vS = verticesStart.get(sequenceId);
		AssemblyVertex vE = verticesEnd.get(sequenceId);
		if(vS==null || vE==null) return false;
		filterEdgesAbnormalFeatures(getEdges(vS));
		filterEdgesAbnormalFeatures(getEdges(vE));
		if(sequenceId == debugIdx) System.out.println("Filtered edges with abnormal features");
		List<AssemblyEdge> edgesS = new ArrayList<AssemblyEdge>();
		if(vS!=null) edgesS.addAll(getEdges(vS));
		double minScoreProportionEdges = 0.3;
		//double minScoreProportionEdges = 0;
		
		double maxScoreSE = 0;
		double maxScoreSF = 0;
		for(AssemblyEdge edge: edgesS) {
			if(edge.isSameSequenceEdge()) continue;
			maxScoreSF = Math.max(maxScoreSF, calculateScoreForEdgeFiltering(edge));
			int connectingLength = getSequenceLength(edge.getConnectingVertex(vS).getSequenceIndex());
			/*if(connectingLength<1.2*sequenceLength && edge.getOverlap() >0.8*sequenceLength)*/ maxScoreSE = Math.max(maxScoreSE, calculateScoreForEmbedded(edge));
		}
		List<AssemblyEdge> edgesE = new ArrayList<AssemblyEdge>();
		if(vE!=null) edgesE.addAll(getEdges(vE));
		double maxScoreEE = 0;
		double maxScoreEF = 0;			
		for(AssemblyEdge edge: edgesE) {
			if(edge.isSameSequenceEdge()) continue;
			maxScoreEF = Math.max(maxScoreEF, calculateScoreForEdgeFiltering(edge));
			int connectingLength = getSequenceLength(edge.getConnectingVertex(vE).getSequenceIndex());
			/*if(connectingLength<1.5*sequenceLength && edge.getOverlap() >0.8*sequenceLength)*/ maxScoreEE = Math.max(maxScoreEE, calculateScoreForEmbedded(edge));
		}
		double minScoreEdges = minScoreProportionEdges*Math.max(maxScoreSF, maxScoreEF);
		for(AssemblyEdge edge: edgesS) {
			if(edge.isSameSequenceEdge()) continue;
			double score = calculateScoreForEdgeFiltering(edge);
			if(sequenceId == debugIdx) System.out.println("Assembly graph. Next edge start "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber()+" overlap: "+edge.getOverlap()+" score: "+score+" Max score start: "+maxScoreSF+" limit: "+minScoreEdges);
			if(score <500 || (score < maxScoreSF && score < minScoreEdges)) {
				if(sequenceId == debugIdx) System.out.println("Assembly graph. Removing edge: "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber());
				removeEdge(edge);
			}
		}
		if(sequenceId == debugIdx) System.out.println("Assembly graph. Initial edges start "+edgesS.size()+" Max scores start: "+maxScoreSF+" "+maxScoreSE+" remaining edges: "+edgesMap.get(vS.getUniqueNumber()).size());
		
		for(AssemblyEdge edge: edgesE) {
			if(edge.isSameSequenceEdge()) continue;
			double score = calculateScoreForEdgeFiltering(edge);
			if(sequenceId == debugIdx) System.out.println("Assembly graph. Next edge end "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber()+" overlap: "+edge.getOverlap()+" score: "+score+" Max score end: "+maxScoreEF+" limit: "+minScoreEdges);
			if(score < 500 || (score < maxScoreEF && score < minScoreEdges)) {
				if(sequenceId == debugIdx) System.out.println("Assembly graph. Removing edge: "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber());
				removeEdge(edge);
			}
		}
		if(sequenceId == debugIdx) System.out.println("Assembly graph. Initial edges end "+edgesE.size()+" Max scores end: "+maxScoreEF+" "+maxScoreEE+" remaining edges: "+edgesMap.get(vE.getUniqueNumber()).size());
		
		double medianRelationship = 1.0*sequenceLength/(double)medianLength;
		
		//double minScoreProportionEmbedded = 0.8;
		//double cumulative = lengthsDistribution.getCumulativeCount(sequenceLength)/lengthsDistribution.getCount();
		double minScoreProportionEmbedded = Math.min(0.8, 0.5*medianRelationship);
		//double minScoreProportionEmbedded = 0.8*cumulative;
		if(minScoreProportionEmbedded<0.5) minScoreProportionEmbedded = 0.5;
		//if(medianRelationship>1 && minScoreProportionEmbedded<0.7) minScoreProportionEmbedded = 0.7;
		
		double maxScoreFilterEmbedded = minScoreProportionEmbedded*Math.max(maxScoreSE, maxScoreEE);
	
		AssemblyEdge sameSequenceEdge = getSameSequenceEdge(sequenceId);
		double sameSeqCSK = sameSequenceEdge.getCoverageSharedKmers();
		double sameSeqWCSK = sameSequenceEdge.getCoverageSharedKmers();
		List<AssemblyEmbedded> embeddedList= new ArrayList<AssemblyEmbedded>();
		embeddedList.addAll(getEmbeddedBySequenceId(sequenceId));
		if(embeddedList.size()==0) return false;
		double maxEvidencePropEmbedded = 0;
		AssemblyEmbedded embeddedMax = null;
		double maxScoreEmbedded = -1;
		int countPass = 0;
		for(AssemblyEmbedded embedded:embeddedList) {
			maxEvidencePropEmbedded = Math.max(maxEvidencePropEmbedded, embedded.calculateEvidenceProportion());
			double CSKprop = (double)embedded.getCoverageSharedKmers()/sameSeqCSK;
			double WCSKprop = (double)embedded.getWeightedCoverageSharedKmers()/sameSeqWCSK;
			double score = calculateScore(embedded);
			if(sequenceId == debugIdx) System.out.println("Assembly graph. Next embedded "+embedded+" score: "+score+" evProp: "+embedded.calculateEvidenceProportion()+" CSK prop "+CSKprop+" WCSK prop: "+WCSKprop+" Indels: "+embedded.getNumIndels());
			
			if(embeddedMax==null || maxScoreEmbedded<score) {
				maxScoreEmbedded = score;
				embeddedMax = embedded;
			}
			
			//if(evidenceProp*CSKprop >=0.25) countPass++;
			if(embedded.calculateEvidenceProportion() >=0.95) countPass++;
		}
			
		//Score proportion filter calculation
		double maxScorePropEmbedded = maxScoreEmbedded/calculateScoreForEmbedded(sameSequenceEdge);
		if(sequenceId == debugIdx) System.out.println("Assembly graph. Median relationship: "+medianRelationship+" evidence proportion "+ maxEvidencePropEmbedded +" max score embedded "+maxScoreEmbedded+" same seq score: "+calculateScoreForEmbedded(sameSequenceEdge)+ " max score prop self: "+maxScorePropEmbedded+" count pass: "+countPass);
		
		//if(countPass==0) {
		if(maxScoreEmbedded<maxScoreFilterEmbedded) {
		//if(maxEvidencePropEmbedded<evidenceProportionThreshold) {
		//if(maxEvidencePropEmbedded<minProportionEmbedded || maxScorePropEmbedded<0.5*minProportionEmbedded) {
		//if(maxScoreEmbedded<maxScoreFilterEmbedded || maxScoreEmbedded < 0.2*calculateScoreForEmbedded(sameSequenceEdge)) {
		//if(maxScoreEmbedded < 0.2*calculateScoreForEmbedded(sameSequenceEdge)) {
			//Replace embedded relationships with edges to make the sequence not embedded
			for(AssemblyEmbedded embedded:embeddedList) {
				removeEmbedded(embedded);
				if(sequenceId == debugIdx) System.out.println("Adding edge replacing embedded "+embedded.getHostId()+" limits: "+embedded.getHostStart()+" "+embedded.getHostEnd()+" host length: "+getSequenceLength(embedded.getHostId())+"score: "+calculateScore(embedded));
				addEdgeFromEmbedded(embedded);
			}
			return false;
		} else {
			if(sequenceId == debugIdx) System.out.println("Sequence is embedded ");
			for(AssemblyEmbedded embedded:embeddedList) {
				if(embedded!=embeddedMax) {
					removeEmbedded(embedded);
					//if(sequenceId == debugIdx) System.out.println("Assembly graph. Removed embedded host: "+embedded.getHostId()+" Embedded relations: "+embeddedMapBySequence.get(sequenceId)+" is embedded: "+isEmbedded(sequenceId));
				}
			}
			return true;
		}
	}

	private void addEdgeFromEmbedded(AssemblyEmbedded embedded) {
		int distanceStart = embedded.getHostStart();
		int distanceEnd = getSequenceLength(embedded.getHostId())-embedded.getHostEnd();
		AssemblyVertex vertexHost=null;
		AssemblyVertex vertexEmbedded=null;
		if(distanceStart<0.5*distanceEnd) {
			vertexHost = getVertex(embedded.getHostId(), true);
			vertexEmbedded = getVertex(embedded.getSequenceId(), embedded.isReverse());
		} else if (distanceEnd<0.5*distanceStart) {
			vertexHost = getVertex(embedded.getHostId(), false);
			vertexEmbedded = getVertex(embedded.getSequenceId(), !embedded.isReverse());
		}
		if(vertexHost==null || vertexEmbedded==null) return;
		int overlap = Math.max(embedded.getCoverageSharedKmers(), embedded.getHostEvidenceEnd()-embedded.getHostEvidenceStart());
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
		addEdge(edge);
	}

	private void filterEdgesAbnormalFeatures(List<AssemblyEdge> edges) {
		List<AssemblyEdge> toRemove = new ArrayList<AssemblyEdge>();
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge()) continue;
			int overlap = edge.getOverlap();
			int v1L = edge.getVertex1().getRead().getLength();
			int v2L = edge.getVertex2().getRead().getLength();
			if(overlap>1.1*v1L || overlap>1.1*v2L) {
				toRemove.add(edge);
			}/* else if (edge.getWeightedCoverageSharedKmers()<0.1*v1L || edge.getWeightedCoverageSharedKmers()<0.1*v2L) {
				toRemove.add(edge);
			}*/
		}
		for(AssemblyEdge edge:toRemove) removeEdge(edge);
	}

	private double calculateScore(AssemblyEmbedded embedded) {
		double evProp = embedded.calculateEvidenceProportion();
		//return embedded.getCoverageSharedKmers();
		//return embedded.getWeightedCoverageSharedKmers();
		//return embedded.getRawKmerHits();
		return embedded.getWeightedCoverageSharedKmers()*evProp*evProp;
		//return embedded.getWeightedCoverageSharedKmers()*embedded.getRawKmerHits();
	}

	private double calculateScoreForEmbedded(AssemblyEdge edge) {
		double evProp = edge.calculateEvidenceProportion();
		//return edge.getCoverageSharedKmers();
		//return edge.getWeightedCoverageSharedKmers();
		//return edge.getRawKmerHits();
		return edge.getWeightedCoverageSharedKmers()*evProp*evProp;
		//return edge.getWeightedCoverageSharedKmers()*edge.getRawKmerHits();
	}
	
	private double calculateScoreForEdgeFiltering(AssemblyEdge edge) {
		double evProp = edge.calculateEvidenceProportion();
		//return edge.getCoverageSharedKmers();
		return edge.getWeightedCoverageSharedKmers()*evProp*evProp;
		//return edge.getRawKmerHits();
		//return edge.getWeightedCoverageSharedKmers()*edge.getRawKmerHits();
	}
}
