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
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import JSci.maths.statistics.NormalDistribution;
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
	 * Private constructor for subgraphs
	 */
	private AssemblyGraph () {
		
	}
	public AssemblyGraph(List<QualifiedSequence> sequences) {
		int n = sequences.size();
		this.sequences = Collections.unmodifiableList(sequences);
		cumulativeReadLength = new long [n];
		initStructures(n);
		
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
	private void initStructures (int n) {
		verticesStart = new HashMap<>(n);
		verticesEnd = new HashMap<>(n);
		verticesByUnique = new HashMap<>(n);
		edgesMap = new HashMap<>(n);
	}
	public AssemblyGraph buildSubgraph(Set<Integer> readIdsCluster) {
		AssemblyGraph subgraph = new AssemblyGraph();
		int n = sequences.size();
		subgraph.sequences = sequences;
		subgraph.cumulativeReadLength = cumulativeReadLength;
		subgraph.initStructures(n);
		//Add vertices
		for(AssemblyVertex vertex:verticesByUnique.values()) {
			if(readIdsCluster==null || readIdsCluster.contains(vertex.getSequenceIndex())) {
				subgraph.verticesByUnique.put(vertex.getUniqueNumber(), vertex);
				subgraph.edgesMap.put(vertex.getUniqueNumber(), new ArrayList<>());
				if(vertex.isStart()) subgraph.verticesStart.put(vertex.getSequenceIndex(),vertex);
				else subgraph.verticesEnd.put(vertex.getSequenceIndex(),vertex);
			}
		}
		//Add edges within the subgraph
		List<AssemblyEdge> edges = getEdges();
		for(AssemblyEdge edge:edges) {
			if(readIdsCluster == null || (readIdsCluster.contains(edge.getVertex1().getSequenceIndex()) && readIdsCluster.contains(edge.getVertex2().getSequenceIndex()))) {
				subgraph.addEdge(edge);
			}
		}
		//Add embedded relationships
		for(List<AssemblyEmbedded> embeddedList:embeddedMapBySequence.values()) {
			for(AssemblyEmbedded embedded:embeddedList) {
				if(readIdsCluster == null || (readIdsCluster.contains(embedded.getSequenceId()) && readIdsCluster.contains(embedded.getHostId()))) {
					subgraph.addEmbedded(embedded);
				}
			}
			
		}
		return subgraph;
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
	
	public long [] estimateNStatisticsFromPaths () {
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
		int numIncompleteEdgesLeft = 0;
		AssemblyVertex vS = verticesStart.get(sequenceId);
		List<AssemblyEdge> edgesS = getEdges(vS);
		AssemblyVertex vE = verticesEnd.get(sequenceId);
		List<AssemblyEdge> edgesE = getEdges(vE);
		int numIncompleteEdgesRight = 0;
		if(hostEvidenceEndsLeft.size()<5 || hostEvidenceStartsRight.size()<5) {
			
			for(AssemblyEdge edge: edgesS) {
				if(edge.isSameSequenceEdge()) continue;
				int nextEvidenceStart = (edge.getVertex1()==vS)?edge.getVertex1EvidenceStart():edge.getVertex2EvidenceStart();
				int nextEvidenceEnd = (edge.getVertex1()==vS)?edge.getVertex1EvidenceEnd():edge.getVertex2EvidenceEnd();
				int unknownLeft = nextEvidenceStart;
				int unknownRight = edge.getOverlap() - nextEvidenceEnd;
				if(unknownRight<0) continue;
				if(sequenceId==idxDebug) System.out.println("Finding chimeras. Edge start "+edge+" evidence end: "+nextEvidenceEnd+" unknown: "+unknownLeft+" "+unknownRight+" count: "+edge.getNumSharedKmers()+" CSK: "+edge.getCoverageSharedKmers());
				if(unknownRight>1000 && unknownLeft<1000) {
					hostEvidenceEndsLeft.add(nextEvidenceEnd);
					hostPredictedEndLeft = Math.max(hostPredictedEndLeft, edge.getOverlap());
				}
				if(unknownLeft>1000 && unknownRight<1000) {
					hostEvidenceStartsRight.add(nextEvidenceStart);
					numIncompleteEdgesLeft++;
					//hostPredictedStartRight = 0;
				}
			}
			for(AssemblyEdge edge: edgesE) {
				if(edge.isSameSequenceEdge()) continue;
				int nextEvidenceStart = (edge.getVertex1()==vE)?edge.getVertex1EvidenceStart():edge.getVertex2EvidenceStart();
				int nextEvidenceEnd = (edge.getVertex1()==vE)?edge.getVertex1EvidenceEnd():edge.getVertex2EvidenceEnd();
				int unknownLeft = edge.getOverlap() - (seqLength-nextEvidenceStart);
				if(unknownLeft<0) continue;
				int unknownRight = seqLength - nextEvidenceEnd;
				if(sequenceId==idxDebug) System.out.println("Finding chimeras. Edge end "+edge+" evidence start: "+nextEvidenceStart+" unknown: "+unknownLeft+" "+unknownRight+" count: "+edge.getNumSharedKmers()+" CSK: "+edge.getCoverageSharedKmers());
				if(unknownRight>1000 && unknownLeft<1000) {
					hostEvidenceEndsLeft.add(nextEvidenceEnd);
					numIncompleteEdgesRight++;
					//hostPredictedEndLeft = seqLength;
				}
				if(unknownLeft>1000 && unknownRight<1000) {
					hostEvidenceStartsRight.add(nextEvidenceStart);
					hostPredictedStartRight = Math.min(hostPredictedStartRight, seqLength-edge.getOverlap());
				}
			}
		}
		
		if(sequenceId==idxDebug) System.out.println("Finding chimeras. Sequence "+sequenceId+". length "+seqLength+" num unknown: "+hostEvidenceEndsLeft.size()+" "+hostEvidenceStartsRight.size()+" starts right: "+hostEvidenceStartsRight);
		if(hostEvidenceStartsRight.size()==0 || hostEvidenceEndsLeft.size()==0 || hostEvidenceStartsRight.size()+numIncompleteEdgesRight<5 || hostEvidenceEndsLeft.size()+numIncompleteEdgesLeft<5 ) return false;
		//if(numIncompleteEdgesLeft<3 && numIncompleteEdgesRight<3 && (hostEvidenceStartsRight.size()<3 || hostEvidenceEndsLeft.size()<3)) return false;
		Collections.sort(hostEvidenceEndsLeft,(n1,n2)->n1-n2);
		int hostEvidenceEndLeft = hostEvidenceEndsLeft.size()>0?hostEvidenceEndsLeft.get(hostEvidenceEndsLeft.size()/2):0;
		Collections.sort(hostEvidenceStartsRight,(n1,n2)->n1-n2);
		int hostEvidenceStartRight = hostEvidenceStartsRight.size()>0?hostEvidenceStartsRight.get(hostEvidenceStartsRight.size()/2):seqLength;
		int minEvidenceStart = Math.min(hostEvidenceStartRight, hostEvidenceEndLeft);
		int maxEvidenceEnd = Math.max(hostEvidenceStartRight, hostEvidenceEndLeft);
		
		int numCrossing = 0;
		for(AssemblyEmbedded embedded:embeddedList) {
			int unknownLeft = embedded.getHostEvidenceStart() - embedded.getHostStart();
			int unknownRight = embedded.getHostEnd() - embedded.getHostEvidenceEnd();
			if(unknownLeft<200 && unknownRight<200 && minEvidenceStart-embedded.getHostEvidenceStart()>100 && embedded.getHostEvidenceEnd()-maxEvidenceEnd>100) {
				if(sequenceId==idxDebug) System.out.println("Embedded crossing. Sequence: "+embedded.getSequenceId());
				numCrossing++;
			}
		}
		for(AssemblyEdge edge: edgesS) {
			if(edge.isSameSequenceEdge()) continue;
			int nextEvidenceEnd = (edge.getVertex1()==vS)?edge.getVertex1EvidenceEnd():edge.getVertex2EvidenceEnd();
			int unknownLeft = (edge.getVertex1()==vS)?edge.getVertex1EvidenceStart():edge.getVertex2EvidenceStart();
			int unknownRight = edge.getOverlap() - nextEvidenceEnd;
			
			if(unknownLeft<200 && unknownRight<200 && nextEvidenceEnd-maxEvidenceEnd>100) {
				if(sequenceId==idxDebug) System.out.println("Edge start crossing. Edge: "+edge);
				numCrossing++;
			}
		}
		for(AssemblyEdge edge: edgesE) {
			if(edge.isSameSequenceEdge()) continue;
			int nextEvidenceStart = (edge.getVertex1()==vE)?edge.getVertex1EvidenceStart():edge.getVertex2EvidenceStart();
			int unknownLeft = edge.getOverlap() - (seqLength-nextEvidenceStart);
			int unknownRight = seqLength - ((edge.getVertex1()==vE)?edge.getVertex1EvidenceEnd():edge.getVertex2EvidenceEnd());
			
			if(unknownLeft<200 && unknownRight<200 && minEvidenceStart - nextEvidenceStart > 100) {
				if(sequenceId==idxDebug) System.out.println("Edge end crossing. Edge: "+edge);
				numCrossing++;
			}
		}
		if(sequenceId==idxDebug) System.out.println("Finding chimeras. Sequence "+sequenceId+". length "+seqLength+" median evidence: "+hostEvidenceEndLeft+" "+hostEvidenceStartRight+" predicted: "+hostPredictedEndLeft+" "+hostPredictedStartRight+" num crossing: "+numCrossing+" incomplete: "+numIncompleteEdgesLeft+" "+numIncompleteEdgesRight);
		int d1 = hostEvidenceEndLeft-hostPredictedStartRight;
		int d2 = hostPredictedEndLeft-hostEvidenceStartRight;
		int d3 = hostPredictedEndLeft-hostEvidenceEndLeft;
		int d4 = hostEvidenceStartRight-hostPredictedStartRight;
		int d5 = hostEvidenceEndLeft - hostEvidenceStartRight;
		if( numCrossing==0  && d1>1000 && d2>1000 && d3>2000 && d4>2000 && d5<10000) {
			System.out.println("Possible chimera identified for sequence "+sequenceId+". length "+seqLength+" num unknown: "+hostEvidenceEndsLeft.size()+" "+hostEvidenceStartsRight.size()+" evidence end : "+hostEvidenceEndLeft+" "+hostEvidenceStartRight+" predicted: "+hostPredictedEndLeft+" "+hostPredictedStartRight+" crossing: "+numCrossing);
			return true;
		}/* else if (numCrossing==0 && hostEvidenceEndLeft==0 && hostEvidenceStartRight>1000 && numIncompleteEdgesLeft>5) {
			System.out.println("Possible chimera identified for start of sequence "+sequenceId+". length "+seqLength+" num unknown: "+hostEvidenceEndsLeft.size()+" "+hostEvidenceStartsRight.size()+" evidence end : "+hostEvidenceEndLeft+" "+hostEvidenceStartRight+" predicted: "+hostPredictedEndLeft+" "+hostPredictedStartRight+" num incomplete: "+numIncompleteEdgesLeft+" "+numIncompleteEdgesRight);
			return true;
		} else if (numCrossing==0 && hostEvidenceEndLeft<seqLength-1000 && hostEvidenceStartRight==seqLength && numIncompleteEdgesRight>5) {
			System.out.println("Possible chimera identified for end of sequence "+sequenceId+". length "+seqLength+" num unknown: "+hostEvidenceEndsLeft.size()+" "+hostEvidenceStartsRight.size()+" evidence end : "+hostEvidenceEndLeft+" "+hostEvidenceStartRight+" predicted: "+hostPredictedEndLeft+" "+hostPredictedStartRight+" num incomplete: "+numIncompleteEdgesLeft+" "+numIncompleteEdgesRight);
			return true;
		}*/
		
		return false;
	}
	
	public List<AssemblyEdge> selectSafeEdges(  ) {
		Set<Integer> repetitiveVertices = predictRepetitiveVertices();
		return selectSafeEdges(repetitiveVertices);
	}
	
	private List<AssemblyEdge> selectSafeEdges( Set<Integer> repetitiveVertices ) {
		List<AssemblyEdge> allEdges = getEdges();
		List<AssemblyEdge> safeEdges = new ArrayList<AssemblyEdge>();
		for(AssemblyEdge edge:allEdges) {
			if(edge.isSameSequenceEdge()) continue;
			if (edge.getEvidenceProportion()<0.9) continue;
			if (edge.getIndelsPerKbp()>20) continue;
			boolean r1 = repetitiveVertices.contains(edge.getVertex1().getUniqueNumber());
			boolean r2 = repetitiveVertices.contains(edge.getVertex2().getUniqueNumber()); 
			int d1 = getEdges(edge.getVertex1()).size();
			int d2 = getEdges(edge.getVertex2()).size();
			//if(logEdge(edge)) System.out.println("Select safe edges. repetitive: "+r1+" "+r2+" initial degrees "+edge.getVertex1().getDegreeUnfilteredGraph()+" "+edge.getVertex2().getDegreeUnfilteredGraph()+" current degrees "+d1+" "+d2+" reciprocal best: "+isRecipocalBest(graph, edge)+" edge: "+edge);
			if((d1==2 && d2==2) || (!r1 && !r2 && isRecipocalBest(edge))) {
				safeEdges.add(edge);
			}
		}
		return safeEdges;
	}
	public Set<Integer> predictRepetitiveVertices() {
		Distribution initialDegreesDist = new Distribution(0, 10000, 1);
		for(AssemblyVertex v:verticesByUnique.values()) {
			if(!isEmbedded(v.getSequenceIndex())) initialDegreesDist.processDatapoint(v.getDegreeUnfilteredGraph());
		}
		System.out.println("Degree average: "+initialDegreesDist.getAverage()+" variance "+initialDegreesDist.getVariance());
		NormalDistribution distDegrees = new NormalDistribution(initialDegreesDist.getAverage(), initialDegreesDist.getVariance());
		
		Set<Integer> repetitiveVertices = new HashSet<Integer>();
		for(AssemblyVertex v:verticesByUnique.values()) {
			if(isRepetivive(v, distDegrees)) repetitiveVertices.add(v.getUniqueNumber());
		}
		System.out.println("Total vertices for layout: "+initialDegreesDist.getCount()+" Number of repetitive vertices: "+repetitiveVertices.size());
		return repetitiveVertices;
	}
	private boolean isRepetivive(AssemblyVertex vertex, NormalDistribution distDegrees) {
		int degree = vertex.getDegreeUnfilteredGraph();
		double pValue = distDegrees.cumulative(degree);
		boolean repetitive = pValue>0.999;
		//if(repetitive) System.out.println("Repetitive vertex: "+vertex+" initial degree: "+vertex.getDegreeUnfilteredGraph()+" average: "+distDegrees.getMean()+" sd "+Math.sqrt(distDegrees.getVariance()));
		return repetitive;
	}
	public boolean isRecipocalBest(AssemblyEdge edge) {
		AssemblyVertex v1 = edge.getVertex1();
		if(!isBestEdge(v1, edge)) return false;
		AssemblyVertex v2 = edge.getVertex2();
		if(!isBestEdge(v2, edge)) return false;
		return true;
	}
	public boolean isBestEdge (AssemblyVertex v, AssemblyEdge edgeT) {
		List<AssemblyEdge> edges = getEdges(v);
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge() || edge==edgeT) continue;
			if(edge.getOverlap()>=edgeT.getOverlap()) return false;
			if(edge.getWeightedCoverageSharedKmers()>=edgeT.getWeightedCoverageSharedKmers()) return false;
		}
		return true;
	}
	
	private Distribution[] calculateStatistics(List<AssemblyEdge> safeEdges, Set<Integer> repetitiveVertices) {
		Distribution overlapDistributionTP = new Distribution(0, 100000, 1);
		Distribution kmerHitCoverageDistributionTP = new Distribution(0, 100000, 1);
		Distribution kmerHitWCovDistributionTP = new Distribution(0, 100000, 1);
		Distribution coverageProportionDistributionTP = new Distribution(0, 1.5, 0.01);
		Distribution evidenceProportionDistributionTP = new Distribution(0, 1.1, 0.01);
		Distribution indelsKbpDistributionTP = new Distribution(0, 300, 1);
		
		for(AssemblyEdge edge:safeEdges) {
			if (edge.isSameSequenceEdge()) continue;
			double overlap = edge.getOverlap();
			overlapDistributionTP.processDatapoint(overlap);
			kmerHitCoverageDistributionTP.processDatapoint(edge.getCoverageSharedKmers());
			kmerHitWCovDistributionTP.processDatapoint(edge.getWeightedCoverageSharedKmers());
			coverageProportionDistributionTP.processDatapoint((double)edge.getWeightedCoverageSharedKmers()/overlap);
			evidenceProportionDistributionTP.processDatapoint(edge.getEvidenceProportion());
			indelsKbpDistributionTP.processDatapoint(edge.getIndelsPerKbp());
		}
		Distribution [] answer = {overlapDistributionTP, kmerHitCoverageDistributionTP, kmerHitWCovDistributionTP, coverageProportionDistributionTP,evidenceProportionDistributionTP,indelsKbpDistributionTP};
		return answer;
	}
	
	public void updateScores (boolean useIndels) {
		updateVertexDegrees();
		Set<Integer> repetitiveVertices = predictRepetitiveVertices();
		List<AssemblyEdge> safeEdges = selectSafeEdges(repetitiveVertices);
		System.out.println("Number of safe edges: "+safeEdges.size());
		Distribution [] edgesStats = calculateStatistics(safeEdges, repetitiveVertices);
		System.out.println("Average overlap TP: "+edgesStats[0].getAverage()+" SD: "+Math.sqrt(edgesStats[0].getVariance())+ " Total: "+edgesStats[0].getCount());
		System.out.println("Average coverage shared kmers TP: "+edgesStats[1].getAverage()+" SD: "+Math.sqrt(edgesStats[1].getVariance())+ " Total: "+edgesStats[1].getCount());
		System.out.println("Average weighted coverage shared kmers TP: "+edgesStats[2].getAverage()+" SD: "+Math.sqrt(edgesStats[2].getVariance())+ " Total: "+edgesStats[2].getCount());
		System.out.println("Average weighted coverage proportion TP: "+edgesStats[3].getAverage()+" SD: "+Math.sqrt(edgesStats[3].getVariance())+ " Total: "+edgesStats[3].getCount());
		System.out.println("Average Evidence proportion TP: "+edgesStats[4].getAverage()+" SD: "+Math.sqrt(edgesStats[4].getVariance())+ " Total: "+edgesStats[4].getCount());
		System.out.println("Average indels kbp TP: "+edgesStats[5].getAverage()+" SD: "+Math.sqrt(edgesStats[5].getVariance())+ " Total: "+edgesStats[5].getCount());
		//System.out.println("Distribution indels kbp");
		//edgesStats[5].printDistribution(System.out);
		AssemblySequencesRelationshipScoresCalculator calculator = new AssemblySequencesRelationshipScoresCalculator();
		calculator.setUseIndels(useIndels);
		List<AssemblyEdge> allEdges = getEdges();
		for(AssemblyEdge edge: allEdges) {
			edge.setScore(calculator.calculateScore(edge,edgesStats));
			edge.setCost(calculator.calculateCost(edge,edgesStats));
		}
		for (List<AssemblyEmbedded> embeddedList:embeddedMapBySequence.values()) {
			for(AssemblyEmbedded embedded:embeddedList) {
				embedded.setScore(calculator.calculateScore(embedded, edgesStats));
				embedded.setCost(calculator.calculateCost(embedded, edgesStats));
			}
		}
	}
	
}
