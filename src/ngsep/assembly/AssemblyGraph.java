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
	public static final int DEF_PLOIDY_ASSEMBLY = 1;
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

	private List<AssemblyPath> paths = new ArrayList<AssemblyPath>();
	
	private Set<Integer> chimericSequenceIds = new HashSet<Integer>();
	
	private int numEdges = 0;
	
	private int ploidy = DEF_PLOIDY_ASSEMBLY;
	
	private long expectedAssemblyLength = 0;

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
	
	public void removeEmbeddedRelations(int sequenceId) {
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

	public synchronized void addPath(AssemblyPath path) {
		paths.add(path);
		for(AssemblyEdge edge:path.getEdges()) edge.setLayoutEdge(true);
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
	
	public List<AssemblyEdge> getEdgesBySequenceId(int seqId) {
		List<AssemblyEdge> answer = new ArrayList<AssemblyEdge>();
		AssemblyVertex v1 = getVertex(seqId, true);
		if(v1!=null) {
			List<AssemblyEdge> edgesV1 = getEdges(v1); 
			if(edgesV1!=null) answer.addAll(edgesV1);
		}
		AssemblyVertex v2 = getVertex(seqId, false);
		if(v2!=null) {
			List<AssemblyEdge> edgesV2 = getEdges(v2);
			if(edgesV2!=null) answer.addAll(getEdges(v2));
		}
		return answer;
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
	public List<AssemblyPath> getPaths() {
		return paths;
	}
	
	public int getPloidy() {
		return ploidy;
	}
	public void setPloidy(int ploidy) {
		this.ploidy = ploidy;
	}
	public long getExpectedAssemblyLength() {
		return expectedAssemblyLength;
	}
	public void setExpectedAssemblyLength(long expectedAssemblyLength) {
		this.expectedAssemblyLength = expectedAssemblyLength;
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
	public List<AssemblyEmbedded> getAllEmbedded() {
		List<AssemblyEmbedded> answer = new ArrayList<AssemblyEmbedded>();
		for(List<AssemblyEmbedded> rels:embeddedMapByHost.values()) {
			answer.addAll(rels);
		}
		return answer;
	}
	public List<AssemblyEmbedded> getAllEmbedded(int sequenceIndex) {
		Map<Integer,AssemblyEmbedded> embeddedSequencesMap = new HashMap<Integer,AssemblyEmbedded>();
		QualifiedSequence rootSequence = getSequence(sequenceIndex);
		LinkedList<Integer> agenda = new LinkedList<Integer>();
		agenda.add(sequenceIndex);
		while (agenda.size()>0) {
			
			int nextSequenceIdx = agenda.removeFirst();
			int nextSequenceLength = getSequenceLength(nextSequenceIdx);
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
					if(parentObject.isReverse()) {
						rootStartSequence = rootStartParent+(nextSequenceLength-embedded.getHostEnd());
						rootEndSequence = rootStartParent+(nextSequenceLength-embedded.getHostStart());
					}
					boolean reverse = parentObject.isReverse()!=embedded.isReverse();
					
					embeddedSequencesMap.put(seqId, new AssemblyEmbedded(seqId, embedded.getRead(), reverse, sequenceIndex, rootSequence, rootStartSequence, rootEndSequence));
				}
				
				agenda.add(embedded.getSequenceId());
			}
					
		}
		List<AssemblyEmbedded> answer = new ArrayList<AssemblyEmbedded>();
		answer.addAll(embeddedSequencesMap.values());
		Collections.sort(answer, (a1,a2)-> a1.getHostStart()-a2.getHostStart());
		
		return answer;
	}
	
	public List<AssemblyPath> buildPaths(List<AssemblyEdge> edges) {
		Map<Integer,AssemblyEdge> edgesByVertex = new HashMap<Integer, AssemblyEdge>(); 
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge()) continue;
			AssemblyVertex v1 = edge.getVertex1();
			if(edgesByVertex.containsKey(v1.getUniqueNumber())) System.err.println("WARN: two edges for vertex. E1 "+edgesByVertex.get(v1.getUniqueNumber())+" E2: "+edge);
			AssemblyVertex v2 = edge.getVertex2();
			if(edgesByVertex.containsKey(v2.getUniqueNumber())) System.err.println("WARN: two edges for vertex. E1 "+edgesByVertex.get(v2.getUniqueNumber())+" E2: "+edge);
			edgesByVertex.put(v1.getUniqueNumber(), edge);
			edgesByVertex.put(v2.getUniqueNumber(), edge);
		}
		int n = getNumSequences();
		List<AssemblyPath> paths = new ArrayList<AssemblyPath>();
		Set<Integer> sequencesInPaths = new HashSet<>();
		for(int i=0;i<n;i++) {
			if(getVertex(i, true)==null || getVertex(i, false)==null) continue;
			if(isEmbedded(i)) continue;
			
			if(sequencesInPaths.contains(i)) continue;
			AssemblyVertex nextVertex = getVertex(i, true);
			if(getEdges(nextVertex).size()==0) continue;
			AssemblyPath path = new AssemblyPath(getSameSequenceEdge(nextVertex));
			nextVertex = path.getVertexLeft();
			//Expand v1
			while(nextVertex!=null) {
				AssemblyEdge nextEdge = edgesByVertex.get(nextVertex.getUniqueNumber());
				if(nextEdge == null) nextVertex=null;
				else {
					nextVertex = nextEdge.getConnectingVertex(nextVertex);
					if(sequencesInPaths.contains(nextVertex.getSequenceIndex())) {
						System.err.println("WARN: Cycle detected building paths. Next edge: "+nextEdge);
						break;
					}
					path.connectEdgeLeft(this, nextEdge);
					sequencesInPaths.add(nextVertex.getSequenceIndex());
					nextVertex = path.getVertexLeft();
				}
			}
			//Expand v2
			nextVertex = path.getVertexRight();
			while(nextVertex!=null) {
				AssemblyEdge nextEdge = edgesByVertex.get(nextVertex.getUniqueNumber());
				if(nextEdge == null) nextVertex=null;
				else {
					nextVertex = nextEdge.getConnectingVertex(nextVertex);
					if(sequencesInPaths.contains(nextVertex.getSequenceIndex())) {
						System.err.println("WARN: Cycle detected building paths. Next edge: "+nextEdge);
						break;
					}
					path.connectEdgeRight(this, nextEdge);
					sequencesInPaths.add(nextVertex.getSequenceIndex());
					nextVertex = path.getVertexRight();
				}
			}
			paths.add(path);
		} 
		return paths;
	}
	
	public long [] estimateNStatisticsFromPaths () {
		List<Integer> lengths = new ArrayList<Integer>(paths.size());
		for(AssemblyPath path:paths) {
			List<AssemblyEdge> edges = path.getEdges();
			
			//Ignore small segment at the end
			int lastLength = edges.get(0).getOverlap();
			int total = lastLength;
			for (AssemblyEdge edge:edges) {
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
			if(calculateChimericStatus(i)) {
				removeVertices(i);
				removeEmbeddedRelations(i);
				chimericSequenceIds.add(i);
			}
		}
	}
	
	public Set<Integer> calculateEmbeddedToChimeric() {
		Set<Integer> chimericIds = new HashSet<>();
		for(int i=0;i<sequences.size();i++) {
			if(calculateChimericStatus(i)) chimericIds.add(i);
		}
		Set<Integer> answer = new HashSet<>();
		for(Map.Entry<Integer,List<AssemblyEmbedded>> entry:embeddedMapBySequence.entrySet()) {
			int i = entry.getKey();
			List<AssemblyEmbedded> relsHosts = entry.getValue();
			boolean orphan = true;
			for(AssemblyEmbedded rel:relsHosts) {
				orphan = chimericIds.contains(rel.getHostId());
				if(!orphan) break;
			}
			if(orphan) answer.add(i);
		}
		return answer;
	}
	
	private boolean calculateChimericStatus(int sequenceId) {
		if(verticesStart.get(sequenceId)==null || verticesEnd.get(sequenceId)==null) return false;
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
			if(embedded.getEvidenceProportion()<0.5) continue;
			int nextEvidenceStart = embedded.getHostEvidenceStart();
			int nextEvidenceEnd = embedded.getHostEvidenceEnd();
			int unknownLeft = nextEvidenceStart - embedded.getHostStart();
			int unknownRight = embedded.getHostEnd() - nextEvidenceEnd;
			if(sequenceId==idxDebug) System.out.println("Finding chimeras. Embedded "+embedded.getSequenceId()+" reverse"+embedded.isReverse()+" limits: "+nextEvidenceStart+" "+nextEvidenceEnd+" unknown: "+unknownLeft+" "+unknownRight+" count: "+embedded.getNumSharedKmers()+" CSK: "+embedded.getCoverageSharedKmers());
			hostPredictedEndLeft = Math.max(hostPredictedEndLeft, embedded.getHostEnd());
			hostPredictedStartRight = Math.min(hostPredictedStartRight, embedded.getHostStart());
			if(unknownRight>1000 && unknownLeft<1000) {
				hostEvidenceEndsLeft.add(nextEvidenceEnd);
			}
			if(unknownLeft>1000 && unknownRight<1000) {
				hostEvidenceStartsRight.add(nextEvidenceStart);
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
				if(edge.getEvidenceProportion()<0.5) continue;
				int nextEvidenceStart = (edge.getVertex1()==vS)?edge.getVertex1EvidenceStart():edge.getVertex2EvidenceStart();
				int nextEvidenceEnd = (edge.getVertex1()==vS)?edge.getVertex1EvidenceEnd():edge.getVertex2EvidenceEnd();
				int unknownLeft = nextEvidenceStart;
				int unknownRight = edge.getOverlap() - nextEvidenceEnd;
				if(sequenceId==idxDebug) System.out.println("Finding chimeras. Edge start "+edge+" evidence end: "+nextEvidenceEnd+" unknown: "+unknownLeft+" "+unknownRight+" count: "+edge.getNumSharedKmers()+" CSK: "+edge.getCoverageSharedKmers());
				hostPredictedEndLeft = Math.max(hostPredictedEndLeft, edge.getOverlap());
				if(unknownRight<0) continue;
				if(unknownRight>1000 && unknownLeft<1000) {
					hostEvidenceEndsLeft.add(nextEvidenceEnd);
					
				}
				if(unknownLeft>1000 && unknownRight<1000) {
					hostEvidenceStartsRight.add(nextEvidenceStart);
					numIncompleteEdgesLeft++;
					//hostPredictedStartRight = 0;
				}
			}
			for(AssemblyEdge edge: edgesE) {
				if(edge.isSameSequenceEdge()) continue;
				if(edge.getEvidenceProportion()<0.5) continue;
				int nextEvidenceStart = (edge.getVertex1()==vE)?edge.getVertex1EvidenceStart():edge.getVertex2EvidenceStart();
				int nextEvidenceEnd = (edge.getVertex1()==vE)?edge.getVertex1EvidenceEnd():edge.getVertex2EvidenceEnd();
				int unknownLeft = edge.getOverlap() - (seqLength-nextEvidenceStart);
				int unknownRight = seqLength - nextEvidenceEnd;
				if(sequenceId==idxDebug) System.out.println("Finding chimeras. Edge end "+edge+" evidence start: "+nextEvidenceStart+" unknown: "+unknownLeft+" "+unknownRight+" count: "+edge.getNumSharedKmers()+" CSK: "+edge.getCoverageSharedKmers());
				hostPredictedStartRight = Math.min(hostPredictedStartRight, seqLength-edge.getOverlap());
				if(unknownLeft<0) continue;
				if(unknownRight>1000 && unknownLeft<1000) {
					hostEvidenceEndsLeft.add(nextEvidenceEnd);
					numIncompleteEdgesRight++;
					//hostPredictedEndLeft = seqLength;
				}
				if(unknownLeft>1000 && unknownRight<1000) {
					hostEvidenceStartsRight.add(nextEvidenceStart);
				}
			}
		}
		
		if(sequenceId==idxDebug) System.out.println("Finding chimeras. Sequence "+sequenceId+". length "+seqLength+" num unknown: "+hostEvidenceEndsLeft.size()+" "+hostEvidenceStartsRight.size()+" ends left: "+hostEvidenceEndsLeft+" starts right: "+hostEvidenceStartsRight);
		//if(hostEvidenceStartsRight.size()==0 || hostEvidenceEndsLeft.size()==0 || hostEvidenceStartsRight.size()+numIncompleteEdgesRight<5 || hostEvidenceEndsLeft.size()+numIncompleteEdgesLeft<5 ) return false;
		//if(numIncompleteEdgesLeft<3 && numIncompleteEdgesRight<3 && (hostEvidenceStartsRight.size()<3 || hostEvidenceEndsLeft.size()<3)) return false;
		Collections.sort(hostEvidenceEndsLeft,(n1,n2)->n1-n2);
		int hostEvidenceEndLeft = hostEvidenceEndsLeft.size()>0?hostEvidenceEndsLeft.get(hostEvidenceEndsLeft.size()/2):0;
		Collections.sort(hostEvidenceStartsRight,(n1,n2)->n1-n2);
		int hostEvidenceStartRight = hostEvidenceStartsRight.size()>0?hostEvidenceStartsRight.get(hostEvidenceStartsRight.size()/2):0;
		if(hostEvidenceEndLeft==0 && hostEvidenceStartRight>0) hostEvidenceEndLeft = hostEvidenceStartRight;
		else if(hostEvidenceEndLeft>0 && hostEvidenceStartRight==0) hostEvidenceStartRight = hostEvidenceEndLeft;
		else if(hostEvidenceStartRight==0) hostEvidenceStartRight = seqLength;
		
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
		double limitGoodOverlap = 0.6*seqLength;
		double limitEvProp = 0.8;
		
		double maxEvPropS = 0;
		int countPassS = 0;
		int countGoodOverlapS = 0;
		for(AssemblyEdge edge: edgesS) {
			if(edge.isSameSequenceEdge()) continue;
			int nextEvidenceEnd = (edge.getVertex1()==vS)?edge.getVertex1EvidenceEnd():edge.getVertex2EvidenceEnd();
			int unknownLeft = (edge.getVertex1()==vS)?edge.getVertex1EvidenceStart():edge.getVertex2EvidenceStart();
			int unknownRight = edge.getOverlap() - nextEvidenceEnd;
			
			if(unknownLeft<200 && unknownRight<200 && nextEvidenceEnd-maxEvidenceEnd>100) {
				if(sequenceId==idxDebug) System.out.println("Edge start crossing. Edge: "+edge);
				numCrossing++;
			}
			if(edge.getOverlap()>limitGoodOverlap) {
				countGoodOverlapS++;
				maxEvPropS = Math.max(maxEvPropS, edge.getEvidenceProportion());
				if(edge.getEvidenceProportion()>=limitEvProp) countPassS++;
			}
		}
		double maxEvPropE = 0;
		int countPassE = 0;
		int countGoodOverlapE = 0;
		for(AssemblyEdge edge: edgesE) {
			if(edge.isSameSequenceEdge()) continue;
			int nextEvidenceStart = (edge.getVertex1()==vE)?edge.getVertex1EvidenceStart():edge.getVertex2EvidenceStart();
			int unknownLeft = edge.getOverlap() - (seqLength-nextEvidenceStart);
			int unknownRight = seqLength - ((edge.getVertex1()==vE)?edge.getVertex1EvidenceEnd():edge.getVertex2EvidenceEnd());
			
			if(unknownLeft<200 && unknownRight<200 && minEvidenceStart - nextEvidenceStart > 100) {
				if(sequenceId==idxDebug) System.out.println("Edge end crossing. Edge: "+edge);
				numCrossing++;
			}
			if(edge.getOverlap()>limitGoodOverlap) {
				countGoodOverlapE++;
				maxEvPropE = Math.max(maxEvPropE, edge.getEvidenceProportion());
				if(edge.getEvidenceProportion()>=limitEvProp) countPassE++;
			}
		}
		if(sequenceId==idxDebug) System.out.println("Finding chimeras. Sequence "+sequenceId+". length "+seqLength+" median evidence: "+hostEvidenceEndLeft+" "+hostEvidenceStartRight+" predicted: "+hostPredictedEndLeft+" "+hostPredictedStartRight+" num crossing: "+numCrossing+" incomplete: "+numIncompleteEdgesLeft+" "+numIncompleteEdgesRight+" edges good overlap: "+countGoodOverlapS+" "+countGoodOverlapE+" countpassEvProp: "+countPassS+" "+countPassE+" embedded size: "+embeddedList.size()+" is emb: "+isEmbedded(sequenceId));
		int d1 = hostEvidenceEndLeft-hostPredictedStartRight;
		int d2 = hostPredictedEndLeft-hostEvidenceStartRight;
		int d3 = hostPredictedEndLeft-hostEvidenceEndLeft;
		int d4 = hostEvidenceStartRight-hostPredictedStartRight;
		int d5 = hostEvidenceEndLeft - hostEvidenceStartRight;
		if( numCrossing<2  && d1>1000 && d2>1000 && d3>1000 && d4>1000 && d5<seqLength/2) {
			System.out.println("Possible chimera identified for sequence "+sequenceId+" "+getSequence(sequenceId).getName()+". length "+seqLength+" num unknown: "+hostEvidenceEndsLeft.size()+" "+hostEvidenceStartsRight.size()+" evidence end : "+hostEvidenceEndLeft+" "+hostEvidenceStartRight+" predicted: "+hostPredictedEndLeft+" "+hostPredictedStartRight+" crossing: "+numCrossing);
			return true;
		} /*else if ((countGoodOverlapS > 5 && countPassS ==0 && hostEvidenceEndsLeft.size()>0) || (countGoodOverlapE>5 && countPassE ==0 && hostEvidenceStartsRight.size()>0)) {
			System.out.println("Possible dangling end identified for sequence "+sequenceId+" "+getSequence(sequenceId).getName()+". length "+seqLength+" num unknown: "+hostEvidenceEndsLeft.size()+" "+hostEvidenceStartsRight.size()+" evidence end : "+hostEvidenceEndLeft+" "+hostEvidenceStartRight+" predicted: "+hostPredictedEndLeft+" "+hostPredictedStartRight+" crossing: "+numCrossing+" edges good overlap: "+countGoodOverlapS+" "+countGoodOverlapE+" countpassEvProp: "+countPassS+" "+countPassE);
			return true;
		} else if (numCrossing==0  && ((hostEvidenceStartsRight.size()>5 && countPassS<=numIncompleteEdgesLeft) || (hostEvidenceEndsLeft.size()>5 && countPassE<=numIncompleteEdgesRight))) {
			System.out.println("No evidence on one side for sequence "+sequenceId+" "+getSequence(sequenceId).getName()+". length "+seqLength+" num unknown: "+hostEvidenceEndsLeft.size()+" "+hostEvidenceStartsRight.size()+" evidence end : "+hostEvidenceEndLeft+" "+hostEvidenceStartRight+" predicted: "+hostPredictedEndLeft+" "+hostPredictedStartRight+" crossing: "+numCrossing+" incomplete: "+numIncompleteEdgesLeft+" "+numIncompleteEdgesRight+" countpassEvProp: "+countPassS+" "+countPassE);
			return true;
		} */
		
		/* else if (numCrossing==0 && hostEvidenceEndLeft==0 && hostEvidenceStartRight>1000 && numIncompleteEdgesLeft>5) {
			System.out.println("Possible chimera identified for start of sequence "+sequenceId+". length "+seqLength+" num unknown: "+hostEvidenceEndsLeft.size()+" "+hostEvidenceStartsRight.size()+" evidence end : "+hostEvidenceEndLeft+" "+hostEvidenceStartRight+" predicted: "+hostPredictedEndLeft+" "+hostPredictedStartRight+" num incomplete: "+numIncompleteEdgesLeft+" "+numIncompleteEdgesRight);
			return true;
		} else if (numCrossing==0 && hostEvidenceEndLeft<seqLength-1000 && hostEvidenceStartRight==seqLength && numIncompleteEdgesRight>5) {
			System.out.println("Possible chimera identified for end of sequence "+sequenceId+". length "+seqLength+" num unknown: "+hostEvidenceEndsLeft.size()+" "+hostEvidenceStartsRight.size()+" evidence end : "+hostEvidenceEndLeft+" "+hostEvidenceStartRight+" predicted: "+hostPredictedEndLeft+" "+hostPredictedStartRight+" num incomplete: "+numIncompleteEdgesLeft+" "+numIncompleteEdgesRight);
			return true;
		}*/
		
		return false;
	}
	
	public boolean isChimeric(int sequenceId) {
		return chimericSequenceIds.contains(sequenceId);
	}
	
	public List<AssemblyEdge> selectSafeEdges(  ) {
		Set<Integer> repetitiveVertices = predictRepetitiveVertices();
		return selectSafeEdges(repetitiveVertices);
	}
	
	private List<AssemblyEdge> selectSafeEdges( Set<Integer> repetitiveVertices ) {
		List<AssemblyEdge> allEdges = getEdges();
		List<AssemblyEdge> rawSafeEdges = new ArrayList<AssemblyEdge>();
		double averageCost = 0;
		double averageIKBP = 0;
		for(AssemblyEdge edge:allEdges) {
			if(isSafeEdge(edge, repetitiveVertices)) {
				rawSafeEdges.add(edge);
				averageCost+=edge.getCost();
				averageIKBP+=edge.getIndelsPerKbp();
			}
		}
		if (rawSafeEdges.size()==0) return rawSafeEdges;
		averageCost/=rawSafeEdges.size();
		averageIKBP/=rawSafeEdges.size();
		List<AssemblyEdge> filteredSafeEdges = new ArrayList<AssemblyEdge>();
		for(AssemblyEdge edge:rawSafeEdges) {
			if(edge.getCost()<=3*averageCost && edge.getIndelsPerKbp()<=5*averageIKBP) filteredSafeEdges.add(edge);
		}
		return filteredSafeEdges;
	}
	private boolean isSafeEdge(AssemblyEdge edge, Set<Integer> repetitiveVertices) {
		if(edge.isSameSequenceEdge()) return false;
		if(isEmbedded(edge.getVertex1().getSequenceIndex())) return false;
		if(isEmbedded(edge.getVertex2().getSequenceIndex())) return false;
		boolean r1 = repetitiveVertices.contains(edge.getVertex1().getUniqueNumber());
		boolean r2 = repetitiveVertices.contains(edge.getVertex2().getUniqueNumber()); 
		int d1 = getEdges(edge.getVertex1()).size();
		int d2 = getEdges(edge.getVertex2()).size();
		//if(edge.getVertex1().getUniqueNumber()==-97856) System.out.println("Select safe edges. repetitive: "+r1+" "+r2+" initial degrees "+edge.getVertex1().getDegreeUnfilteredGraph()+" "+edge.getVertex2().getDegreeUnfilteredGraph()+" current degrees "+d1+" "+d2+" reciprocal best: "+isRecipocalBest(edge)+" edge: "+edge);
		if((d1==2 && d2==2) || (!r1 && !r2 && isRecipocalBest(edge))) return true;
		return false;
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
	
	public boolean isRecipocalBestByCost(AssemblyEdge edge) {
		AssemblyVertex v1 = edge.getVertex1();
		if(!isBestEdgeByCost(v1, edge)) return false;
		AssemblyVertex v2 = edge.getVertex2();
		if(!isBestEdgeByCost(v2, edge)) return false;
		return true;
	}
	public boolean isBestEdgeByCost (AssemblyVertex v, AssemblyEdge edgeT) {
		List<AssemblyEdge> edges = getEdges(v);
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge() || edge==edgeT) continue;
			if(edge.getCost()<edgeT.getCost()) return false;
		}
		return true;
	}
	public AssemblyEdge getEdgeBestOverlap(AssemblyVertex vertex) {
		List<AssemblyEdge> edges = getEdges(vertex);
		AssemblyEdge maxOverlapEdge = null;
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge()) continue;
			if(maxOverlapEdge==null || maxOverlapEdge.getOverlap()<edge.getOverlap()) maxOverlapEdge = edge;
		}
		return maxOverlapEdge;
	}
	public AssemblyEdge getEdgeMinCost(AssemblyVertex vertex) {
		List<AssemblyEdge> edges = getEdges(vertex);
		AssemblyEdge minCostEdge = null;
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge()) continue;
			if(minCostEdge==null || minCostEdge.getCost()>edge.getCost()) minCostEdge = edge;
		}
		return minCostEdge;
	}
	
	public NormalDistribution[] estimateDistributions(List<AssemblyEdge> edges, Set<Integer> repetitiveVertices) {
		Distribution overlapDistributionSafe = new Distribution(0, 100000, 1000);
		Distribution wcskDistributionSafe = new Distribution(0, 100000, 1000);
		Distribution wcskPropOverlapSafe = new Distribution(0, 1.1, 0.01);
		Distribution overlapSDDistributionSafe = new Distribution(0, 500, 10);
		Distribution evPropDistributionSafe = new Distribution(0, 1.1, 0.02);
		Distribution indelsKbpDistributionSafe = new Distribution(0, 100, 1);
		Distribution [] distsSafe = {overlapDistributionSafe, wcskDistributionSafe, wcskPropOverlapSafe, overlapSDDistributionSafe,evPropDistributionSafe,indelsKbpDistributionSafe};
		Distribution overlapDistributionAll = new Distribution(0, 100000, 1000);
		Distribution wcskDistributionAll = new Distribution(0, 100000, 1000);
		Distribution wcskPropOverlapAll = new Distribution(0, 1.1, 0.01);
		Distribution overlapSDDistributionAll = new Distribution(0, 50, 10);
		Distribution evPropDistributionAll = new Distribution(0, 1.1, 0.02);
		Distribution indelsKbpDistributionAll = new Distribution(0, 100, 1);
		Distribution [] distsAll = {overlapDistributionAll, wcskDistributionAll, wcskPropOverlapAll, overlapSDDistributionAll,evPropDistributionAll,indelsKbpDistributionAll};
		for(AssemblyEdge edge:edges) {
			if (edge.isSameSequenceEdge()) continue;
			double overlap = edge.getOverlap();
			overlapDistributionAll.processDatapoint(overlap);
			wcskDistributionAll.processDatapoint(edge.getWeightedCoverageSharedKmers());
			wcskPropOverlapAll.processDatapoint((double)edge.getWeightedCoverageSharedKmers()/(edge.getOverlap()+1));
			overlapSDDistributionAll.processDatapoint((double)edge.getOverlapStandardDeviation());
			evPropDistributionAll.processDatapoint(edge.getEvidenceProportion());
			indelsKbpDistributionAll.processDatapoint(edge.getIndelsPerKbp());
			if (isSafeEdge(edge, repetitiveVertices)) {
				overlapDistributionSafe.processDatapoint(overlap);
				wcskDistributionSafe.processDatapoint(edge.getWeightedCoverageSharedKmers());
				wcskPropOverlapSafe.processDatapoint((double)edge.getWeightedCoverageSharedKmers()/(edge.getOverlap()+1));
				overlapSDDistributionSafe.processDatapoint((double)edge.getOverlapStandardDeviation());
				evPropDistributionSafe.processDatapoint(edge.getEvidenceProportion());
				indelsKbpDistributionSafe.processDatapoint(edge.getIndelsPerKbp());
			}
		}
		/*System.out.println("Dist overlap");
		overlapDistributionAll.printDistribution(System.out);
		System.out.println("Dist WCSK");
		wcskDistributionAll.printDistribution(System.out);
		System.out.println("Dist WCSK prop");
		wcskPropOverlapAll.printDistribution(System.out);
		System.out.println("Dist ev prop");
		evPropDistributionAll.printDistribution(System.out);*/
		double numSafe = overlapDistributionSafe.getCount();
		System.out.println("Number of safe edges: "+numSafe);
		NormalDistribution [] answer = new NormalDistribution[distsAll.length];
		for(int i=0;i<distsAll.length;i++) {
			//TODO: Better mean when there are no safe edges
			double mean;
			if(numSafe>20) {
				mean = distsSafe[i].getLocalMode(distsSafe[i].getAverage()/2, distsSafe[i].getAverage()*2);
			} else if (i<5) {
				mean = distsAll[i].getLocalMode(distsAll[i].getAverage(), distsAll[i].getMaxValueDistribution());
				if(i<2) mean*=2;
			} else {
				mean = distsAll[i].getLocalMode(distsAll[i].getMinValueDistribution(), distsAll[i].getAverage());
			}
			if(i==5 ) mean = Math.max(mean, 0.5);
			double variance;
			if(i!=5) variance = Math.pow(distsAll[i].getEstimatedStandardDeviationPeak(mean),2);
			else if (numSafe>20) variance = distsSafe[i].getVariance();
			else variance = distsAll[i].getVariance();
			if(numSafe>20) variance = Math.max(variance, distsSafe[i].getVariance());
			// For overlap, WCSK and IKBP the variance should be at least the mean
			if((i<2 || i==5) && variance < mean) variance = mean;
			//For proportions, the stdev should be at least 0.1
			if(i==2 || i==4) variance = Math.max(variance, 0.01);
			 
			answer[i] = new NormalDistribution(mean,variance);
		}
		return answer;
	}
	
	public void updateScores (double weightsSecondaryFeatures) {
		updateVertexDegrees();
		Set<Integer> repetitiveVertices = predictRepetitiveVertices();
		NormalDistribution [] edgesDists = estimateDistributions(getEdges(),repetitiveVertices);
		System.out.println("Average overlap: "+edgesDists[0].getMean()+" SD: "+Math.sqrt(edgesDists[0].getVariance()));
		System.out.println("Average weighted coverage shared kmers: "+edgesDists[1].getMean()+" SD: "+Math.sqrt(edgesDists[1].getVariance()));
		System.out.println("Average WCSK proportion overlap: "+edgesDists[2].getMean()+" SD: "+Math.sqrt(edgesDists[2].getVariance()));
		System.out.println("Average overlap standard deviation: "+edgesDists[3].getMean()+" SD: "+Math.sqrt(edgesDists[3].getVariance()));
		System.out.println("Average Evidence proportion: "+edgesDists[4].getMean()+" SD: "+Math.sqrt(edgesDists[4].getVariance()));
		System.out.println("Average indels kbp: "+edgesDists[5].getMean()+" SD: "+Math.sqrt(edgesDists[5].getVariance()));
		AssemblySequencesRelationshipScoresCalculator calculator = new AssemblySequencesRelationshipScoresCalculator();
		calculator.setWeightsSecondaryFeatures(weightsSecondaryFeatures);
		List<AssemblyEdge> allEdges = getEdges();
		for(AssemblyEdge edge: allEdges) {
			edge.setScore(calculator.calculateScore(edge,edgesDists));
			edge.setCost(calculator.calculateCost(edge,edgesDists));
		}
		for (List<AssemblyEmbedded> embeddedList:embeddedMapBySequence.values()) {
			for(AssemblyEmbedded embedded:embeddedList) {
				embedded.setScore(calculator.calculateScore(embedded, edgesDists));
				embedded.setCost(calculator.calculateCost(embedded, edgesDists));
			}
		}
	}
	public void addRelationship(AssemblySequencesRelationship rel) {
		if(rel instanceof AssemblyEmbedded) addEmbedded((AssemblyEmbedded)rel);
		if(rel instanceof AssemblyEdge) addEdge((AssemblyEdge)rel);
		
	}
	public void removeRelationship(AssemblySequencesRelationship rel) {
		if(rel instanceof AssemblyEdge) removeEdge((AssemblyEdge)rel);
		if(rel instanceof AssemblyEmbedded) removeEmbedded((AssemblyEmbedded)rel);
		
		
	}
	
}
