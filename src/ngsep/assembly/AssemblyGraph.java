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
		for (int i=0;i<sequences.size();i++) {
			QualifiedSequence seq = sequences.get(i);
			AssemblyVertex vS = new AssemblyVertex(seq, true, i);
			verticesStart.put(i,vS);
			verticesByUnique.put(vS.getUniqueNumber(), vS);
			edgesMap.put(vS.getUniqueNumber(), new ArrayList<>());
			AssemblyVertex vE = new AssemblyVertex(seq, false, i);
			verticesEnd.put(i,vE);
			verticesByUnique.put(vE.getUniqueNumber(), vE);
			edgesMap.put(vE.getUniqueNumber(), new ArrayList<>());
			AssemblyEdge edge = new AssemblyEdge(vS, vE, seq.getLength());
			edge.setAverageOverlap(seq.getLength());
			edge.setMedianOverlap(seq.getLength());
			edge.setFromLimitsOverlap(seq.getLength());
			edge.setCoverageSharedKmers(seq.getLength());
			edge.setWeightedCoverageSharedKmers(seq.getLength());
			edge.setNumSharedKmers(seq.getLength());
			edge.setOverlapStandardDeviation(0);
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
		int firstEvidence = -1;
		int lastEvidence = -1;
		int lastUnknownRight = -1;
		/*AssemblyVertex vS = verticesStart.get(sequenceId);
		AssemblyVertex vE = verticesEnd.get(sequenceId);
		for(AssemblyEdge edge:edgesMap.get(vS.getUniqueNumber())) {
			KmerHitsCluster cluster = edge.getEvidence();
			if(edge.isSameSequenceEdge() || cluster==null) continue;
			if(firstEvidence==-1) {
				firstEvidence = 1;
			}
			int unknown = cluster.getLast()-cluster.getSubjectEvidenceEnd();
			
			lastEvidence = Math.max(lastEvidence, cluster.getSubjectEvidenceEnd());
		}*/
		if(sequenceId==idxDebug) System.out.println("Finding chimeras. Left edges evidence limits: "+firstEvidence+" "+lastEvidence);
		List<AssemblyEmbedded> embeddedList = new ArrayList<AssemblyEmbedded>();
		List<AssemblyEmbedded> emb = embeddedMapByHost.get(sequenceId);
		if(emb==null) return false;
		if(sequenceId==idxDebug) System.out.println("Finding chimeras. Embedded sequences "+emb.size());
		embeddedList.addAll(emb);
		Collections.sort(embeddedList,(e1,e2)->e1.getHostEvidenceStart()-e2.getHostEvidenceStart());
		
		for(AssemblyEmbedded embedded:embeddedList) {
			
			int nextLeft = embedded.getHostEvidenceStart();
			int nextRight = embedded.getHostEvidenceEnd();
			int unknownLeft = nextLeft - embedded.getHostStart();
			int unknownRight = embedded.getHostEnd() - nextRight;
			if(sequenceId==idxDebug) System.out.println("Finding chimeras. Last evidence: "+lastEvidence+" Embedded "+embedded.getSequenceId()+" reverse"+embedded.isReverse()+" limits: "+nextLeft+" "+nextRight+" unknown: "+unknownLeft+" "+unknownRight+" count: "+embedded.getNumSharedKmers());
			if(firstEvidence==-1) {
				firstEvidence = nextLeft;
				lastEvidence = nextRight;
			}
			if(nextLeft>lastEvidence) {
				int distance = nextLeft-lastEvidence;
				if(lastEvidence>1000 && seqLength-nextLeft >1000 && lastUnknownRight>2*distance && unknownLeft>2*distance) {
					System.out.println("Possible chimera identified for sequence "+sequenceId+". length "+seqLength+" last evidence: "+lastEvidence+" next: "+nextLeft+" unknown limits : "+lastUnknownRight+" "+unknownLeft);
					return true;
				}
			}
			if(lastEvidence<nextRight) {
				lastEvidence = nextRight;
				lastUnknownRight = unknownRight;
			}	
		}
		//TODO: check end
		
		return false;
	}
	
	public void filterEdgesAndEmbedded() {
		int medianLength = getMedianLength();
		System.out.println("Median read length: "+medianLength);
		for (int seqId = sequences.size()-1; seqId >=0; seqId--) {
			filterEdgesAndEmbedded(seqId,medianLength);
		}
		System.out.println("Filtered edges. Prunning embedded");
		pruneEmbeddedSequences();
		System.out.println("Prunned embedded sequences");
		//filterEdgesCloseRelationships();
	}

	public void filterEdgesAndEmbedded(int sequenceId,int medianLength) {
		int debugIdx = -1;
		int sequenceLength = getSequenceLength(sequenceId);
		AssemblyVertex vS = verticesStart.get(sequenceId);
		AssemblyVertex vE = verticesEnd.get(sequenceId);
		if(vS==null || vE==null) return;
		filterEdgesAbnormalFeatures(getEdges(vS));
		filterEdgesAbnormalFeatures(getEdges(vE));
		if(sequenceId == debugIdx) System.out.println("Filtered edges with abnormal features");
		List<AssemblyEdge> edgesS = new ArrayList<AssemblyEdge>();
		if(vS!=null) edgesS.addAll(getEdges(vS));
		double minScoreProportionEdges = 0.5;
		
		double maxScoreSE = 0;
		double maxScoreSF = 0;
		for(AssemblyEdge edge: edgesS) {
			if(edge.isSameSequenceEdge()) continue;
			maxScoreSE = Math.max(maxScoreSE, calculateScoreForEmbedded(edge));
			maxScoreSF = Math.max(maxScoreSF, calculateScoreForEdgeFiltering(edge));
		}
		for(AssemblyEdge edge: edgesS) {
			if(edge.isSameSequenceEdge()) continue;
			double score = calculateScoreForEdgeFiltering(edge);
			if(sequenceId == debugIdx) System.out.println("Assembly graph. Next edge start "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber()+" overlap: "+edge.getOverlap()+" score: "+score+" Max score start: "+maxScoreSF);
			if(score < minScoreProportionEdges*maxScoreSF) {
				if(sequenceId == debugIdx) System.out.println("Assembly graph. Removing edge: "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber());
				removeEdge(edge);
			}
		}
		if(sequenceId == debugIdx) System.out.println("Assembly graph. Initial edges start "+edgesS.size()+" Max scores start: "+maxScoreSF+" "+maxScoreSE+" remaining edges: "+edgesMap.get(vS.getUniqueNumber()).size());
		List<AssemblyEdge> edgesE = new ArrayList<AssemblyEdge>();
		if(vE!=null) edgesE.addAll(getEdges(vE));
		double maxScoreEE = 0;
		double maxScoreEF = 0;			
		for(AssemblyEdge edge: edgesE) {
			if(edge.isSameSequenceEdge()) continue;
			maxScoreEE = Math.max(maxScoreEE, calculateScoreForEmbedded(edge));
			maxScoreEF = Math.max(maxScoreEF, calculateScoreForEdgeFiltering(edge));
		}
		for(AssemblyEdge edge: edgesE) {
			if(edge.isSameSequenceEdge()) continue;
			double score = calculateScoreForEdgeFiltering(edge);
			if(sequenceId == debugIdx) System.out.println("Assembly graph. Next edge end "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber()+" overlap: "+edge.getOverlap()+" score: "+score+" Max score end: "+maxScoreEF);
			if(score < minScoreProportionEdges*maxScoreEF) {
				if(sequenceId == debugIdx) System.out.println("Assembly graph. Removing edge: "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber());
				removeEdge(edge);
			}
		}
		if(sequenceId == debugIdx) System.out.println("Assembly graph. Initial edges end "+edgesE.size()+" Max scores end: "+maxScoreEF+" "+maxScoreEE+" remaining edges: "+edgesMap.get(vE.getUniqueNumber()).size());
		
		//double minScoreProportionEmbedded = 0.8;
		//double minScoreProportionEmbedded = Math.min(0.9, (double)getSequenceLength(sequenceId)/50000.0);
		double minScoreProportionEmbedded = Math.min(0.9, 0.4*sequenceLength/(double)medianLength);
		if(minScoreProportionEmbedded<0.5) minScoreProportionEmbedded = 0.5;
		
		double maxScoreFilterEmbedded = minScoreProportionEmbedded*Math.max(maxScoreSE, maxScoreEE);
		List<AssemblyEmbedded> embeddedList= new ArrayList<AssemblyEmbedded>();
		embeddedList.addAll(getEmbeddedBySequenceId(sequenceId));
		if(embeddedList.size()==0) return;
		double maxScoreEmbedded = -1;
		for(AssemblyEmbedded embedded:embeddedList) {
			if(sequenceId == debugIdx) System.out.println("Assembly graph. Next embedded "+embedded.getHostId()+" limits: "+embedded.getHostStart()+" "+embedded.getHostEnd()+" score: "+calculateScore(embedded));
			maxScoreEmbedded = Math.max(maxScoreEmbedded, calculateScore(embedded));
		}
		
		AssemblyEdge sameSequenceEdge = getSameSequenceEdge(sequenceId);
		if(maxScoreEmbedded<maxScoreFilterEmbedded) {
		//if(maxScoreEmbedded<maxScoreFilterEmbedded || maxScoreEmbedded < 0.2*calculateScore(sameSequenceEdge)) {
		//if(maxScoreEmbedded < 0.2*calculateScore(sameSequenceEdge)) {
			//Replace embedded relationships with edges to make the sequence not embedded
			for(AssemblyEmbedded embedded:embeddedList) {
				removeEmbedded(embedded);
				if(sequenceId == debugIdx) System.out.println("Adding edge replacing embedded "+embedded.getHostId()+" limits: "+embedded.getHostStart()+" "+embedded.getHostEnd()+" host length: "+getSequenceLength(embedded.getHostId())+"score: "+calculateScore(embedded));
				addEdgeFromEmbedded(embedded);
			}
		} else {
			if(sequenceId == debugIdx) System.out.println("Assembly graph. Sequence is embedded. Max score to filter "+ maxScoreFilterEmbedded +" max score embedded "+maxScoreEmbedded+" minprop "+minScoreProportionEmbedded);
			filterEmbedded(sequenceId, 0.9, 1);
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
		AssemblyEdge edge = new AssemblyEdge(vertexHost, vertexEmbedded, getSequenceLength(embedded.getSequenceId()-1));
		edge.setOverlapStandardDeviation(embedded.getHostStartStandardDeviation());
		edge.setWeightedCoverageSharedKmers(embedded.getWeightedCoverageSharedKmers());
		edge.setCoverageSharedKmers(embedded.getCoverageSharedKmers());
		edge.setNumSharedKmers(embedded.getNumSharedKmers());
		edge.setRawKmerHits(embedded.getRawKmerHits());
		edge.setRawKmerHitsSubjectStartSD(embedded.getRawKmerHitsSubjectStartSD());
		addEdge(edge);
	}

	public void filterEmbedded(int sequenceId, double minProportion, int numRetain) {
		List<AssemblyEmbedded> embeddedList= new ArrayList<AssemblyEmbedded>();
		embeddedList.addAll(getEmbeddedBySequenceId(sequenceId));
		if(embeddedList.size()==0) return;
		Collections.sort(embeddedList,(e1,e2)->(int)Math.round(calculateScore(e2)-calculateScore(e1)));
		double maxScore = calculateScore(embeddedList.get(0));
		int count = 0;
		for(AssemblyEmbedded embedded:embeddedList) {
			count++;
			if(count>numRetain || calculateScore(embedded)<minProportion*maxScore) {
				removeEmbedded(embedded);
				//if(sequenceId == debugIdx) System.out.println("Assembly graph. Removed embedded host: "+embedded.getHostId()+" Embedded relations: "+embeddedMapBySequence.get(sequenceId)+" is embedded: "+isEmbedded(sequenceId));
			}
		}
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
		return embedded.getCoverageSharedKmers();
		//return embedded.getRawKmerHits();
		//return embedded.getWeightedCoverageSharedKmers();
		//return embedded.getWeightedCoverageSharedKmers()*embedded.getRawKmerHits();
	}

	private double calculateScoreForEmbedded(AssemblyEdge edge) {
		return edge.getCoverageSharedKmers();
		//return edge.getRawKmerHits();
		//return edge.getWeightedCoverageSharedKmers();
		//return edge.getWeightedCoverageSharedKmers()*edge.getRawKmerHits();
	}
	
	private double calculateScoreForEdgeFiltering(AssemblyEdge edge) {
		//return edge.getCoverageSharedKmers();
		//return edge.getRawKmerHits();
		return edge.getWeightedCoverageSharedKmers();
		//return edge.getWeightedCoverageSharedKmers()*edge.getRawKmerHits();
	}

	public void filterEdgesCloseRelationships() {
		for (int seqId = 0; seqId <sequences.size(); seqId++) {
			AssemblyVertex v1 = getVertex(seqId, true);
			if (v1!=null) filterEdgesCloseRelationships(v1);
			AssemblyVertex v2 = getVertex(seqId, false);
			if (v2!=null) filterEdgesCloseRelationships(v2);
		}
		
	}

	private void filterEdgesCloseRelationships(AssemblyVertex vertex) {
		int debugIdx = -1;
		List<AssemblyEdge> edges = getEdges(vertex);
		if(vertex.getSequenceIndex()==debugIdx) System.out.println("Filter edges close. Vertex: "+vertex+" Total edges: "+edges.size());
		if(edges.size()<3) return;
		AssemblyEdge maxOverlapEdge = null;
		AssemblyEdge maxScoreEdge = null;
		for(AssemblyEdge edge: edges) {
			if(edge.isSameSequenceEdge()) continue;
			if(maxOverlapEdge==null || maxOverlapEdge.getOverlap() < edge.getOverlap()) {
				maxOverlapEdge = edge;
			}
			if(maxScoreEdge==null || calculateScoreForEdgeFiltering(maxScoreEdge)<calculateScoreForEdgeFiltering(edge)) {
				maxScoreEdge = edge;
			}
		}
		AssemblyEdge secondOverlapEdge = null;
		for(AssemblyEdge edge: edges) {
			if(edge.isSameSequenceEdge() || edge == maxOverlapEdge) continue;
			if(secondOverlapEdge==null || secondOverlapEdge.getOverlap() < edge.getOverlap()) {
				secondOverlapEdge = edge;
			}
		}
		if(secondOverlapEdge==null) System.err.println("Error. Vertex "+vertex+" has three edges but no second overlap. Max overlap: "+maxOverlapEdge+" edges: "+edges);
		if(vertex.getSequenceIndex()==debugIdx) System.out.println("Filter edges close. Best edge: "+maxOverlapEdge.getConnectingVertex(vertex).getUniqueNumber()+" second: "+secondOverlapEdge.getConnectingVertex(vertex).getUniqueNumber());
		int ov1 = maxOverlapEdge.getOverlap();
		int ov2 = secondOverlapEdge.getOverlap();
		int diff = ov1-ov2;
		if(vertex.getSequenceIndex()==debugIdx) System.out.println("Filter edges close. overlap 1: "+ov1+" overlap2: "+ov2+" diff: "+diff);
		if(maxOverlapEdge==maxScoreEdge && diff > 0.1*ov1) return;
		if(badTransitivity(vertex, maxOverlapEdge, secondOverlapEdge)) {
			if(vertex.getSequenceIndex()==debugIdx) System.out.println("Filter edges close. removing edge "+maxOverlapEdge);
			removeEdge(maxOverlapEdge);
		} else if (maxOverlapEdge!=maxScoreEdge && badTransitivity(vertex,maxScoreEdge,maxOverlapEdge)) {
			if(vertex.getSequenceIndex()==debugIdx) System.out.println("Filter edges close. removing edge "+maxScoreEdge);
			removeEdge(maxScoreEdge);
		}
	}

	private boolean badTransitivity(AssemblyVertex vertex, AssemblyEdge edge1, AssemblyEdge edge2) {
		AssemblyVertex v1 = edge1.getConnectingVertex(vertex);
		AssemblyEdge seqEdge1 = getSameSequenceEdge(v1);
		AssemblyVertex v12 = seqEdge1.getConnectingVertex(v1);
		AssemblyVertex v2 = edge2.getConnectingVertex(vertex);
		AssemblyEdge seqEdge2 = getSameSequenceEdge(v2);
		AssemblyVertex v22 = seqEdge2.getConnectingVertex(v2);
		AssemblyEdge transitiveGood = getEdge(v12, v2);
		AssemblyEdge transitiveBad = getEdge(v22, v1);
		//if(vertex.getSequenceIndex()==debugIdx) System.out.println("Filter edges close. transitive good: "+transitiveGood+" transitive bad: "+transitiveBad);
		return transitiveBad!=null && transitiveGood==null;
	}
	
}
