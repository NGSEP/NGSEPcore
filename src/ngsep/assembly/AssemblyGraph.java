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

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ngsep.math.Distribution;
import ngsep.sequences.KmerHitsCluster;
import ngsep.sequences.QualifiedSequence;

/**
 * @author Jorge Duitama
 * @author Juan Camilo Bojaca
 * @author David Guevara
 */
public class AssemblyGraph implements Serializable {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -4780291335510617705L;
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
			edge.setMismatches(0);
			edge.setCoverageSharedKmers(seq.getLength());
			addEdge(edge);
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

	public void setReadNames(List<String> readNames) {
		this.readNames = readNames;
	}

	public void addEdge(AssemblyEdge edge) {
		edgesMap.get(edge.getVertex1().getUniqueNumber()).add(edge);
		edgesMap.get(edge.getVertex2().getUniqueNumber()).add(edge);
		numEdges++;
	}
	
	public void removeEdge (AssemblyEdge edge) {
		edgesMap.get(edge.getVertex1().getUniqueNumber()).remove(edge);
		edgesMap.get(edge.getVertex2().getUniqueNumber()).remove(edge);
		numEdges--;
	}

	public AssemblyVertex getVertex(int indexSequence, boolean start) {
		if(start) return verticesStart.get(indexSequence);
		return verticesEnd.get(indexSequence);
	}
	
	public AssemblyVertex getVertexByUniqueId(int uniqueId) {
		return verticesByUnique.get(uniqueId);
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
	
	public void pruneEmbeddedSequences() {
		for(int i:embeddedMapBySequence.keySet()) {
			if(verticesStart.get(i)!=null) {
				removeEdges(i);
			}
		}
	}
	
	public void selectEmbeddedSequences() {
		List<Integer> currentEmbedded = new ArrayList<Integer>();
		currentEmbedded.addAll(embeddedMapBySequence.keySet());
		Collections.sort(currentEmbedded);
		Set<Integer> freezedEmbedded = new HashSet<Integer>();
		for(int i:currentEmbedded) {
			List<AssemblyEmbedded> embeddedRelations = embeddedMapBySequence.get(i);
			AssemblyEmbedded toRemove = null;
			if(!freezedEmbedded.contains(i) && embeddedRelations.size()==1) {
				AssemblyEmbedded embedded = embeddedRelations.get(0);
				if (makeSpecialEdge (embedded)) {
					toRemove = embedded;
				}
			}
			if(toRemove==null) {
				if (verticesStart.get(i)!=null) removeVertices(i);
			} else {
				removeEmbedded(toRemove);
				List<AssemblyEmbedded> children = embeddedMapByHost.get(i);
				if(children==null) return;
				for(AssemblyEmbedded child:children) {
					freezedEmbedded.add(child.getSequenceId());
				}
			}
		}
	}

	private boolean makeSpecialEdge(AssemblyEmbedded embedded) {
		//TODO: Parameters instead of hard filters
		int seqIdHost = embedded.getHostId();
		CharSequence host = sequences.get(seqIdHost).getCharacters();
		int seqIdEmb = embedded.getSequenceId();
		CharSequence seq = sequences.get(seqIdEmb).getCharacters();
		if(seq.length()<5000) return false;
		if(seq==null || host == null) return false;
		if (seq.length()<0.7*host.length()) return false;
		if (seq.length()>0.95*host.length()) return false;
		if(embeddedMapBySequence.containsKey(embedded.getHostId())) return false;
		int left = embedded.getStartPosition();
		int right = host.length()-(left+seq.length());
		if(Math.abs(left-right)<100) return false;
		AssemblyVertex vHost = null;
		AssemblyVertex vSeq = null;
		if (left < right) {
			if (left > 100) return false;
			vHost = verticesStart.get(seqIdHost);
			if (embedded.isReverse()) vSeq = verticesStart.get(seqIdEmb);
			else vSeq = verticesEnd.get(seqIdEmb);
		} else {
			if (right > 100) return false;
			vHost = verticesEnd.get(seqIdHost);
			if (embedded.isReverse()) vSeq = verticesEnd.get(seqIdEmb);
			else vSeq = verticesStart.get(seqIdEmb);
		}
		if(vHost==null || vSeq == null ) return false;
		if(getEdge(vHost, vSeq)==null) {
			AssemblyEdge specialEdge = new AssemblyEdge(vHost, vSeq, seq.length());
			specialEdge.setMismatches(embedded.getMismatches());
			specialEdge.setCoverageSharedKmers(embedded.getCoverageSharedKmers());
			addEdge(specialEdge);
		}
		System.out.println("Special edge created for embedded sequence: "+seqIdEmb+" with length "+seq.length()+" in host: "+seqIdHost+" length "+host.length()+" leftover first: "+left+" leftover right: "+ right);
		return true;
	}

	public void removeEdges(int sequenceId) {
		AssemblyVertex v1 = getVertex(sequenceId, true);
		List<AssemblyEdge> edgesToRemove = new ArrayList<>(); 
		edgesToRemove.addAll(edgesMap.get(v1.getUniqueNumber()));
		AssemblyVertex v2 = getVertex(sequenceId, false);
		edgesToRemove.addAll(edgesMap.get(v2.getUniqueNumber()));
		for(AssemblyEdge edge:edgesToRemove) {
			removeEdge(edge);
		}
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

	public void addPath(List<AssemblyEdge> path) {
		paths.add(path);
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
		throw new RuntimeException("Same sequence edge not found for vertex: "+vertex.getSequenceIndex()+"-"+vertex.isStart());
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


	public void serialize(String outFileGraph) throws IOException {
		try (FileOutputStream fos = new FileOutputStream( outFileGraph );
			 ObjectOutputStream oos = new ObjectOutputStream( fos );) {
			oos.writeObject( this );
		}
	}
	
	public static AssemblyGraph load(String filename) throws IOException {
		try (FileInputStream fis = new FileInputStream(filename);
			 ObjectInputStream ois = new ObjectInputStream(fis);) {
			return (AssemblyGraph) ois.readObject();
		} catch (ClassNotFoundException e) {
			throw new RuntimeException("AssemblyGraph class not found",e);
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
		Collections.sort(embeddedList,(e1,e2)->e1.getEvidence().getSubjectEvidenceStart()-e2.getEvidence().getSubjectEvidenceStart());
		
		for(AssemblyEmbedded embedded:embeddedList) {
			KmerHitsCluster cluster = embedded.getEvidence();
			if(cluster==null) continue;
			
			int nextLeft = cluster.getSubjectEvidenceStart();
			int nextRight = cluster.getSubjectEvidenceEnd();
			int unknownLeft = nextLeft - cluster.getSubjectPredictedStart();
			int unknownRight = cluster.getSubjectPredictedEnd() - nextRight;
			if(sequenceId==idxDebug) System.out.println("Finding chimeras. Last evidence: "+lastEvidence+" Embedded "+embedded.getSequenceId()+" reverse"+embedded.isReverse()+" limits: "+nextLeft+" "+nextRight+" unknown: "+unknownLeft+" "+unknownRight+" count: "+cluster.getNumDifferentKmers()+" weighted: "+cluster.getWeightedCount());
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


	public void filterEdgesAndEmbedded(int sequenceId) {
		int debugIdx = -1;
		AssemblyVertex vS = verticesStart.get(sequenceId);
		AssemblyVertex vE = verticesEnd.get(sequenceId);
		filterEdgesAbnormalOverlap(getEdges(vS));
		filterEdgesAbnormalOverlap(getEdges(vE));
		List<AssemblyEdge> edgesS = new ArrayList<AssemblyEdge>();
		if(vS!=null) edgesS.addAll(getEdges(vS));
		double minScoreProportion = 0.5;
		double maxScoreS = 0;			
		for(AssemblyEdge edge: edgesS) {
			if(edge.isSameSequenceEdge()) continue;
			double score = calculateScore(edge);
			if(score > maxScoreS) {
				maxScoreS = score;
			}
		}
		for(AssemblyEdge edge: edgesS) {
			if(edge.isSameSequenceEdge()) continue;
			double score = calculateScore(edge);
			if(sequenceId == debugIdx) System.out.println("Assembly graph. Next edge start "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber()+" overlap: "+edge.getOverlap()+" score: "+score+" Max score start: "+maxScoreS);
			if(score < minScoreProportion*maxScoreS) {
				if(sequenceId == debugIdx) System.out.println("Assembly graph. Removing edge: "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber());
				removeEdge(edge);
			}
		}
		if(sequenceId == debugIdx) System.out.println("Assembly graph. Initial edges start "+edgesS.size()+" Max score end: "+maxScoreS+" remaining edges: "+edgesMap.get(vS.getUniqueNumber()).size());
		List<AssemblyEdge> edgesE = new ArrayList<AssemblyEdge>();
		if(vE!=null) edgesE.addAll(getEdges(vE));
		double maxScoreE = 0;			
		for(AssemblyEdge edge: edgesE) {
			if(edge.isSameSequenceEdge()) continue;
			double score = calculateScore(edge);
			if(score > maxScoreE) {
				maxScoreE = score;
			}
		}
		for(AssemblyEdge edge: edgesE) {
			if(edge.isSameSequenceEdge()) continue;
			double score = calculateScore(edge);
			if(sequenceId == debugIdx) System.out.println("Assembly graph. Next edge end "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber()+" overlap: "+edge.getOverlap()+" score: "+score+" Max score end: "+maxScoreE);
			if(score < minScoreProportion*maxScoreE) {
				if(sequenceId == debugIdx) System.out.println("Assembly graph. Removing edge: "+edge.getVertex1().getUniqueNumber()+" "+edge.getVertex2().getUniqueNumber());
				removeEdge(edge);
			}
		}
		if(sequenceId == debugIdx) System.out.println("Assembly graph. Initial edges end "+edgesE.size()+" Max score end: "+maxScoreE+" remaining edges: "+edgesMap.get(vE.getUniqueNumber()).size());
		
		double maxScore = Math.max(maxScoreS, maxScoreE);
		List<AssemblyEmbedded> embeddedList= new ArrayList<AssemblyEmbedded>();
		embeddedList.addAll(getEmbeddedBySequenceId(sequenceId));
		double maxE = 0;
		AssemblyEmbedded maxEmbedded = null;
		for(AssemblyEmbedded embedded:embeddedList) {
			double score = calculateScore(embedded);
			if(sequenceId == debugIdx) System.out.println("Assembly graph. Max score: "+maxScore+" embedded host: "+embedded.getHostId()+" embedded score: "+score+" max embedded score: "+maxE);
			if(score > 0.7*maxScore && score > maxE) {
				maxE = score;
				maxEmbedded = embedded;
			}
		}
		for(AssemblyEmbedded embedded:embeddedList) {
			if(embedded!=maxEmbedded) {
				removeEmbedded(embedded);
				if(sequenceId == debugIdx) System.out.println("Assembly graph. Removed embedded host: "+embedded.getHostId()+" Embedded relations: "+embeddedMapBySequence.get(sequenceId)+" is embedded: "+isEmbedded(sequenceId));
			} else {
				removeEdges(embedded.getSequenceId());
				if(sequenceId == debugIdx) System.out.println("Assembly graph. Removed edges for sequence: "+embedded.getSequenceId());
			}
		}
	}

	private void filterEdgesAbnormalOverlap(List<AssemblyEdge> edges) {
		List<AssemblyEdge> toRemove = new ArrayList<AssemblyEdge>();
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge()) continue;
			int overlap = edge.getOverlap();
			if(overlap>edge.getVertex1().getRead().getLength() || overlap>edge.getVertex2().getRead().getLength()) {
				toRemove.add(edge);
			}
		}
		for(AssemblyEdge edge:toRemove) removeEdge(edge);
	}

	private double calculateScore(AssemblyEmbedded embedded) {
		//KmerHitsCluster cluster = embedded.getEvidence();
		//return 1.0*(cluster.getQueryEvidenceEnd()-cluster.getQueryEvidenceStart())*cluster.getWeightedCount()/cluster.getQuery().length();
		return embedded.getCoverageSharedKmers();
	}

	private double calculateScore(AssemblyEdge edge) {
		//KmerHitsCluster cluster = edge.getEvidence();
		//return (cluster.getQueryEvidenceEnd()-cluster.getQueryEvidenceStart())*cluster.getWeightedCount()/(edge.getOverlap()+1);
		return edge.getCoverageSharedKmers();
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
					int rootStartParent = parentObject.getStartPosition();
					int rootStartSequence = rootStartParent+embedded.getStartPosition();
					boolean reverse = parentObject.isReverse()!=embedded.isReverse();
					embeddedSequencesMap.put(seqId, new AssemblyEmbedded(seqId, embedded.getRead(), reverse, sequenceIndex, rootStartSequence));
				}
				
				agenda.add(embedded.getSequenceId());
			}
					
		}
		List<AssemblyEmbedded> answer = new ArrayList<AssemblyEmbedded>();
		answer.addAll(embeddedSequencesMap.values());
		Collections.sort(answer, (a1,a2)-> a1.getStartPosition()-a2.getStartPosition());
		
		return answer;
	}

	public void filterEdgesAndEmbedded() {
		for (int seqId = sequences.size()-1; seqId >=0; seqId--) {
			filterEdgesAndEmbedded(seqId);
		}
		pruneEmbeddedSequences();	
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
		if(vertex.getSequenceIndex()==debugIdx) System.out.println("Filter edges close. Total edges: "+edges.size());
		if(edges.size()<3) return;
		AssemblyEdge maxOverlapEdge = null;
		for(AssemblyEdge edge: edges) {
			if(edge.isSameSequenceEdge()) continue;
			if(maxOverlapEdge==null || maxOverlapEdge.getOverlap() < edge.getOverlap()) {
				maxOverlapEdge = edge;
			}
		}
		AssemblyEdge secondOverlapEdge = null;
		for(AssemblyEdge edge: edges) {
			if(edge.isSameSequenceEdge() || edge == maxOverlapEdge) continue;
			if(secondOverlapEdge==null || secondOverlapEdge.getOverlap() < edge.getOverlap()) {
				secondOverlapEdge = edge;
			}
		}
		if(vertex.getSequenceIndex()==debugIdx) System.out.println("Filter edges close. Best edge: "+maxOverlapEdge.getConnectingVertex(vertex).getUniqueNumber()+" second: "+secondOverlapEdge.getConnectingVertex(vertex).getUniqueNumber());
		int ov1 = maxOverlapEdge.getOverlap();
		int ov2 = secondOverlapEdge.getOverlap();
		int diff = ov1-ov2;
		if(vertex.getSequenceIndex()==debugIdx) System.out.println("Filter edges close. overlap 1: "+ov1+" overlap2: "+ov2+" diff: "+diff);
		if(diff> 0.01*ov1) return;
		AssemblyVertex vMax = maxOverlapEdge.getConnectingVertex(vertex);
		AssemblyEdge seqEdgeMax = getSameSequenceEdge(vMax);
		AssemblyVertex vMax2 = seqEdgeMax.getConnectingVertex(vMax);
		AssemblyVertex vSecond = secondOverlapEdge.getConnectingVertex(vertex);
		AssemblyEdge seqEdgeSecond = getSameSequenceEdge(vSecond);
		AssemblyVertex vSecond2 = seqEdgeSecond.getConnectingVertex(vSecond);
		AssemblyEdge transitiveGood = getEdge(vMax2, vSecond);
		AssemblyEdge transitiveBad = getEdge(vSecond2, vMax);
		if(vertex.getSequenceIndex()==debugIdx) System.out.println("Filter edges close. transitive good: "+transitiveGood+" transitive bad: "+transitiveBad);
		if(transitiveBad!=null && transitiveGood==null) {
			//System.out.println("Filter edges close. removing edge: "+maxOverlapEdge.getVertex1().getUniqueNumber()+" "+maxOverlapEdge.getVertex2().getUniqueNumber()+" transitive good: "+transitiveGood+" transitive bad: "+transitiveBad);
			removeEdge(maxOverlapEdge);
		}
	}
	
}
