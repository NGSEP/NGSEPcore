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

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPOutputStream;

import ngsep.main.io.ConcatGZIPInputStream;
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
			edge.setMismatches(0);
			edge.setCoverageSharedKmers(seq.getLength());
			edge.setNumSharedKmers(seq.getLength());
			edge.setOverlapStandardDeviation(0);
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
		int left = embedded.getHostStart();
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
		List<AssemblyEdge> toRemove = new ArrayList<AssemblyEdge>();
		List<AssemblyEdge> edges1 = edgesMap.get(v1.getUniqueNumber()); 
		if(edges1!=null) toRemove.addAll(edges1);
		AssemblyVertex v2 = getVertex(sequenceId, false);
		List<AssemblyEdge> edges2 = edgesMap.get(v2.getUniqueNumber());
		if(edges2!=null) {
			for(AssemblyEdge edge:edges2) {
				if(!edge.isSameSequenceEdge()) toRemove.add(edge);
			}
		}
		for(AssemblyEdge edge:toRemove) {
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
			if(edgesVertex==null) continue;
			for(AssemblyEdge edge:edgesVertex) {
				if(edge==null) continue;
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


	public void save(String outFileGraph) throws IOException {
		try (OutputStream os = new GZIPOutputStream(new FileOutputStream(outFileGraph));
			 PrintStream out = new PrintStream(os)) {
			out.println("#SEQUENCES");
			for(QualifiedSequence seq:sequences) {
				out.println(seq.getName()+"\t"+seq.getLength());
			}
			out.println("#EMBEDDED");
			for(List<AssemblyEmbedded> embeddedList:embeddedMapBySequence.values()) {
				for(AssemblyEmbedded embedded:embeddedList) {
					int reverse = embedded.isReverse()?1:0;
					out.print(""+embedded.getSequenceId()+"\t"+embedded.getHostId()+"\t"+embedded.getHostStart()+"\t"+embedded.getHostEnd()+"\t"+reverse);
					out.print("\t"+embedded.getNumSharedKmers()+"\t"+embedded.getCoverageSharedKmers()+"\t"+embedded.getWeightedCoverageSharedKmers());
					out.println("\t"+embedded.getMismatches()+"\t"+embedded.getHostEvidenceStart()+"\t"+embedded.getHostEvidenceEnd());
				}
			}
			
			out.println("#EDGES");
			List<AssemblyEdge> edges = getEdges();
			for(AssemblyEdge edge:edges) {
				if(edge.isSameSequenceEdge()) continue;
				out.print(""+edge.getVertex1().getUniqueNumber()+"\t"+edge.getVertex2().getUniqueNumber()+"\t"+edge.getOverlap());
				out.print("\t"+edge.getOverlapStandardDeviation()+"\t"+edge.getNumSharedKmers()+"\t"+edge.getCoverageSharedKmers());
				out.println("\t"+edge.getWeightedCoverageSharedKmers()+"\t"+edge.getMismatches());
			}
			
		}
	}
	
	public static AssemblyGraph load(List<QualifiedSequence> sequences, String graphFilename) throws IOException {
		AssemblyGraph graph = new AssemblyGraph(sequences);
		String line = null;
		try (ConcatGZIPInputStream gzs = new ConcatGZIPInputStream(new FileInputStream(graphFilename));
			 BufferedReader in = new BufferedReader(new InputStreamReader(gzs))) {
			line = in.readLine();
			if(!"#SEQUENCES".equals(line)) throw new IOException("Graph file misses sequence names. First line: "+line);
			int seqId = 0;
			line=in.readLine();
			while(line!=null && seqId<sequences.size() && !line.startsWith("#")) {
				String [] items = line.split("\t");
				QualifiedSequence seq = sequences.get(seqId);
				if(!seq.getName().equals(items[0]))  throw new IOException("Unexpected sequence name at index " +(seqId+2)+". Double check that the graph was built from the given sequences or build again the graph. Expected: "+seq.getName()+" "+seq.getLength()+" loaded: "+line);
				if(seq.getLength()!=Integer.parseInt(items[1])) throw new IOException("Unexpected sequence length at index" +(seqId+2)+". Sequence name: "+seq.getName()+". Double check that the graph was built from the given sequences or build again the graph. Expected: "+seq.getLength()+" loaded: "+items[1]);
				seqId++;
				line=in.readLine();
			}
			if(line==null) throw new IOException("Unexpected end of file reading sequences");
			//if(seqId<sequences.size()) log.info("Missing sequences in graph. Expected: "+sequences.size()+". Loaded: "+seqId);
			if (!line.equals("#EMBEDDED")) throw new IOException("Unexpected line after loading sequences. Expected: #EMBEDDED. Line: "+line);
			line=in.readLine();
			while(line!=null && !line.startsWith("#")) {
				String [] items = line.split("\t");
				int embSeqId = Integer.parseInt(items[0]);
				int hostId = Integer.parseInt(items[1]);
				boolean reverse = Integer.parseInt(items[4])==1;
				QualifiedSequence embeddedSeq = sequences.get(embSeqId);
				AssemblyEmbedded embedded = new AssemblyEmbedded(embSeqId, embeddedSeq, reverse, hostId, Integer.parseInt(items[2]), Integer.parseInt(items[3]));
				embedded.setNumSharedKmers(Integer.parseInt(items[5]));
				embedded.setCoverageSharedKmers(Integer.parseInt(items[6]));
				embedded.setWeightedCoverageSharedKmers(Integer.parseInt(items[7]));
				embedded.setMismatches(Integer.parseInt(items[8]));
				embedded.setHostEvidenceStart(Integer.parseInt(items[9]));
				embedded.setHostEvidenceEnd(Integer.parseInt(items[10]));
				graph.addEmbedded(embedded);
				line=in.readLine();
			}
			if(line==null) throw new IOException("Unexpected end of file reading embedded relationships");
			if (!line.equals("#EDGES")) throw new IOException("Unexpected line after loading embedded. Expected: #EDGES. Line: "+line);
			line=in.readLine();
			while(line!=null && !line.startsWith("#")) {
				String [] items = line.split("\t");
				int v1Idx = Integer.parseInt(items[0]);
				int v2Idx = Integer.parseInt(items[1]);
				int overlap = Integer.parseInt(items[2]);
				AssemblyVertex v1 = graph.getVertexByUniqueId(v1Idx);
				AssemblyVertex v2 = graph.getVertexByUniqueId(v2Idx);
				AssemblyEdge edge = new AssemblyEdge(v1, v2, overlap);
				edge.setOverlapStandardDeviation(Double.parseDouble(items[3]));
				edge.setNumSharedKmers(Integer.parseInt(items[4]));
				edge.setCoverageSharedKmers(Integer.parseInt(items[5]));
				edge.setWeightedCoverageSharedKmers(Integer.parseInt(items[6]));
				edge.setMismatches(Integer.parseInt(items[7]));
				graph.addEdge(edge);
				line=in.readLine();
			}
		} catch (NumberFormatException e) {
			throw new IOException("Error loading number at line: "+line,e);
		}
		graph.updateVertexDegrees();
		return graph;
	}
	
	public static List<QualifiedSequence> loadSequenceNamesFromGraphFile(String graphFilename) throws IOException {
		List<QualifiedSequence> sequenceNames = new ArrayList<QualifiedSequence>();
		String line = null;
		try (ConcatGZIPInputStream gzs = new ConcatGZIPInputStream(new FileInputStream(graphFilename));
			 BufferedReader in = new BufferedReader(new InputStreamReader(gzs))) {
			line = in.readLine();
			if(!"#SEQUENCES".equals(line)) throw new IOException("Graph file misses sequence names. First line: "+line);
			line=in.readLine();
			while(line!=null && !line.startsWith("#")) {
				String [] items = line.split("\t");
				QualifiedSequence seq = new QualifiedSequence(items[0]);
				seq.setLength(Integer.parseInt(items[1]));
				sequenceNames.add(seq);
				line=in.readLine();
			}
		} catch (NumberFormatException e) {
			throw new IOException("Error loading number at line: "+line,e);
		}
		return sequenceNames;
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

	public void removeEdgesChimericReads () {
		for(int i=0;i<sequences.size();i++) {
			if(isChimeric(i)) {
				removeEdges(i);
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
		if(vS==null || vE==null) return;
		filterEdgesAbnormalFeatures(getEdges(vS));
		filterEdgesAbnormalFeatures(getEdges(vE));
		List<AssemblyEdge> edgesS = new ArrayList<AssemblyEdge>();
		if(vS!=null) edgesS.addAll(getEdges(vS));
		double minScoreProportion = 0.4;
		//double minScoreProportion = 0.5;
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
			if(sequenceId == debugIdx) System.out.println("Assembly graph. Length "+getSequenceLength(sequenceId)+" coverageKmers: "+embedded.getCoverageSharedKmers()+" Max score: "+maxScore+" embedded host: "+embedded.getHostId()+" embedded score: "+score+" max embedded score: "+maxE);
			if(score > minScoreProportion*maxScore && score > maxE) {
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

	private void filterEdgesAbnormalFeatures(List<AssemblyEdge> edges) {
		List<AssemblyEdge> toRemove = new ArrayList<AssemblyEdge>();
		for(AssemblyEdge edge:edges) {
			if(edge.isSameSequenceEdge()) continue;
			int overlap = edge.getOverlap();
			int v1L = edge.getVertex1().getRead().getLength();
			int v2L = edge.getVertex2().getRead().getLength();
			if(overlap>1.1*v1L || overlap>1.1*v2L) {
				toRemove.add(edge);
			} else if (edge.getWeightedCoverageSharedKmers()<0.1*v1L || edge.getWeightedCoverageSharedKmers()<0.1*v2L) {
				toRemove.add(edge);
			}
		}
		for(AssemblyEdge edge:toRemove) removeEdge(edge);
	}

	private double calculateScore(AssemblyEmbedded embedded) {
		//KmerHitsCluster cluster = embedded.getEvidence();
		//return 1.0*(cluster.getQueryEvidenceEnd()-cluster.getQueryEvidenceStart())*cluster.getWeightedCount()/cluster.getQuery().length();
		//return embedded.getCoverageSharedKmers();
		return embedded.getWeightedCoverageSharedKmers();
	}

	private double calculateScore(AssemblyEdge edge) {
		//KmerHitsCluster cluster = edge.getEvidence();
		//return (cluster.getQueryEvidenceEnd()-cluster.getQueryEvidenceStart())*cluster.getWeightedCount()/(edge.getOverlap()+1);
		//return edge.getCoverageSharedKmers();
		return edge.getWeightedCoverageSharedKmers();
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

	public void filterEdgesAndEmbedded() {
		for (int seqId = sequences.size()-1; seqId >=0; seqId--) {
			filterEdgesAndEmbedded(seqId);
		}
		pruneEmbeddedSequences();
		filterEdgesCloseRelationships();
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
		AssemblyEdge maxScoreEdge = null;
		for(AssemblyEdge edge: edges) {
			if(edge.isSameSequenceEdge()) continue;
			if(maxOverlapEdge==null || maxOverlapEdge.getOverlap() < edge.getOverlap()) {
				maxOverlapEdge = edge;
			}
			if(maxScoreEdge==null || calculateScore(maxScoreEdge)<calculateScore(edge)) {
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
			removeEdge(maxOverlapEdge);
		} else if (maxOverlapEdge!=maxScoreEdge && badTransitivity(vertex,maxScoreEdge,maxOverlapEdge)) {
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
