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
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import ngsep.math.Distribution;

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
	private List<CharSequence> sequences;
	/**
	 * Start vertices indexed by sequence index
	 */
	private Map<Integer,AssemblyVertex> verticesStart;
	/**
	 * End vertices indexed by sequence index
	 */
	private Map<Integer,AssemblyVertex> verticesEnd;

	private boolean [] embedded;
	private Map<AssemblyVertex,List<AssemblyEdge>> edgesMap = new HashMap<>();
	// Map with sequence ids as keys and embedded sequences as objects
	private Map<Integer, List<AssemblyEmbedded>> embeddedMap = new HashMap<>();

	private List<List<AssemblyEdge>> paths = new ArrayList<List<AssemblyEdge>>();
	
	/**
	 * Optional attribute to store the sequence names. Useful for the gold standard graph
	 */
	private List<String> readNames;

	public AssemblyGraph(List<CharSequence> sequences) {
		int n = sequences.size();
		this.sequences = Collections.unmodifiableList(sequences);
		embedded = new boolean[n];
		Arrays.fill(embedded, false);
		verticesStart = new HashMap<>(n);
		verticesEnd = new HashMap<>(n);
		edgesMap = new HashMap<>(n);
		for (int i=0;i<sequences.size();i++) {
			CharSequence seq = sequences.get(i);
			AssemblyVertex vS = new AssemblyVertex(seq, true, i);
			verticesStart.put(i,vS);
			edgesMap.put(vS, new ArrayList<>());
			AssemblyVertex vE = new AssemblyVertex(seq, false, i);
			verticesEnd.put(i,vE);
			edgesMap.put(vE, new ArrayList<>());
			addEdge(vS, vE, seq.length(), seq.length());
			embeddedMap.put(i, new ArrayList<>());
		}
	}

	/**
	 * @return the sequences
	 */
	public List<CharSequence> getSequences() {
		return sequences;
	}
	public CharSequence getSequence(int sequenceIdx) {
		return sequences.get(sequenceIdx);
	}
	public int getSequenceLength(int sequenceIdx) {
		return sequences.get(sequenceIdx).length();
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

	public void addEdge(AssemblyVertex v1, AssemblyVertex v2, int cost, int overlap) {
		AssemblyEdge edge = new AssemblyEdge(v1, v2, cost, overlap);
		edgesMap.get(v1).add(edge);
		edgesMap.get(v2).add(edge);
	}
	
	public void removeEdge (AssemblyEdge edge) {
		edgesMap.get(edge.getVertex1()).remove(edge);
		edgesMap.get(edge.getVertex2()).remove(edge);
	}

	public AssemblyVertex getVertex(int indexSequence, boolean start) {
		if(start) return verticesStart.get(indexSequence);
		return verticesEnd.get(indexSequence);
	}

	public void addEmbedded(int ind, AssemblyEmbedded embeddedObject) {
		List<AssemblyEmbedded> list = embeddedMap.computeIfAbsent(ind, key -> new LinkedList<>());
		list.add(embeddedObject);
		embedded[embeddedObject.getSequenceId()] = true;
	}
	
	public void pruneEmbeddedSequences() {
		for(int i=0;i<embedded.length;i++) {
			if(embedded[i] && verticesStart.get(i)!=null) {
				removeVertices(i);
			}
		}
	}

	public void removeVertices(int sequenceId) {
		AssemblyVertex v1 = getVertex(sequenceId, true);
		List<AssemblyEdge> edgesToRemove = new ArrayList<>(); 
		edgesToRemove.addAll(edgesMap.get(v1));
		AssemblyVertex v2 = getVertex(sequenceId, false);
		edgesToRemove.addAll(edgesMap.get(v2));
		for(AssemblyEdge edge:edgesToRemove) {
			removeEdge(edge);
		}
		edgesMap.remove(v1);
		edgesMap.remove(v2);
		verticesStart.remove(sequenceId);
		verticesEnd.remove(sequenceId);
	}
	

	/**
	 * Return the list of embedded sequences for the given read
	 * 
	 * @param index of the read
	 * @return list of embedded sequences
	 */
	public List<AssemblyEmbedded> getEmbedded(int index) {
		List<AssemblyEmbedded> answer = embeddedMap.get(index);
		if(answer!=null) return answer;
		return new ArrayList<AssemblyEmbedded>();
	}
	
	public boolean isEmbedded(int sequenceId) {
		return embedded[sequenceId];
	}
	public int getEmbeddedCount () {
		int count = 0;
		for(int i=0;i<embedded.length;i++) {
			if(embedded[i]) count++;
		}
		return count;
	}

	public void addPath(List<AssemblyEdge> path) {
		paths.add(path);
	}
	
	public List<AssemblyVertex> getVertices() {
		List<AssemblyVertex> vertices = new ArrayList<>();
		vertices.addAll(edgesMap.keySet());
		return vertices;
	}

	/**
	 * @return the edges
	 */
	public List<AssemblyEdge> getEdges() {
		List<AssemblyEdge> edges = new ArrayList<>();
		for(AssemblyVertex v:edgesMap.keySet()) {
			List<AssemblyEdge> edgesVertex = edgesMap.get(v);
			for(AssemblyEdge edge:edgesVertex) {
				//Avoid adding twice the same edge
				if(edge.getVertex1()==v) edges.add(edge);
			}
		}
		return edges;
	}
	public List<AssemblyEdge> getEdges(AssemblyVertex vertex) {
		return edgesMap.get(vertex);
	}
	/**
	 * Returns the edge connecting the given vertex with the corresponding vertex in the same sequence
	 * @param vertex
	 * @return AssemblyEdge
	 */
	public AssemblyEdge getSameSequenceEdge(AssemblyVertex vertex) {
		List<AssemblyEdge> edges = edgesMap.get(vertex);
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
		List<AssemblyEdge> edgesV1 = edgesMap.get(v1);
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
		Distribution answer = new Distribution(0, edgesMap.size(), 1);
		for(List<AssemblyEdge> edges:edgesMap.values()) {
			answer.processDatapoint(edges.size());
		}
		return answer;
	}

	

	

	

	
}
