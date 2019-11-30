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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * @author Jorge Duitama
 * @author Juan Camilo Bojaca
 * @author David Guevara
 */
public class AssemblyGraph {
	
	private List<CharSequence> sequences;
	private List<AssemblyVertex> vertices;

	private boolean [] embedded;
	private List<AssemblyEdge> edges = new ArrayList<>();
	// Map with sequence ids as keys and embedded sequences as objects
	private Map<Integer, List<AssemblyEmbedded>> embeddedMap = new HashMap<>();

	// Indexes in the vertices list
	private List<List<AssemblyEdge>> paths = new ArrayList<List<AssemblyEdge>>();

	public AssemblyGraph(List<CharSequence> sequences) {
		int n = sequences.size();
		this.sequences = Collections.unmodifiableList(sequences);
		embedded = new boolean[n];
		Arrays.fill(embedded, false);
		vertices = new ArrayList<>(2*n);
		edges = new ArrayList<>(2*n);
		for (CharSequence seq : sequences) {
			AssemblyVertex vS = new AssemblyVertex(seq, true, vertices.size());
			vertices.add(vS);
			AssemblyVertex vE = new AssemblyVertex(seq, false, vertices.size());
			vertices.add(vE);
			addEdge(vS, vE, seq.length());
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

	public void addEdge(AssemblyVertex v1, AssemblyVertex v2, int overlap) {
		AssemblyEdge edge = new AssemblyEdge(v1, v2, overlap); 
		edges.add(edge);
		v1.addEdge(edge);
		v2.addEdge(edge);
	}

	public AssemblyVertex getVertex(int indexSequence, boolean start) {
		return vertices.get(2 * indexSequence + (start ? 0 : 1));
	}

	public void addEmbedded(int ind, AssemblyEmbedded embeddedObject) {
		List<AssemblyEmbedded> list = embeddedMap.computeIfAbsent(ind, key -> new LinkedList<>());
		list.add(embeddedObject);
		embedded[embeddedObject.getSequenceId()] = true;
	}
	

	/**
	 * Return the list of embedded sequences for the given read
	 * 
	 * @param index of the read
	 * @return list of embedded sequences
	 */
	public List<AssemblyEmbedded> getEmbedded(int index) {
		return embeddedMap.get(index);
	}
	
	public boolean isEmbedded(int sequenceId) {
		return embedded[sequenceId];
	}

	public void addPath(List<AssemblyEdge> path) {
		paths.add(path);
	}

	/**
	 * @return the vertices
	 */
	public List<AssemblyVertex> getAllVertices() {
		return vertices;
	}
	
	public List<AssemblyVertex> getNotEmbeddedVertices() {
		List<AssemblyVertex> vertices = new ArrayList<>();
		for(int i=0;i<embedded.length;i++) {
			if(!embedded[i]) {
				vertices.add(getVertex(i, true));
				vertices.add(getVertex(i, false));
			}
		}
		return vertices;
	}

	/**
	 * @return the edges
	 */
	public List<AssemblyEdge> getEdges() {
		return edges;
	}

	/**
	 * @return the paths
	 */
	public List<List<AssemblyEdge>> getPaths() {
		return paths;
	}


	public void serialize(String outFileGraph) {
		// TODO : Implement
		
	}
	
	public static AssemblyGraph load(String filename) throws IOException {
		//TODO: Implement
		return null;
	}

	
}
