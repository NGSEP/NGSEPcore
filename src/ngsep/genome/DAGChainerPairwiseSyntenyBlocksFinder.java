package ngsep.genome;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class DAGChainerPairwiseSyntenyBlocksFinder  implements PairwiseSyntenyBlocksFinder {
	
	private int minBlockLength = DEF_MIN_BLOCK_LENGTH;
	private int maxDistance = DEF_MAX_DISTANCE_BETWEEN_UNITS;
	private int maxMatch;
	private int gapUnitLen;
	private int gapOpen;
	private int gapExtend;
	
	public DAGChainerPairwiseSyntenyBlocksFinder() {
		maxMatch = 100;
		gapUnitLen = 10000;
		gapOpen = 0;
		gapExtend = -3;
		maxDistance = 100000;
		minBlockLength = 6;
	}
	
	public int getMinBlockLength() {
		return minBlockLength;
	}
	public void setMinBlockLength(int minBlockLength) {
		this.minBlockLength = minBlockLength;
	}

	public int getMaxDistance() {
		return maxDistance;
	}
	
	public void setMaxDistance(int maxDistance) {
		this.maxDistance = maxDistance;
	}
	
	public int getMaxMatch() {
		return maxMatch;
	}
	
	public void setMaxMatch(int maxMatch) {
		this.maxMatch = maxMatch;
	}
	
	public int getGapUnitLen() {
		return gapUnitLen;
	}
	public void setGapUnitLen(int gapUnitLen) {
		this.gapUnitLen = gapUnitLen;
	}
	public int getGapOpen() {
		return gapOpen;
	}
	public void setGapOpen(int gapOpen) {
		this.gapOpen = gapOpen;
	}
	public int getGapExtend() {
		return gapExtend;
	}
	public void setGapExtend(int gapExtend) {
		this.gapExtend = gapExtend;
	}
	
	/**
	 * 
	 */
	public List<PairwiseSyntenyBlock> findSyntenyBlocks(AnnotatedReferenceGenome g1, AnnotatedReferenceGenome g2, List<HomologyCluster> clusters) {

		List<PairwiseSyntenyBlock> syntenyBlocks = new ArrayList<PairwiseSyntenyBlock>();
		List<SyntenyVertex> vertices = buildVertices(g1, g2, clusters);
		
		System.out.println("Built graph with " + vertices.size() +" vertices");
		
		//Reconstruct synteny blocks in the same orientation
		List<DAGChainerEdge> edges = buildEdgesVertex(vertices, false);
		
		System.out.println("Built graph with " + edges.size() +" edges");
	
		Set<Integer> verticesInBlocks = new HashSet<>();
		
		List<Integer> nextPathIndexes = findPaths(vertices,edges,verticesInBlocks);
		while (nextPathIndexes!=null) {
			//System.out.println("Found next path of size: "+nextPathIndexes.size());
			List<SyntenyVertex> path = new ArrayList<>(nextPathIndexes.size());
			for(int i:nextPathIndexes) {
				path.add(vertices.get(i));
				verticesInBlocks.add(i);
			}
			PairwiseSyntenyBlock block = new PairwiseSyntenyBlock(path);
			syntenyBlocks.add(block);
			
			nextPathIndexes = findPaths(vertices,edges,verticesInBlocks);
		}
		
		//Reconstruct synteny blocks in the opposite orientation
		edges = buildEdgesVertex(vertices, true);
		
		System.out.println("Built graph with " + edges.size() +" edges");
	
		verticesInBlocks = new HashSet<>();
		
		nextPathIndexes = findPaths(vertices,edges,verticesInBlocks);
		while (nextPathIndexes!=null) {
			//System.out.println("Found next path of size: "+nextPathIndexes.size());
			List<SyntenyVertex> path = new ArrayList<>(nextPathIndexes.size());
			for(int i:nextPathIndexes) {
				path.add(vertices.get(i));
				verticesInBlocks.add(i);
			}
			PairwiseSyntenyBlock block = new PairwiseSyntenyBlock(path);
			syntenyBlocks.add(block);

			nextPathIndexes = findPaths(vertices,edges,verticesInBlocks);
		}
		
		
		return syntenyBlocks;
	}
	
	/**
	 * Build vertices of the graph from Homology Clusters. One vertex per each pair of matches.
	 * @param g1 Annotated Reference Genome 1
	 * @param g2 Annotated Reference Genome 2
	 * @param clusters List of Homology Clusters
	 * @return List of Sorted Synteny Vertices
	 */
	private List<SyntenyVertex> buildVertices(AnnotatedReferenceGenome g1, AnnotatedReferenceGenome g2, List<HomologyCluster> clusters) {
		List<SyntenyVertex>  vertices = new ArrayList<SyntenyVertex>(clusters.size());
		for(HomologyCluster cluster:clusters) {
			List<LocalHomologyCluster> clustersG1 = cluster.getLocalClusters(g1.getId());
			List<LocalHomologyCluster> clustersG2 = cluster.getLocalClusters(g2.getId());
			
			for(LocalHomologyCluster c1:clustersG1) {
				for(LocalHomologyCluster c2: clustersG2) {
					SyntenyVertex v = new SyntenyVertex(c1, c2);
					vertices.add(v);
				}
			}	
		}
		GenomicRegionComparator cmp = new GenomicRegionComparator(g1.getSequencesMetadata());
		vertices.sort((a,b)->cmp.compare(a.getLocalRegion1(), b.getLocalRegion1()));
		
		return vertices;
	}
	
	/**
	 * Build edges of the graph if two vertices are in order within the genome
	 * @param vertices
	 * @return List of Edges
	 */
	private List<DAGChainerEdge> buildEdgesVertex(List<SyntenyVertex> vertices,boolean inverse) {
		List<DAGChainerEdge> edges = new ArrayList<>();
		
		for (int i=0; i<vertices.size();i++) {
			SyntenyVertex v1 = vertices.get(i);
			for (int j=i+1;j<vertices.size();j++) {
				SyntenyVertex v2 = vertices.get(j);
				LocalHomologyCluster v1g1 = v1.getLocalRegion1();
				LocalHomologyCluster v2g1 = v2.getLocalRegion1();
				LocalHomologyCluster v1g2 = v1.getLocalRegion2();
				LocalHomologyCluster v2g2 = v2.getLocalRegion2();
				
				if (!v1g1.getSequenceName().equals(v2g1.getSequenceName()) || !v1g2.getSequenceName().equals(v2g2.getSequenceName())) {
					continue;
				} 
				else if (!calculatePossibleEdge(v1, v2, inverse)) {
					continue;
				}
				else {
					DAGChainerEdge edge = new DAGChainerEdge(i, j);
					edges.add(edge);
				}		
			}
		}	
		return edges;
	}
	
	private boolean calculatePossibleEdge(SyntenyVertex v1, SyntenyVertex v2, boolean inverse) {
		
		int v1g1pos = calculateMidPositiion(v1.getLocalRegion1());
		int v1g2pos = calculateMidPositiion(v1.getLocalRegion2());
		int v2g1pos = calculateMidPositiion(v2.getLocalRegion1());
		int v2g2pos = calculateMidPositiion(v2.getLocalRegion2());
		
		
		if(!inverse && v1g1pos<v2g1pos && v1g2pos < v2g2pos)
			return true;
		else if(inverse && v1g1pos<v2g1pos && v1g2pos > v2g2pos)
			return true;
		return false;
	}
	
	/**
	 * 
	 * @param cluster
	 * @return
	 */
	public int calculateMidPositiion(LocalHomologyCluster cluster)
	{
		return (cluster.getLast() + cluster.getFirst())/2;
	}
	/**
	 * 
	 * @param vertices
	 * @param edges
	 * @return
	 */
	private List<Integer> findPaths(List<SyntenyVertex> vertices, List<DAGChainerEdge> edges, Set<Integer> verticesInBlocks)
	{
		double[] paths = new double[vertices.size()];
		int[] predecessor = new int[vertices.size()];
		Arrays.fill(predecessor, -1);
				
		HashMap<Integer, ArrayList<Integer>> mapPredecessor = buildMapPredecessors(edges);
		
		//Calculate scores using dynamic programming
		for (int i = 0; i < paths.length; i++) 
		{
			if(verticesInBlocks.contains(i)) continue;
			double maxScore = 0;
			int bestPredecessor = -1;
			
			ArrayList<Integer> listPredecessors = mapPredecessor.get(i);
			
			if(listPredecessors!=null){
				for(Integer vOrigen : listPredecessors){
					double scoreVertexOrigen = 	Math.max(paths[vOrigen] + calculateGapPenalty(vertices.get(vOrigen),vertices.get(i)), 0);
					if(scoreVertexOrigen>maxScore){
						maxScore = scoreVertexOrigen;
						bestPredecessor = vOrigen;
					}
				}
				predecessor[i] = bestPredecessor;
			}
			paths[i] = vertices.get(i).getMaximumEdgePCTSharedKmers() + maxScore;
			//System.out.println(i + " score: " + paths[i] + " contig: " + vertices.get(i).getLocalRegion1().getSequenceName() + " " + vertices.get(i).getLocalRegion2().getSequenceName() + " predecesor " + predecessor[i]);
		}
		
		List<Integer> homologies = extractBestPath(vertices,edges,paths,predecessor);
		
		return homologies;
	}
	
	
	private List<Integer> extractBestPath(List<SyntenyVertex> vertices, List<DAGChainerEdge> edges, double[] paths,
			int[] predecessor) {
		
		List<Integer> homologies = new ArrayList<>();
		int bestVertex = -1;
		double  maxScore = 0;
		
		//Find best path
		for (int i = 0; i < paths.length; i++) {
			if(paths[i]>maxScore){
				bestVertex = i;
				maxScore = paths[i];
			}
		}
		System.out.println("The best path score is  " + maxScore + " and ends in the vertex " + bestVertex);
		
		int countVertex = 1;
		int i = bestVertex;
		while(predecessor[i] >= 0){
			homologies.add(i);
			i = predecessor[i];
			countVertex ++;
		}
		
		Collections.reverse(homologies);
		
		//Calculate block Size (Distance or genes)
		//int end = calculateMidPositiion(vertices.get(bestVertex).getLocalRegion1());
		//int start = calculateMidPositiion(vertices.get(predecessor[i]).getLocalRegion1());		
		//if(end-start >= minBlockLength)
		if(countVertex>=minBlockLength)
			return homologies;
		
		return null;
	}

	/**
	 * Builds the map of possible predecessors of a SyntenyVertex
	 * @param edges
	 * @return
	 */
	private HashMap<Integer,ArrayList<Integer>> buildMapPredecessors(List<DAGChainerEdge> edges)
	{
		HashMap<Integer, ArrayList<Integer>> mapPredecessor = new HashMap<>();
		
		for (DAGChainerEdge edge : edges){
			if(mapPredecessor.get(edge.v2Index)==null){
				mapPredecessor.put(edge.v2Index, new ArrayList<Integer>());
			}
			mapPredecessor.get(edge.v2Index).add(edge.v1Index);
		}
		
		return mapPredecessor;
	}

	public double calculateNumGaps(int distanceG1,int distanceG2) {
		
		double numGaps = (distanceG1 + distanceG2 + Math.abs(distanceG1 - distanceG2))/(2 * gapUnitLen) + 0.5;
		return numGaps;
	}
	
	
	public double calculateGapPenalty(SyntenyVertex v1, SyntenyVertex v2) {
	
		int v1g1pos = calculateMidPositiion(v1.getLocalRegion1());
		int v1g2pos = calculateMidPositiion(v1.getLocalRegion2());
		int v2g1pos = calculateMidPositiion(v2.getLocalRegion1());
		int v2g2pos = calculateMidPositiion(v2.getLocalRegion2());
		
		double gapPenalty = 0;
		
		int distanceG1 = v2g1pos - v1g1pos;
		int distanceG2 = v2g2pos - v1g2pos;
		
		if(v2g2pos<v1g2pos)
			distanceG2 = v1g2pos - v2g2pos;
		
		double numGaps = calculateNumGaps(distanceG1, distanceG2);
		
		
		if (Math.max(distanceG1, distanceG2) > maxDistance) {
			return -1000000000;
		}
		else if(numGaps > 0) {
			gapPenalty = gapOpen + (numGaps * gapExtend);
		}
		
		return gapPenalty;
	}
}

class DAGChainerEdge {
	int v1Index;
	int v2Index;
	public DAGChainerEdge(int v1Index, int v2Index) {
		super();
		this.v1Index = v1Index;
		this.v2Index = v2Index;
	}
	public String toString() {
		return " i1: "+v1Index+" i2: "+v2Index;
	}
	
}
