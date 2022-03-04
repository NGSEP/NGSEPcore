package ngsep.genome;

import java.util.ArrayList;
import java.util.List;


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
		maxDistance = 10000;
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
		
		List<DAGChainerEdge> edges = buildEdgesVertex(vertices);
		
		System.out.println("Built graph with " + edges.size() +" edges");
	
		syntenyBlocks.addAll(findPaths(vertices, edges));
		
		
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
	private List<DAGChainerEdge> buildEdgesVertex(List<SyntenyVertex> vertices) {
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
					break;
				} 
				else if (!calculatePossibleEdge(v1, v2)) {
					break;
				}
				else {
						DAGChainerEdge edge = new DAGChainerEdge(i, j);
						edges.add(edge);
				}
							
			}
		}	
		return edges;
	}
	
	private boolean calculatePossibleEdge(SyntenyVertex v1, SyntenyVertex v2) {
		
		double v1g1pos = calculateMidPositiion(v1.getLocalRegion1());
		double v1g2pos = calculateMidPositiion(v1.getLocalRegion2());
		double v2g1pos = calculateMidPositiion(v2.getLocalRegion1());
		double v2g2pos = calculateMidPositiion(v2.getLocalRegion2());
		
		if(v1g1pos<v2g1pos && v1g2pos < v2g2pos)
			return true;
		
		return false;
	}
	
	/**
	 * 
	 * @param cluster
	 * @return
	 */
	public double calculateMidPositiion(LocalHomologyCluster cluster)
	{
		return (cluster.getLast()-cluster.getFirst())/2;
	}
	/**
	 * 
	 * @param vertices
	 * @param edges
	 * @return
	 */
	public List<PairwiseSyntenyBlock> findPaths(List<SyntenyVertex> vertices, List<DAGChainerEdge> edges)
	{
		List<PairwiseSyntenyBlock> syntenyBlocks = new ArrayList<PairwiseSyntenyBlock>();
		
		double[] paths = new double[vertices.size()];
		int[] predecesor = new int[vertices.size()];
		paths[0] = vertices.get(0).getMaximumEdgePCTSharedKmers();
		predecesor[0] = -1;
		
		for (int i = 1; i < paths.length; i++) 
		{
			predecesor[i]=-1;
			paths[i] = vertices.get(i).getMaximumEdgePCTSharedKmers();

			if(vertices.get(i).getLocalRegion1().getFirst() != vertices.get(i-1).getLocalRegion1().getFirst() && vertices.get(i).getLocalRegion2().getFirst() != vertices.get(i-1).getLocalRegion2().getFirst())
			{
				if(vertices.get(i).getLocalRegion1().getSequenceName().equals(vertices.get(i-1).getLocalRegion1().getSequenceName()))
				{
					paths[i] += Math.max(paths[i-1] + calculateGapPenalty(vertices.get(i-1),vertices.get(i)), 0);
					predecesor[i] = i-1;
				}
			}
			
			System.out.println(i + " score: " + paths[i] + " contig: " + vertices.get(i).getLocalRegion1().getSequenceName() + " predecesor " + predecesor[i]);
		}
		
		//Add blocks
		
		return syntenyBlocks;

	}

	public double calculateNumGaps(double v1g1pos,double v1g2pos,double v2g1pos,double v2g2pos) {
		
		double numGaps = ((v2g1pos-v1g1pos) + (v2g2pos-v1g2pos) + Math.abs((v2g1pos-v1g1pos) + (v2g2pos-v1g2pos)))/(2*gapUnitLen) + 0.5;
		return numGaps;
	}
	
	
	public double calculateGapPenalty(SyntenyVertex v1, SyntenyVertex v2) {
	
		double v1g1pos = calculateMidPositiion(v1.getLocalRegion1());
		double v1g2pos = calculateMidPositiion(v1.getLocalRegion2());
		double v2g1pos = calculateMidPositiion(v2.getLocalRegion1());
		double v2g2pos = calculateMidPositiion(v2.getLocalRegion2());
		
		double gapPenalty = 0;
		
		double numGaps = calculateNumGaps(v1g1pos,v1g2pos,v2g1pos,v2g2pos);
		
		if (Math.max(v2g1pos-v1g1pos, v2g2pos-v1g2pos) > maxDistance) {
			return -1000000000;
		}
		else if(numGaps > 0) {
			gapPenalty = gapOpen + (numGaps*gapExtend);
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
