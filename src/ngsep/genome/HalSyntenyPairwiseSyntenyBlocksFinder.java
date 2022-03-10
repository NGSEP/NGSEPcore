package ngsep.genome;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;



public class HalSyntenyPairwiseSyntenyBlocksFinder implements PairwiseSyntenyBlocksFinder {
	
	private int minBlockLength = DEF_MIN_BLOCK_LENGTH;
	private int minHomologUnitsBlock = DEF_MIN_HOMOLOGY_UNITS_BLOCK;
	private int maxDistance = DEF_MAX_DISTANCE_BETWEEN_UNITS;
	
	public int getMinBlockLength() {
		return minBlockLength;
	}
	public void setMinBlockLength(int minBlockLength) {
		this.minBlockLength = minBlockLength;
	}
	
	public int getMinHomologUnitsBlock() {
		return minHomologUnitsBlock;
	}
	public void setMinHomologUnitsBlock(int minHomologUnitsBlock) {
		this.minHomologUnitsBlock = minHomologUnitsBlock;
	}

	public int getMaxDistance() {
		return maxDistance;
	}
	public void setMaxDistance(int maxDistance) {
		this.maxDistance = maxDistance;
	}

	@Override
	public List<PairwiseSyntenyBlock> findSyntenyBlocks(AnnotatedReferenceGenome g1, AnnotatedReferenceGenome g2, List<HomologyCluster> clusters) {
		List<PairwiseSyntenyBlock> syntenyBlocks = new ArrayList<PairwiseSyntenyBlock>();
		//System.out.println("Building graph. Clusters: "+clusters.size());
		List<SyntenyVertex> vertices = buildVertices(g1, g2, clusters);
		//System.out.println("Built graph with "+vertices.size()+" vertices");
		int n = vertices.size();
		int [] vertexWeights = new int [n];
		for (int i=0;i<n;i++) {
			vertexWeights[i] = vertices.get(i).getLocalRegion2().length();
		}
		Set<Integer> verticesInBlocks = new HashSet<>();
		List<List<HalSyntenyEdge>> edges = new ArrayList<>(n);
		for (int i=0;i<n;i++) {
			edges.add(buildEdgesVertex(vertices, vertexWeights, i));
		}
		System.out.println("Built graph with "+vertices.size()+" vertices and "+edges.size()+" edges");
		List<Integer> nextPathIndexes = findNextPath(vertexWeights,edges,verticesInBlocks);
		while (nextPathIndexes!=null) {
			//System.out.println("Found next path of size: "+nextPathIndexes.size());
			List<SyntenyVertex> path = new ArrayList<>(nextPathIndexes.size());
			for(int i:nextPathIndexes) {
				path.add(vertices.get(i));
				verticesInBlocks.add(i);
			}
			//Discard if it has too few units but add to vertices in blocks to avoid infinite loop
			if(nextPathIndexes.size()>minHomologUnitsBlock) {
				PairwiseSyntenyBlock block = new PairwiseSyntenyBlock(path);
				syntenyBlocks.add(block);
			}
			nextPathIndexes = findNextPath(vertexWeights,edges,verticesInBlocks);
		}
		return syntenyBlocks;
	}

	private List<SyntenyVertex> buildVertices(AnnotatedReferenceGenome g1, AnnotatedReferenceGenome g2, List<HomologyCluster> clusters) {
		List<SyntenyVertex>  vertices = new ArrayList<SyntenyVertex>(clusters.size());
		for(HomologyCluster cluster:clusters) {
			List<LocalHomologyCluster> clustersG1 = cluster.getLocalClusters(g1.getId());
			List<LocalHomologyCluster> clustersG2 = cluster.getLocalClusters(g2.getId());
			if(clustersG1.size()*clustersG2.size()>100) continue;
			for(LocalHomologyCluster c1:clustersG1) {
				for(LocalHomologyCluster c2: clustersG2) {
					SyntenyVertex v = new SyntenyVertex(c1, c2);
					vertices.add(v);
				}
			}	
		}
		GenomicRegionComparator cmp = new GenomicRegionComparator(g1.getSequencesMetadata());
		vertices.sort((a,b)->cmp.compare(a.getLocalRegion1(), b.getLocalRegion1()));
		//DEBUG
		for(int i=0;i<vertices.size();i++) {
			SyntenyVertex v = vertices.get(i);
			LocalHomologyCluster c1 = v.getLocalRegion1();
			LocalHomologyCluster c2 = v.getLocalRegion2();
			//if(c1.getFirst()>240000 && c1.getLast()<255000) System.out.println("Next vertex g1 "+c1.getSequenceName()+":"+c1.getFirst()+"-"+c1.getLast()+" g2 "+c2.getSequenceName()+":"+c2.getFirst()+"-"+c2.getLast()+" idx: "+i);
		}
		return vertices;
	}
	
	private List<HalSyntenyEdge> buildEdgesVertex(List<SyntenyVertex> vertices, int [] vertexWeights, int i) {
		List<HalSyntenyEdge> edges = new ArrayList<>();
		SyntenyVertex v1 = vertices.get(i);
		for (int j=i+1;j<vertices.size();j++) {
			SyntenyVertex v2 = vertices.get(j);
			if (!v1.getLocalRegion1().getSequenceName().equals(v2.getLocalRegion1().getSequenceName())) {
				break;
			} else if (v2.getLocalRegion1().getFirst()-v1.getLocalRegion1().getLast()> maxDistance) {
				break;
			}
			int d = calculatePossibleEdgeDistance(v1, v2);
			if(d>=0) {
				HalSyntenyEdge edge = new HalSyntenyEdge(i, j, d+vertexWeights[j], v2.getLocalRegion2().getFirst()>v1.getLocalRegion2().getFirst());
				edges.add(edge);
			} 
		}
		return edges;
	}
	
	private int calculatePossibleEdgeDistance(SyntenyVertex v1, SyntenyVertex v2) {
		//LocalHomologyCluster lc11 = v1.getLocalRegion1();
		//LocalHomologyCluster lc21 = v2.getLocalRegion1();
		LocalHomologyCluster lc12 = v1.getLocalRegion2();
		LocalHomologyCluster lc22 = v2.getLocalRegion2();
		if(!lc12.getSequenceName().equals(lc22.getSequenceName())) return -1;
		int d1 = lc22.getFirst()-lc12.getLast();
		int d2 = lc12.getFirst()-lc22.getLast();
		int d = d1>0?d1:d2;
		if(d>maxDistance) return -1;
		if(d<0) d=0;
		return d;
	}
	private List<Integer> findNextPath(int[] vertexWeights, List<List<HalSyntenyEdge>> edges, Set<Integer> verticesInBlocks) {
		int [] pathWeights = Arrays.copyOf(vertexWeights, vertexWeights.length);
		HalSyntenyEdge [] predecessors = new HalSyntenyEdge[vertexWeights.length];
		Arrays.fill(predecessors, null);
		int maxWeightIdx = -1;
		int maxWeight = 0;
		for (int i=0;i<vertexWeights.length;i++) {
			if(verticesInBlocks.contains(i)) continue;
			//if(i>=1027 && i<=1035) System.out.println("Vertex: "+i+" initialW: "+vertexWeights[i]+" pathW "+pathWeights[i]+" edges: "+edges.get(i).size());
			int weight = pathWeights[i];
			if(weight>maxWeight) {
				maxWeight = weight;
				maxWeightIdx = i;
			}
			HalSyntenyEdge pI = predecessors[i];
			List<HalSyntenyEdge> edgesV1= edges.get(i);
			for (HalSyntenyEdge se : edgesV1) {
				int j = se.v2Index;
				//if(j>=1027 && j<=1035) System.out.println("Next edge: "+se);
				if(pI!=null && pI.v2Positive!=se.v2Positive) continue;
				
				if(verticesInBlocks.contains(j)) continue;
				int newWeight = weight + se.weight;
				int currentWeight = pathWeights[j];
				if (newWeight > currentWeight) {
					pathWeights[j] = newWeight;
					predecessors[j] = se;
				}				
			}
		}
		if(maxWeight<minBlockLength) return null;
		//Traceback from maximum
		List<Integer> answer = new ArrayList<>();
		answer.add(maxWeightIdx);
		HalSyntenyEdge predEdge = predecessors[maxWeightIdx];
		while(predEdge!=null) {
			int i = predEdge.v1Index;
			answer.add(i);
			predEdge = predecessors[i];
		}
		Collections.reverse(answer);
		return answer;
	}
}
class HalSyntenyEdge {
	int v1Index;
	int v2Index;
	int weight;
	boolean v2Positive;
	public HalSyntenyEdge(int v1Index, int v2Index, int weight, boolean v2Positive) {
		super();
		this.v1Index = v1Index;
		this.v2Index = v2Index;
		this.weight = weight;
		this.v2Positive = v2Positive;
	}
	public String toString() {
		return " i1: "+v1Index+" i2: "+v2Index+" w: "+weight+" p: "+v2Positive;
	}
	
}
