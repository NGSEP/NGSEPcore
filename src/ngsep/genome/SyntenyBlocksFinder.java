package ngsep.genome;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class SyntenyBlocksFinder {
	
	private int minBlockLength;
	
	private int maxDistance;
	
	public SyntenyBlocksFinder(int minBlockLength, int maxDistance) {
		this.minBlockLength = minBlockLength;
		this.maxDistance = maxDistance;
	}
	
	/**
	 * PRE: the query unit of all homology edges belong to the same sequence within the same genome
	 * and the subject units belong to the same genome but can belong to different sequences
	 * @param homologyEdges
	 * @return List<SyntenyBlock> Blocks of synteny between the query genome and the subject genome
	 */
	public List<PairwiseSyntenyBlock> findSyntenyBlocks(List<HomologyEdge> homologyEdges) {
		List<PairwiseSyntenyBlock> syntenyBlocks = new ArrayList<>();
		List<SyntenyVertex> vertices = buildGraph(homologyEdges);
		Map<SyntenyVertex, List<SyntenyVertex>> orthologsTraceback = traverseVertices(vertices);
		SyntenyVertex maxVertex = findMaximalWeight(vertices);
		while (maxVertex!=null) {
			List<SyntenyVertex> path = orthologsTraceback.getOrDefault(maxVertex, new ArrayList<SyntenyVertex>());
			PairwiseSyntenyBlock block = buildSyntenyBlock(path);
			if(block == null) break;
			//TODO: Check how length is validated
			int length = 0;
			for (SyntenyVertex v : path) length += v.getHomologyCluster().getHomologyUnitsCluster().get(0).length();
			if (length >= minBlockLength) syntenyBlocks.add(block);
			maxVertex = findMaximalWeight(vertices);
		}
		return syntenyBlocks;
	}
	
	
	
	private List<SyntenyVertex> buildGraph(List<HomologyEdge> homologyRelationships) {
		int n = homologyRelationships.size();
		List<SyntenyVertex>  vertices = new ArrayList<SyntenyVertex>(n);
		homologyRelationships.sort((a,b)->a.getQueryUnit().getFirst()-b.getQueryUnit().getFirst());
		for (int i=0;i<n;i++) {
			HomologyEdge rel = homologyRelationships.get(i);
			List<HomologyUnit> units = new ArrayList<HomologyUnit>();
			units.add(rel.getQueryUnit());
			units.add(rel.getSubjectUnit());
			HomologyCluster cluster = new HomologyCluster(i, units);
			SyntenyVertex sv = new SyntenyVertex(cluster);
			sv.setWeight(rel.getSubjectUnit().length());
			vertices.add(sv);
		}
		for (int i=0;i<n;i++) {
			SyntenyVertex v1 = vertices.get(i);
			HomologyCluster c1 = v1.getHomologyCluster();
			HomologyUnit q1 = c1.getHomologyUnitsCluster().get(0);
			int centerQ1 = (q1.getFirst()+q1.getLast())/2;
			HomologyUnit s1 = c1.getHomologyUnitsCluster().get(1);
			int centerS1 = (s1.getFirst()+s1.getLast())/2;
			for (int j=i+1;j<n;j++) {
				SyntenyVertex v2 = vertices.get(i);
				HomologyCluster c2 = v2.getHomologyCluster();
				HomologyUnit q2 = c2.getHomologyUnitsCluster().get(0);
				int centerQ2 = (q2.getFirst()+q2.getLast())/2;
				HomologyUnit s2 = c2.getHomologyUnitsCluster().get(1);
				int centerS2 = (q2.getFirst()+q2.getLast())/2;
				int dq = centerQ2-centerQ1;
				if(dq>maxDistance) break;
				int ds = Math.abs(centerS1-centerS2);
				if(!s1.getSequenceName().equals(s2.getSequenceName()) || ds>maxDistance) continue; 
				//TODO: Check weight
				SyntenyEdge se = new SyntenyEdge(v1, v2, s2.length());
				v1.addEdge(se);
			}
		}
		return vertices;
	}
	
	private Map<SyntenyVertex, List<SyntenyVertex>> traverseVertices(List<SyntenyVertex> vertices) {
		Map<SyntenyVertex, List<SyntenyVertex>> orthologsTraceback = new HashMap<SyntenyVertex, List<SyntenyVertex>>();
		for (SyntenyVertex vertex:vertices) {
			List<SyntenyEdge> edges = vertex.getEdges();
			for (SyntenyEdge se : edges) {
				SyntenyVertex target = se.getTarget();
				int newWeight = se.getWeight() + vertex.getWeight();
				int currentWeight = target.getWeight();
				if (newWeight > currentWeight) {
					target.setWeight(newWeight);
					List<SyntenyVertex> parentsSource = orthologsTraceback.getOrDefault(vertex, new ArrayList<SyntenyVertex>());
					List<SyntenyVertex> parentsTarget = new ArrayList<SyntenyVertex>();
					parentsTarget.addAll(parentsSource);
					parentsTarget.add(vertex);
					orthologsTraceback.put(target, parentsTarget);
				}				
			}
		}
		return orthologsTraceback;
	}
	
	private SyntenyVertex findMaximalWeight(List<SyntenyVertex> vertices) {
		SyntenyVertex maxVertex = null;
		int maxWeight = 0;
		for (SyntenyVertex vertex:vertices) {
			int weight = vertex.getWeight();
			if (weight > maxWeight) {
				maxVertex = vertex;
				maxWeight = weight;
			}
		}
		return maxVertex;
	}
	
	private PairwiseSyntenyBlock buildSyntenyBlock(List<SyntenyVertex> path) {
		if (path.isEmpty()) return null;
		
		// remove vertices
		int genomeId1 = -1;
		int genomeId2 = -1;
		for (SyntenyVertex v : path) {
			v.setWeight(-1);
			if(genomeId1<0) {
				List<HomologyUnit> units = v.getHomologyCluster().getHomologyUnitsCluster(); 
				genomeId1 = units.get(0).getGenomeId();
				genomeId2 = units.get(1).getGenomeId();
			}
		}
		PairwiseSyntenyBlock sb = new PairwiseSyntenyBlock(genomeId1,genomeId2,path); 
		
		return sb;

	}
	
	/*private List<PairwiseSyntenyBlock> collapseOverlapped(List<PairwiseSyntenyBlock> prev) {
		List<PairwiseSyntenyBlock> temp = prev;
		while (selfOverlap(temp)) {
			List<PairwiseSyntenyBlock> collapsed = new ArrayList<>();
			int i = 0;
			while (!temp.isEmpty()) {
				List<PairwiseSyntenyBlock> toJoin = new ArrayList<>();
				PairwiseSyntenyBlock actual = temp.get(i);
				toJoin.add(actual);
				for (int j = 0; j < temp.size(); j++) {
					PairwiseSyntenyBlock other = temp.get(j);
					if (i != j && actual.getRegionGenome1().getSequenceName().equals(other.getRegionGenome1().getSequenceName()) && areOverlapping(actual, other)) {
						toJoin.add(other);
					}
				}
				temp.removeAll(toJoin);
				collapsed.add(joinCollapsed(toJoin));
			}
			temp = collapsed;
		}
		return temp;
	}
	
	private PairwiseSyntenyBlock joinCollapsed(List<PairwiseSyntenyBlock> toJoin) {
		if (toJoin.size() == 1)
			return toJoin.get(0);
		List<SyntenyEdge> joinedHomologies = new ArrayList<>();
		for (PairwiseSyntenyBlock e : toJoin) {
			joinedHomologies.addAll(e.getHomologies());
		}
		PairwiseSyntenyBlock joined = new PairwiseSyntenyBlock(joinedHomologies);
		return joined;
	}*/
	
	public boolean selfOverlap(List<PairwiseSyntenyBlock> test) {
		for (int i = 0; i < test.size(); i++) {
			PairwiseSyntenyBlock actual = test.get(i);
			for (int j = 0; j < test.size(); j++) {
				PairwiseSyntenyBlock other = test.get(j);
				if (i != j && actual.getRegionGenome1().getSequenceName().equals(other.getRegionGenome1().getSequenceName()) && areOverlapping(actual, other)) {
					return true;
				}
			}
		}
		return false;
	}

	private boolean areOverlapping(PairwiseSyntenyBlock e1, PairwiseSyntenyBlock e2) {
		long start1 = e1.getRegionGenome1().getFirst();
		long end1 = e1.getRegionGenome1().getLast();
		long start2 = e2.getRegionGenome1().getFirst();
		long end2 = e2.getRegionGenome1().getLast();
		return (start1 >= start2 && start1 <= end2) || (end1 <= end2 && end1 >= start2)
				|| (start2 >= start1 && start2 <= end1) || (end2 <= end1 && end2 >= start1);
	}

	
}
