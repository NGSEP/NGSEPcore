package ngsep.genome;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;


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
	public List<SyntenyBlock> findSyntenyBlocks(List<HomologyEdge> homologyEdges) {
		List<SyntenyBlock> syntenyBlocks = new ArrayList<>();
		List<SyntenyVertex> vertices = buildGraph(homologyEdges);
		HashMap<SyntenyVertex, List<SyntenyEdge>> orthologsTraceback = traverseVertices(vertices);
		SyntenyVertex maxVertex = findMaximalWeight(vertices);
		while (maxVertex!=null) {
			List<SyntenyEdge> path = orthologsTraceback.getOrDefault(maxVertex, new ArrayList<SyntenyEdge>());
			SyntenyBlock block = buildSyntenyBlock(path);
			if(block == null) break;
			//TODO: Check how length is validated
			int length = 0;
			for (SyntenyEdge se : path) length += se.getSource().getHomologyRelationship().getQueryUnit().length() + se.getTarget().getHomologyRelationship().getQueryUnit().length();
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
			SyntenyVertex sv = new SyntenyVertex(rel, rel.getSubjectUnit().length());
			vertices.add(sv);
		}
		for (int i=0;i<n;i++) {
			SyntenyVertex v1 = vertices.get(i);
			HomologyEdge rel1 = v1.getHomologyRelationship();
			HomologyUnit q1 = rel1.getQueryUnit();
			int centerQ1 = (q1.getFirst()+q1.getLast())/2;
			HomologyUnit s1 = rel1.getSubjectUnit();
			int centerS1 = (s1.getFirst()+s1.getLast())/2;
			for (int j=i+1;j<n;j++) {
				SyntenyVertex v2 = vertices.get(i);
				HomologyEdge rel2 = v2.getHomologyRelationship();
				HomologyUnit q2 = rel2.getQueryUnit();
				int centerQ2 = (q2.getFirst()+q2.getLast())/2;
				HomologyUnit s2 = rel2.getSubjectUnit();
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
	
	private HashMap<SyntenyVertex, List<SyntenyEdge>> traverseVertices(List<SyntenyVertex> vertices) {
		HashMap<SyntenyVertex, List<SyntenyEdge>> orthologsTraceback = new HashMap<SyntenyVertex, List<SyntenyEdge>>();
		for (SyntenyVertex vertex:vertices) {
			List<SyntenyEdge> edges = vertex.getEdges();
			for (SyntenyEdge se : edges) {
				SyntenyVertex target = se.getTarget();
				int newWeight = se.getWeight() + vertex.getWeight();
				int currentWeight = target.getWeight();
				if (newWeight > currentWeight) {
					target.setWeight(newWeight);
					List<SyntenyEdge> parentsSource = orthologsTraceback.getOrDefault(vertex, new ArrayList<SyntenyEdge>());
					List<SyntenyEdge> parentsTarget = new ArrayList<SyntenyEdge>();
					parentsTarget.addAll(parentsSource);
					parentsTarget.add(se);
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
	
	private SyntenyBlock buildSyntenyBlock(List<SyntenyEdge> path) {
		if (path.isEmpty()) return null;
		
		// remove vertices
		for (SyntenyEdge se : path) {
			SyntenyVertex s = se.getSource();
			SyntenyVertex t = se.getTarget();
			s.setWeight(-1);
			t.setWeight(-1);
		}
		SyntenyBlock sb = new SyntenyBlock(path); 
		
		return sb;

	}
	
	private List<SyntenyBlock> collapseOverlapped(List<SyntenyBlock> prev) {
		List<SyntenyBlock> temp = prev;
		while (selfOverlap(temp)) {
			List<SyntenyBlock> collapsed = new ArrayList<>();
			int i = 0;
			while (!temp.isEmpty()) {
				List<SyntenyBlock> toJoin = new ArrayList<>();
				SyntenyBlock actual = temp.get(i);
				toJoin.add(actual);
				for (int j = 0; j < temp.size(); j++) {
					SyntenyBlock other = temp.get(j);
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
	
	private SyntenyBlock joinCollapsed(List<SyntenyBlock> toJoin) {
		if (toJoin.size() == 1)
			return toJoin.get(0);
		List<SyntenyEdge> joinedHomologies = new ArrayList<>();
		for (SyntenyBlock e : toJoin) {
			joinedHomologies.addAll(e.getHomologies());
		}
		SyntenyBlock joined = new SyntenyBlock(joinedHomologies);
		return joined;
	}
	
	public boolean selfOverlap(List<SyntenyBlock> test) {
		for (int i = 0; i < test.size(); i++) {
			SyntenyBlock actual = test.get(i);
			for (int j = 0; j < test.size(); j++) {
				SyntenyBlock other = test.get(j);
				if (i != j && actual.getRegionGenome1().getSequenceName().equals(other.getRegionGenome1().getSequenceName()) && areOverlapping(actual, other)) {
					return true;
				}
			}
		}
		return false;
	}

	private boolean areOverlapping(SyntenyBlock e1, SyntenyBlock e2) {
		long start1 = e1.getRegionGenome1().getFirst();
		long end1 = e1.getRegionGenome1().getLast();
		long start2 = e2.getRegionGenome1().getFirst();
		long end2 = e2.getRegionGenome1().getLast();
		return (start1 >= start2 && start1 <= end2) || (end1 <= end2 && end1 >= start2)
				|| (start2 >= start1 && start2 <= end1) || (end2 <= end1 && end2 >= start1);
	}

	
}
