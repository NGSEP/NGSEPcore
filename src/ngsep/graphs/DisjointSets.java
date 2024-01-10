package ngsep.graphs;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class DisjointSets {
	private int [] parents;
	private int [] h;
	private int numSubsets;
	public DisjointSets (int n) {
		parents= new int[n];
		h = new int[n];
		for(int i=0;i<n;i++) {
			parents[i] = i;
			h[i]=1;
		}
		numSubsets = n;
	}
	public void union (int i, int j) {
		int s1 = find(i);
		int s2 = find(j);
		if(h[s1]<h[s2]) parents[s1] = s2;
		else {
			parents[s2] = s1;
			if(h[s1]==h[s2]) h[s1]++;
		}
		numSubsets--;
	}
	public int find (int i) {
		int p = parents[i];
		if(p==i) return p;
		int r = find(p);
		parents[i]= r;
		return r;
	}
	
	public boolean sameSubsets(int i, int j) {
		return find(i)==find(j);
	}
	public int getNumSubsets() {
		return numSubsets;
	}
	public Map<Integer,Set<Integer>> getSubsets() {
		Map<Integer,Set<Integer>> subsetsMap = new HashMap<>();
		for(int i=0;i<parents.length;i++) {
			int s = find(i);
			Set<Integer> set = subsetsMap.computeIfAbsent(s,l->new HashSet<>());
			set.add(i);
		}
		return subsetsMap;
	}
	
}
