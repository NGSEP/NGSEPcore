package ngsep.assembly;

import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.Set;

public class LayoutBuilderGreedy implements LayourBuilder 
{	
	@Override
	public void findPaths(AssemblyGraph graph) {
		//Hashmap of edges per vertex
		Map<AssemblyVertex, List<AssemblyEdge>> edges = new HashMap<AssemblyVertex, List<AssemblyEdge>>(); 
		//Adds the current edge to both vertices if they don't have it already in their list of edges
		for (AssemblyEdge assemblyEdge : graph.getEdges()) 
		{
			AssemblyVertex v1 = assemblyEdge.getVertex1();
			AssemblyVertex v2 = assemblyEdge.getVertex2();
			edges.computeIfAbsent(v1, (AssemblyVertex x) -> new ArrayList<AssemblyEdge>()).add(assemblyEdge);
			edges.computeIfAbsent(v2, (AssemblyVertex x) -> new ArrayList<AssemblyEdge>()).add(assemblyEdge);
		}
		
		//List of globally used vertices 
		Set<Integer> usedVerticesGlobal = new HashSet<Integer>(); 
		//Creates contigs until all the vertices are used or don't have available edges
		while(usedVerticesGlobal.size() < graph.getVertices().size() && edges.values().stream().filter(e -> e.size() > 0).count() > 0)
		{			
			//Look for the vertices with the least linked vertices
			List<AssemblyVertex> origins = new ArrayList<AssemblyVertex>();
			for(int i = 1; i < Integer.MAX_VALUE && origins.size() == 0; i++)
			{
				for(Entry<AssemblyVertex, List<AssemblyEdge>> entry : edges.entrySet())
				{
					//Set of linked vertices to the vertex
					Set<Integer> vertices = new HashSet<Integer>();
					for(AssemblyEdge edge : entry.getValue())
					{
						vertices.add(edge.getVertex1().getIndex() == entry.getKey().getIndex() ? edge.getVertex2().getIndex() : edge.getVertex1().getIndex());
					}
					if(vertices.size() == i)
					{
						origins.add(entry.getKey());
					}
				}
			}
			
			//Add path for every origin, picks the largest overlap with an unused node.
			//List of lists of used vertices
			List<Set<Integer>> usedVerticesList = new ArrayList<Set<Integer>>(); 
			//List of lists of possible paths. The one with the most edges visited is chosen.
			List<List<AssemblyEdge>> paths = new ArrayList<List<AssemblyEdge>>(); 
			//Check every origin found and look for the one with most edges visited
			for(AssemblyVertex origin : origins) 
			{
				//Chosen path for origin
				List<AssemblyEdge> ans = new LinkedList<AssemblyEdge>();
				//List of used vertices is initialized with the already accepted used vertices
				Set<Integer> usedVerticesLocal = usedVerticesGlobal;
				//First used vertex is the origin
				usedVerticesLocal.add(origin.getIndex());
				boolean end = false;
				AssemblyVertex pre = origin;
				while(!end)
				{
					AssemblyEdge max = null;
					//Picks the edge with the longest overlap for the current vertex
					for(AssemblyEdge e : edges.get(pre))
					{
						AssemblyVertex post = e.getVertex1().getIndex() == pre.getIndex() ? e.getVertex2() : e.getVertex1();
						if((max == null || e.getOverlap() > max.getOverlap()) && !usedVerticesLocal.contains(post.getIndex()))
							max = e;
					}
					//If there's no edge found, the path ends
					if(max == null)
						end = true;
					else
					{
						ans.add(max);
						usedVerticesLocal.add(pre.getIndex());
						pre = max.getVertex1().getIndex() == pre.getIndex() ? max.getVertex2() : max.getVertex1();
					}
				}
				paths.add(ans);
				usedVerticesList.add(usedVerticesLocal);
			}
			
			//Picks the longest path among the possible origins
			List<AssemblyEdge> maxPath = null;
			Set<Integer> maxUsedVertices = null;
			for(int i = 0; i < paths.size(); i++)
			{
				List<AssemblyEdge> path = paths.get(i);
				if(maxPath == null || maxPath.size() < path.size())
				{
					maxPath = path;
					maxUsedVertices = usedVerticesList.get(i);
				}
			}
			//Adds the longest path to the list of paths of the graph
			if(maxPath != null)
			{
				graph.addPath(maxPath);
				usedVerticesGlobal.addAll(maxUsedVertices);
			}

			//Remove from map all vertices that exist in the used vertices list
			edges.entrySet().removeIf(e -> usedVerticesGlobal.contains(e.getKey().getIndex()));
			//Remove from map all edges that have an used vertex
			for(Entry<AssemblyVertex, List<AssemblyEdge>> entry : edges.entrySet())
			{
				entry.getValue().removeIf(e -> 
					usedVerticesGlobal.contains(e.getVertex1().getIndex()) ||
					usedVerticesGlobal.contains(e.getVertex2().getIndex())
				);
			}
		}
	}
}
