package ngsep.genome;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.logging.Logger;

import ngsep.graphs.MCLJob;
import ngsep.graphs.SparseMatrix;
import ngsep.graphs.ValuePair;
import ngsep.math.Distribution;

public class HomologClustersCalculator {
	
	private int countSkipped = 0;
	private int countClique = 0;
	private int countMCL = 0;
 
	private Logger log.info;
	
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}

	public List<List<HomologyUnit>> clusterHomologs(List<AnnotatedReferenceGenome> genomes, List<HomologyEdge> homologyEdges) {
		log.info("Clustering orthologs and paralogs");
		
		//Separate homologs into smaller tasks
		List<HomologyUnit> units = new ArrayList<>();
		for(AnnotatedReferenceGenome genome : genomes) units.addAll(genome.getHomologyUnits()); 
		List<List<HomologyUnit>> partitions = divideUnits(units);
		
		Distribution distPartitionSizes = new Distribution(0, 50, 1);
		for(List<HomologyUnit> partition : partitions) distPartitionSizes.processDatapoint(partition.size());
		log.info("===== Partitions Distribution =====");
		log.info(String.format("AVG %f || MIN %f || MAX %f || COUNT %f", distPartitionSizes.getAverage(), distPartitionSizes.getMinValueData(), distPartitionSizes.getMaxValueData(), distPartitionSizes.getCount()));
		
		List<List<HomologyUnit>> clusters = new ArrayList<List<HomologyUnit>>();
		//Infer clusters from each resulting task
		for(List<HomologyUnit> partition : partitions) {
			List<List<HomologyUnit>> result = this.processPartition(partition);
			if(result.size() > 0) clusters.addAll(result);
		}
		
		Distribution distClusterSizes = new Distribution(0, 50, 1);
		for(List<HomologyUnit> cluster : clusters) distClusterSizes.processDatapoint(cluster.size());
		
		log.info("==== SIZE STATISTICS =====");
		log.info(String.format("SKIPPED: %d || CLIQUES: %d || MCL: %d", countSkipped, countClique, countMCL, clusters.size()));
		log.info("==== CLUSTERING STATISTICS =====");
		log.info(String.format("AVG %f || MIN %f || MAX %f || COUNT %f", distClusterSizes.getAverage(), distClusterSizes.getMinValueData(), distClusterSizes.getMaxValueData(), distClusterSizes.getCount()));
		
		return partitions;
	}
	
	public List<List<HomologyUnit>> divideUnits(List<HomologyUnit> units) {
		List<List<HomologyUnit>> partitions = new ArrayList<>();
		
		List<HomologyUnit> filteredUnits = new ArrayList<>();
		for(HomologyUnit unit : units) {
			if (unit.getTotalHomologs() > 0) {
				filteredUnits.add(unit);
			}
		}
		
		HashMap<String, Boolean> marked = new HashMap<>();
		for (int i = 0; i < filteredUnits.size(); i++) {
			if (marked.get(filteredUnits.get(i).getUniqueKey()) != null) continue;
			Queue<HomologyUnit> queue = new LinkedList<>();
			List<HomologyUnit> currentPartition = new ArrayList<>();
			
			queue.add(filteredUnits.get(i));
			while(!queue.isEmpty()) {
				HomologyUnit currentUnit = queue.poll();
				if(currentPartition.contains(currentUnit)) {
					//Element already inside current partition.
					continue;
				} else {
					if(marked.get(currentUnit.getUniqueKey()) != null) {
						//Element inside different partition, merge with current partition.
						boolean merged = false; 
						for(int j = 0; j < partitions.size() && !merged; j++) {
							List<HomologyUnit> set = partitions.get(j);
							if(set.contains(currentUnit)) {
								merged = true;
								currentPartition.addAll(set);
								partitions.remove(j);
							}
						}
						
						if(!merged) log.warning(String.format("Did not find cluster to merge, but unit was marked. ID: %s", currentUnit.getUniqueKey()));
					} else {
						//Add element to partition and add its edges to the queue
						marked.put(currentUnit.getUniqueKey(), true);
						currentPartition.add(currentUnit);
						Collection<HomologyEdge> edges = currentUnit.getAllHomologyRelationships();
						for(HomologyEdge edge : edges) {
							queue.add(edge.getSubjectUnit());
						}
					}
				}
			}
			
			if (currentPartition.size() > 0) partitions.add(currentPartition);
		}
			
		return partitions;
	}
	
	public List<List<HomologyUnit>> processPartition(List<HomologyUnit> partition) {
		List<List<HomologyUnit>> clusters = new ArrayList<>();
		
		if (partition.size() <= 10) {
			//Clique
			countClique++;
		} else if (partition.size() <= 1000){
			//MCL
			countMCL++;
			
			HashMap<String, Integer> indexOf = new HashMap<>();
			for(int i = 0; i < partition.size(); i++) { 
				indexOf.put(partition.get(i).getUniqueKey(), i);
			}
			
			double[][] matrix = new double[partition.size()][partition.size()];
			for(int i = 0; i < partition.size(); i++) {
				HomologyUnit currentUnit = partition.get(i);
				for(HomologyEdge edge : currentUnit.getAllHomologyRelationships()) {
					matrix[i][indexOf.get(edge.getQueryUnit().getUniqueKey())] = edge.getScore();
				}
			}
			
			MCLJob job = new MCLJob(matrix);
			job.run();
			
			List<List<Integer>> results = job.getResults();
			for(List<Integer> indexList : results) {
				List<HomologyUnit> cluster = new ArrayList<>();
				for(Integer k : indexList) cluster.add(partition.get(k));
				clusters.add(cluster);
			}
		} else {
			countSkipped++;
		}
		
		return clusters;
	}
}
