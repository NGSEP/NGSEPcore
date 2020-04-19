package ngsep.genome;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.graphs.MCLJob;
import ngsep.math.Distribution;

public class HomologClustersCalculator {
	
	private int countLarge = 0;
	private int countSmall = 0;
	private int countMedium = 0;
 
	private Logger log;
	
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}

	public List<List<HomologyUnit>> clusterHomologs(List<AnnotatedReferenceGenome> genomes, List<HomologyEdge> homologyEdges) {
		log.info("Clustering orthologs and paralogs");
		
		//Separate homologs into smaller partitions
		List<HomologyUnit> units = new ArrayList<>();
		for(AnnotatedReferenceGenome genome : genomes) units.addAll(genome.getHomologyUnits()); 
		List<List<HomologyUnit>> partitions = divideUnits(units);
		
		List<List<HomologyUnit>> clusters = new ArrayList<List<HomologyUnit>>();
		//Infer clusters from each resulting partition
		for(List<HomologyUnit> partition : partitions) {
			List<List<HomologyUnit>> result = this.processPartition(partition);
			if(result.size() > 0) clusters.addAll(result);
		}
		generateStatistics(clusters);
		
		return clusters;
	}
	
	public List<List<HomologyUnit>> clusterHomologsOrganisms(List<HomologyCatalog> organisms, List<HomologyEdge> homologyEdges) {
		log.info("Clustering orthologs and paralogs");
		
		//Separate homologs into smaller partitions
		List<HomologyUnit> units = new ArrayList<>();
		for(HomologyCatalog organism : organisms) units.addAll(organism.getHomologyUnits()); 
		List<List<HomologyUnit>> partitions = divideUnits(units);
		
		List<List<HomologyUnit>> clusters = new ArrayList<List<HomologyUnit>>();
		//Infer clusters from each resulting partition
		for(List<HomologyUnit> partition : partitions) {
			List<List<HomologyUnit>> result = this.processPartition(partition);
			if(result.size() > 0) clusters.addAll(result);
		}
		generateStatistics(clusters);
		
		return clusters;
	}
	
	/**
	 * Prints out useful statistics from the generated clusters to the logger.
	 * @param clusters resulting clusters.
	 */
	private void generateStatistics(List<List<HomologyUnit>> clusters) {
		Distribution distClusterSizes = new Distribution(0, 1000, 1);
		for(List<HomologyUnit> cluster : clusters) distClusterSizes.processDatapoint(cluster.size());
		
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		distClusterSizes.printDistributionInt(out);
		
		log.info("OrthoGroup Results");
		log.info("Size Statistics");
		log.info(String.format("SMALL (2-10): %d || MEDIUM (11-1000): %d || LARGE (1000+): %d", countSmall, countMedium, countLarge));
		log.info("Clustering Statistics");
		log.info(String.format("%s", os.toString()));
	}

	/**
	 * Takes a set of homology units and their edges and divides them into connected components.
	 * @param units units to be used.
	 * @return Lists of each division made.
	 */
	public List<List<HomologyUnit>> divideUnits(List<HomologyUnit> units) {
		List<List<HomologyUnit>> partitions = new ArrayList<>();
		
		List<HomologyUnit> filteredUnits = new ArrayList<>();
		for(HomologyUnit unit : units) {
			if (unit.getTotalHomologs() > 0) {
				filteredUnits.add(unit);
			}
		}
		
		Set<String> marked = new HashSet<>();
		for (int i = 0; i < filteredUnits.size(); i++) {
			if (marked.contains(filteredUnits.get(i).getUniqueKey())) continue;
			Queue<HomologyUnit> queue = new LinkedList<>();
			List<HomologyUnit> currentPartition = new ArrayList<>();
			
			queue.add(filteredUnits.get(i));
			while(!queue.isEmpty()) {
				HomologyUnit currentUnit = queue.poll();
				if(currentPartition.contains(currentUnit)) {
					//Element already inside current partition.
					continue;
				} else {
					if(marked.contains(currentUnit.getUniqueKey())) {
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
						
						if(!merged) log.warning(String.format("Did not find partition to merge, but unit was marked. ID: %s", currentUnit.getUniqueKey()));
					} else {
						//Add element to partition and add its edges to the queue
						marked.add(currentUnit.getUniqueKey());
						currentPartition.add(currentUnit);
						Collection<HomologyEdge> edges = currentUnit.getAllHomologyRelationships();
						for(HomologyEdge edge : edges) {
							queue.add(edge.getSubjectUnit());
						}
					}
				}
			}
			
			if (currentPartition.size() > 1) partitions.add(currentPartition);
		}
			
		return partitions;
	}
	
	/**
	 * Takes a partition of any size and tries to break it into clusters
	 * @param partition the partition to be used
	 * @return a list of clusters which add up to the starting partition
	 */
	public List<List<HomologyUnit>> processPartition(List<HomologyUnit> partition) {
		List<List<HomologyUnit>> clusters = new ArrayList<>();
		
		if (partition.size() <= 10) {
			//Clique (needs to be verified)
			countSmall++;
			clusters.add(partition);
		} else if (partition.size() <= 1000){
			//MCL
			countMedium++;
			
			log.info(String.format("MCL RUN #%d", countMedium));
			HashMap<String, Integer> indexOf = new HashMap<>();
			for(int i = 0; i < partition.size(); i++) { 
				indexOf.put(partition.get(i).getUniqueKey(), i);
			}
			
			double[][] matrix = new double[partition.size()][partition.size()];
			for(int i = 0; i < partition.size(); i++) {
				HomologyUnit currentUnit = partition.get(i);
				for(HomologyEdge edge : currentUnit.getAllHomologyRelationships()) {
					matrix[i][indexOf.get(edge.getSubjectUnit().getUniqueKey())] = edge.getScore();
				}
			}
			
			MCLJob job = new MCLJob(matrix);
			job.run();
			
			List<List<Integer>> results = job.getResults();
			for(List<Integer> indexList : results) {
				List<HomologyUnit> cluster = new ArrayList<>();
				for(Integer k : indexList) cluster.add(partition.get(k));
				if(cluster.size() > 1) clusters.add(cluster);
			}
		} else {
			//Too large for MCL, to be defined later.
			countLarge++;
			clusters.add(partition);
		}
		
		return clusters;
	}
}
