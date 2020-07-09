package ngsep.genome;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.graphs.MCLJob;
import ngsep.math.Distribution;

public class HomologClustersCalculator {
	//Possible program arguments
	private static final int PREFERRED_ORTHOGROUP_SIZE = 50;
	private static final int MAX_SIZE_MCL = 5000;
	
	//Statistics
	private int countLarge = 0;
	private int countSmall = 0;
	private int countMedium = 0;
	private int reProcessedClusters = 0;
	
	//Run parameters
	private boolean skipMCL;
	private Distribution distClusterSizes = new Distribution(0, PREFERRED_ORTHOGROUP_SIZE, 1);
	private int sampleSize = 0;
	
	private Logger log;
	
	/**
	 * @param skipMCL skips the MCL pipeline during clusterHomologs if true
	 */
	public HomologClustersCalculator(boolean skipMCL) {
		this.skipMCL = skipMCL;
	}
	
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}

	public List<List<HomologyUnit>> clusterHomologs(List<AnnotatedReferenceGenome> genomes, List<HomologyEdge> homologyEdges) {
		List<HomologyCatalog> catalogs = new ArrayList<>();
		for(AnnotatedReferenceGenome genome : genomes) catalogs.add(genome.getHomologyCatalog()); 
		return clusterHomologsCatalogs(catalogs, homologyEdges);
	}
	
	public List<List<HomologyUnit>> clusterHomologsCatalogs(List<HomologyCatalog> catalogs, List<HomologyEdge> homologyEdges) {
		log.info("Clustering orthologs and paralogs");
		
		//Set sample size for MCL skip for short clusters
		sampleSize = catalogs.size();
		//Divide homologs into smaller partitions
		List<HomologyUnit> units = new ArrayList<>();
		for(HomologyCatalog catalog : catalogs) units.addAll(catalog.getHomologyUnits()); 
		log.info("Dividing homologs via connected components");
		List<List<HomologyUnit>> partitions = divideUnits(units);
		log.info("Finished dividing homolog units");
		printPartitionsResults(partitions);
		
		//Skip mcl, return connected components
		if(skipMCL) {
			log.info("Skipped MCL, returning connected components instead.");
			return partitions;
		}
		
		log.info("Starting processing partitions");
		ArrayList<PartitionTask> tasks = new ArrayList<PartitionTask>();
		for(List<HomologyUnit> partition : partitions) tasks.add(new PartitionTask(partition, false));
		
		//Infer clusters from each resulting partition task
		List<List<HomologyUnit>> clusters = new ArrayList<List<HomologyUnit>>();
		while(!tasks.isEmpty()) {
			PartitionTask task = tasks.remove(0);
			task = this.processPartition(task);
			if(task.getNewTasks() != null) tasks.addAll(task.getNewTasks());
			if(task.getResults().size() > 0) clusters.addAll(task.getResults());
		}
		
		//Cluster statistics
		for(List<HomologyUnit> cluster : clusters) distClusterSizes.processDatapoint(cluster.size());
		generateStatistics();
		
		return clusters;
	}
	
	private void printPartitionsResults(List<List<HomologyUnit>> partitions) {
		Distribution partitionSizes = new Distribution(0, PREFERRED_ORTHOGROUP_SIZE, 1);
		for(List<HomologyUnit> p : partitions) partitionSizes.processDatapoint(p.size());
		
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		partitionSizes.printDistributionInt(out);
		
		log.info("Partition Size Statistics");
		log.info(String.format("%s", os.toString()));
	}

	/**
	 * Prints out useful statistics from the generated clusters to the logger.
	 */
	private void generateStatistics() {
		log.info("OrthoGroup Results");
		log.info("Size Statistics");
		log.info(String.format("SMALL (2-%d): %d || MEDIUM (%d-MAX): %d || LARGE (%d+): %d", sampleSize, countSmall, sampleSize, countMedium, MAX_SIZE_MCL, countLarge));
		
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		distClusterSizes.printDistributionInt(out);
		
		log.info("Clusters Sizes Statistics");
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
		int markedCount = 0;
		int totalCount = units.size();
		for (int i = 0; i < filteredUnits.size(); i++) {
			if (marked.contains(filteredUnits.get(i).getUniqueKey())) continue;
			Queue<HomologyUnit> queue = new LinkedList<>();
			List<HomologyUnit> currentPartition = new ArrayList<>();
			
			queue.add(filteredUnits.get(i));
			while(!queue.isEmpty()) {
				if(markedCount%1000 == 0) {
					log.info(String.format("Dividing units, current progress: %d/%d. Queue size: %d", markedCount, totalCount, queue.size()));
				}
				
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
						markedCount++;
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
	public PartitionTask processPartition(PartitionTask task) {
		List<HomologyUnit> partition = task.getPartition();
		if (partition.size() <= sampleSize) {
			//Clique
			countSmall++;
			List<List<HomologyUnit>> clusters = new ArrayList<>();
			clusters.add(partition);
			task.setResults(clusters);
		} else if (partition.size() <= MAX_SIZE_MCL){
			//MCL
			countMedium++;
			task = dispatchMCL(task);
		} else {
			log.info(String.format("Passed cluster too large for MCL. Size: %d", partition.size()));
			//Too large for MCL
			countLarge++;
			List<List<HomologyUnit>> clusters = new ArrayList<>();
			clusters.add(partition);
			task.setResults(clusters);
		}
		
		return task;
	}
	
	public PartitionTask dispatchMCL(PartitionTask task) {
		List<HomologyUnit> partition = task.getPartition();
		List<List<HomologyUnit>> clusters = new ArrayList<>();
		
		log.info(String.format("%s partition of size %d. MCL RUN #%d", task.isRequeed()? "Re-Processing" : "Processing", partition.size(), countMedium));
		//Reference index for matrix creation
		HashMap<String, Integer> indexOf = new HashMap<>();
		for(int i = 0; i < partition.size(); i++) { 
			indexOf.put(partition.get(i).getUniqueKey(), i);
		}
		
		//Generating score matrix for MCL
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
		
		//Clustering statistics
		log.info(String.format("Finished MCL RUN #%d. Created %d cluster(s).", countMedium, clusters.size()));
		ArrayList<Integer> shares = new ArrayList<Integer>();
		for(List<HomologyUnit> cluster : clusters) shares.add(cluster.size());
		log.info(Arrays.toString(shares.toArray()));
		
		//Verifying cluster sizes for re-processing
		if(clusters.size() > 1) {
			List<PartitionTask> newTasks = new ArrayList<>();
			Iterator<List<HomologyUnit>> it = clusters.iterator();
			while (it.hasNext()) {
				List<HomologyUnit> cluster = it.next();
				if(cluster.size() > PREFERRED_ORTHOGROUP_SIZE) {
					newTasks.add(new PartitionTask(cleanEdgesCluster(cluster), true));
					reProcessedClusters++;
					it.remove();
					log.info(String.format("Re-Processing scheduled cluster of size %d. Total Re-Processed: %d.", cluster.size(), reProcessedClusters));
				}
			}
			
			task.setResults(clusters);
			task.setNewTasks(newTasks);
		} else {
			task.setResults(clusters);
		}
		
		return task;
	}
	 
	private List<HomologyUnit> cleanEdgesCluster(List<HomologyUnit> cluster) {
		for(HomologyUnit unit : cluster) {
			Collection<HomologyEdge> edges = unit.getAllHomologyRelationships();
			unit.removeAllHomologyRelationships();
			
			for(HomologyEdge edge : edges){
				if (cluster.contains(edge.getSubjectUnit())) {
					unit.addHomologRelationship(new HomologyEdge(edge.getQueryUnit(), edge.getSubjectUnit(), edge.getScore()));
				}
			}
			
			for(HomologyEdge edge : unit.getAllHomologyRelationships())
				if(!cluster.contains(edge.getSubjectUnit())) System.err.println("Well something is clearly not working here, duh");
		}
		return cluster;
	}

	private class PartitionTask {
		private List<HomologyUnit> partition;
		private List<List<HomologyUnit>> results;
		private List<PartitionTask> newTasks;
		private boolean requeed;
		
		public PartitionTask(List<HomologyUnit> partition, boolean requeed) {
			super();
			this.partition = partition;
			this.requeed = requeed;
		}

		public List<HomologyUnit> getPartition() {
			return partition;
		}

		public void setPartition(List<HomologyUnit> partition) {
			this.partition = partition;
		}

		public List<List<HomologyUnit>> getResults() {
			return results;
		}

		public void setResults(List<List<HomologyUnit>> results) {
			this.results = results;
		}

		public List<PartitionTask> getNewTasks() {
			return newTasks;
		}

		public void setNewTasks(List<PartitionTask> newTasks) {
			this.newTasks = newTasks;
		}

		public boolean isRequeed() {
			return requeed;
		}

		public void setRequeed(boolean requeed) {
			this.requeed = requeed;
		}
		
		
	}
}
