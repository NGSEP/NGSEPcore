/*******************************************************************************
 * NGSEP - Next Generation Sequencing Experience Platform
 * Copyright 2016 Jorge Duitama
 *
 * This file is part of NGSEP.
 *
 *     NGSEP is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     NGSEP is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with NGSEP.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
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
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.logging.Logger;

import ngsep.graphs.MCLJob;
import ngsep.math.Distribution;

/**
 * 
 * @author Jorge Gomez
 *
 */
public class HomologClustersCalculator {
	//Possible program arguments
	private static final int PREFERRED_ORTHOGROUP_SIZE = 50;
	private static final int MAX_SIZE_MCL = 1000;
	
	//Statistics
	private int countLarge = 0;
	private int countSmall = 0;
	private int countMedium = 0;
	private int reProcessedClusters = 0;
	
	//Run parameters
	private boolean skipMCL;
	private Distribution distClusterSizes = new Distribution(0, PREFERRED_ORTHOGROUP_SIZE, 1);
	private int sampleSize = 0;
	private HomologRelationshipsFinder finder;
	private boolean removeEdgesOutsideClusters = false;
	
	private Logger log = Logger.getAnonymousLogger();
	
	/**
	 * @param skipMCL skips the MCL pipeline during clusterHomologs if true
	 */
	public HomologClustersCalculator(HomologRelationshipsFinder finder, boolean skipMCL) {
		this.finder = finder;
		this.skipMCL = skipMCL;
	}
	
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}
	

	public List<HomologyCluster> clusterHomologs(List<AnnotatedReferenceGenome> genomes) {
		List<HomologyCatalog> catalogs = new ArrayList<>();
		for(AnnotatedReferenceGenome genome : genomes) catalogs.add(genome.getHomologyCatalog()); 
		return clusterHomologsCatalogs(catalogs);
	}
	
	public List<HomologyCluster> clusterHomologsCatalogs(List<HomologyCatalog> catalogs) {
		log.info("Clustering orthologs and paralogs");
		
		//Set sample size for MCL skip for short clusters
		sampleSize = catalogs.size();
		//Divide homologs into smaller partitions
		List<HomologyUnit> units = new ArrayList<>();
		for(HomologyCatalog catalog : catalogs) units.addAll(catalog.getHomologyUnits());
		makeBidirectional(units);
		log.info("Dividing homologs via connected components");
		List<List<HomologyUnit>> partitions = divideUnits(units);
		//for(List<HomologyUnit> partition : partitions) updateEdges(partition);
		//Second round with better scores
		//makeBidirectional(units);
		//partitions = divideUnits(units);
		//for(List<HomologyUnit> partition : partitions) updateEdges(partition);
		log.info("Finished dividing homolog units. Number of clusters: "+partitions.size());
		printPartitionsResults(partitions);
		
		//partitions = breakPartitionsByDegree (partitions);
		//Skip mcl, return connected components
		if(skipMCL) {
			log.info("Skipped MCL, returning connected components instead.");
			
			return convertListHUtoHC(partitions);
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
		if (removeEdgesOutsideClusters) {
			for(List<HomologyUnit> partition : clusters) updateEdges(partition);
		}
		//Cluster statistics
		for(List<HomologyUnit> cluster : clusters) distClusterSizes.processDatapoint(cluster.size());
		generateStatistics();
		
		return convertListHUtoHC(clusters);
	}
	
	private void makeBidirectional(List<HomologyUnit> units) {
		for(HomologyUnit unit1:units) {
			Collection<HomologyEdge> edges = unit1.getAllHomologyRelationships();
			for(HomologyEdge edge:edges) {
				HomologyUnit u2 = edge.getSubjectUnit();
				if(u2.getHomologyEdge(unit1)==null) {
					u2.addHomologRelationship(new HomologyEdge(u2, unit1, 0.01));
				}
			}
		}
		
	}

	private List<List<HomologyUnit>> breakPartitionsByDegree(List<List<HomologyUnit>> partitions) {
		List<List<HomologyUnit>> answer = new ArrayList<List<HomologyUnit>>();
		for(List<HomologyUnit> cluster : partitions) {
			int n = cluster.size();
			if(n<20 || n<=2*sampleSize) {
				answer.add(cluster);
				continue;
			}
			//TODO: Divide by node degree peaks
			List<HomologyUnit> highDegreeNodes = new ArrayList<HomologyUnit>();
			List<HomologyUnit> lowDegreeNodes = new ArrayList<HomologyUnit>();
			for(HomologyUnit unit:cluster) {
				double degree = unit.getTotalHomologs();
				if(degree < n/10) {
					lowDegreeNodes.add(unit);
				} else {
					highDegreeNodes.add(unit);
				}
			}
			if(highDegreeNodes.size()>0) {
				answer.add(highDegreeNodes);
				if(lowDegreeNodes.size()>0) {
					List<List<HomologyUnit>> lowDegreepartitions = divideUnits(lowDegreeNodes);
					answer.addAll(lowDegreepartitions);
					log.info("Breaking by degree cluster with "+cluster.size()+" units. High degree: "+highDegreeNodes.size()+" low degree: "+lowDegreeNodes.size()+" partitions: "+lowDegreepartitions.size());
				} else {
					log.info("Dense cluster with "+cluster.size()+" units.");
				}
				
				
			} else {
				log.info("Sparse cluster with "+cluster.size()+" units.");
				answer.add(cluster);
			}
			
			
		}
		
		return answer;
	}

	private void updateEdges(List<HomologyUnit> partition) { 
		Map<Long,Set<Integer>> unitsByKmerCode = finder.indexKmersHomologyUnits(partition);
		for(HomologyUnit unit:partition) unit.removeAllHomologyRelationships();
		finder.calculateHomologs(partition, partition, unitsByKmerCode);
	}

	/**
	 * Convert a list of lists of homology units to a list of homology clusters
	 * @param clusters
	 * @return
	 */
	private List<HomologyCluster> convertListHUtoHC(List<List<HomologyUnit>> clusters){
		 List<HomologyCluster> homologyClusters = new ArrayList<>();
		 int count = 0;
		 for(List<HomologyUnit> cluster : clusters) {
			 HomologyCluster newCluster = new HomologyCluster(count, cluster);
			 count++;
			 homologyClusters.add(newCluster);
		 }
		 return homologyClusters;
		
	}
	
	private void printPartitionsResults(List<List<HomologyUnit> > partitions) {
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
		
		Set<String> inputUnitIds = new HashSet<String>();
		for(HomologyUnit unit : units) inputUnitIds.add(unit.getUniqueKey());
		
		Map<String,Integer> clustersMap = new HashMap<String, Integer>();
		Set<String> marked = new HashSet<String>();
		int totalCount = units.size();
		for (int i = 0; i < units.size(); i++) {
			HomologyUnit nextStart = units.get(i); 
			if (clustersMap.containsKey(nextStart.getUniqueKey())) continue;
			//System.out.println("Next start "+nextStart.getUniqueKey());
			Queue<HomologyUnit> queue = new LinkedList<>();
			List<HomologyUnit> currentPartition = new ArrayList<>();
			
			queue.add(nextStart);
			marked.add(nextStart.getUniqueKey());
			if(marked.size()%1000 == 0) log.info(String.format("Dividing units, current progress: %d/%d. Queue size: %d", marked.size(), totalCount, queue.size()));
			while(!queue.isEmpty()) {
				HomologyUnit currentUnit = queue.poll();
				String key = currentUnit.getUniqueKey();
				if(clustersMap.containsKey(key)) {
					//Element inside different partition, merge with current partition.
					log.warning(String.format("Processing unit ID: %s assigned to cluster %d. Current cluster: %d", key,clustersMap.get(key),partitions.size()));
					continue;
				}
				//Add element to partition and add its edges to the queue
				currentPartition.add(currentUnit);
				clustersMap.put(key,partitions.size());
				//System.out.println("Adding key "+key + " to cluster:" + partitions.size() );
				Collection<HomologyEdge> edges = currentUnit.getAllHomologyRelationships();
				for(HomologyEdge edge : edges) {
					HomologyUnit su = edge.getSubjectUnit();
					String skey = su.getUniqueKey();
					if(inputUnitIds.contains(skey) && !marked.contains(skey)) {
						queue.add(su);
						marked.add(skey);
						if(marked.size()%1000 == 0) log.info(String.format("Dividing units, current progress: %d/%d. Queue size: %d", marked.size(), totalCount, queue.size()));
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
				String key = edge.getSubjectUnit().getUniqueKey();
				Integer j = indexOf.get(key);
				if(j!=null)matrix[i][j] = edge.getScore();
			}
		}
		log.info(String.format("Matrix built for MCL RUN #%d.", countMedium, clusters.size()));
		
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
		
	}
}
