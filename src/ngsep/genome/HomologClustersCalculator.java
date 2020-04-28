package ngsep.genome;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
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
	//Statistics
	private int countLarge = 0;
	private int countSmall = 0;
	private int countMedium = 0;
	
	//Run parameters
	private boolean skipMCL;
	private Distribution distMCLSpread = new Distribution(0, 1, 0.01);
	private Distribution distMCLCount = new Distribution(1, 100, 1);
	private Distribution distClusterSizes = new Distribution(0, 1000, 1);
	
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
		
		//Divide homologs into smaller partitions
		List<HomologyUnit> units = new ArrayList<>();
		for(HomologyCatalog catalog : catalogs) units.addAll(catalog.getHomologyUnits()); 
		log.info("Dividing homologs via connected components");
		List<List<HomologyUnit>> partitions = divideUnits(units);
		log.info("Finished dividing homolog units");
		
		//Skip mcl, return connected components
		if(skipMCL) {
			log.info("Skipped MCL, returning connected components instead.");
			return partitions;
		}
		
		log.info("Starting processing partitions");
		List<List<HomologyUnit>> clusters = new ArrayList<List<HomologyUnit>>();
		//Infer clusters from each resulting partition
		for(List<HomologyUnit> partition : partitions) {
			log.info(String.format("Processing partition of size %d", partition.size()));
			List<List<HomologyUnit>> result = this.processPartition(partition);
			if(result.size() > 0) clusters.addAll(result);
		}
		
		//Cluster statistics
		for(List<HomologyUnit> cluster : clusters) distClusterSizes.processDatapoint(cluster.size());
		generateStatistics(clusters);
		
		return clusters;
	}
	
	/**
	 * Prints out useful statistics from the generated clusters to the logger.
	 * @param clusters resulting clusters.
	 */
	private void generateStatistics(List<List<HomologyUnit>> clusters) {
		log.info("OrthoGroup Results");
		log.info("Size Statistics");
		log.info(String.format("SMALL (2-10): %d || MEDIUM (11-5000): %d || LARGE (5000+): %d", countSmall, countMedium, countLarge));
		
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		PrintStream out = new PrintStream(os);
		distMCLCount.printDistributionInt(out);
		
		log.info("MCL Sizes Statistics");
		log.info(String.format("%s", os.toString()));
		
		os = new ByteArrayOutputStream();
		out = new PrintStream(os);
		distMCLSpread.printDistributionInt(out);
		
		log.info("MCL Spread Statistics");
		log.info(String.format("%s", os.toString()));
		
		os = new ByteArrayOutputStream();
		out = new PrintStream(os);
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
	public List<List<HomologyUnit>> processPartition(List<HomologyUnit> partition) {
		List<List<HomologyUnit>> clusters = new ArrayList<>();
		
		if (partition.size() <= 10) {
			//Clique (needs to be verified)
			countSmall++;
			clusters.add(partition);
		} else if (partition.size() <= 5000){
			//MCL
			countMedium++;
			
			log.info(String.format("Processing partition of size %d. MCL RUN #%d", partition.size(), countMedium));
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
			
			//Clustering statistics
			log.info(String.format("Finished MCL RUN #%d. Created %d clusters.", countMedium, clusters.size()));
			distMCLCount.processDatapoint(clusters.size());
			ArrayList<Integer> shares = new ArrayList<Integer>();
			for(List<HomologyUnit> cluster : clusters) {
				double val = ((double)cluster.size())/((double)clusters.size());
				shares.add(cluster.size());
				distMCLSpread.processDatapoint(val);
			}
			log.info(Arrays.toString(shares.toArray()));
		} else {
			//Too large for MCL
			countLarge++;
			clusters.add(partition);
		}
		return clusters;
	}
}
