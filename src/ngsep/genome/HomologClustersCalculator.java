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
 
	private Logger log;
	
	public Logger getLog() {
		return log;
	}

	public void setLog(Logger log) {
		this.log = log;
	}

	public List<List<HomologyUnit>> clusterHomologs(List<AnnotatedReferenceGenome> genomes, List<HomologyEdge> homologyEdges) {
		log.info("Clustering orthologs and paralogs");
		
		List<HomologyUnit> units = new ArrayList<>();
		for(AnnotatedReferenceGenome genome : genomes) units.addAll(genome.getHomologyUnits()); 
		List<List<HomologyUnit>> partitions = divideUnits(units);
		
		Distribution distClusterSizes = new Distribution(0, 50, 1);
		for(List<HomologyUnit> partition : partitions) distClusterSizes.processDatapoint(partition.size());
		log.info("===== Cluster Distribution =====");
		log.info(String.format("AVG %f || MIN %f || MAX %f || COUNT %f", distClusterSizes.getAverage(), distClusterSizes.getMinValueData(), distClusterSizes.getMaxValueData(), distClusterSizes.getCount()));
		
		
		
		return partitions;
	}
	
	public List<List<HomologyUnit>> oldClusterHomologs(List<AnnotatedReferenceGenome> genomes, List<HomologyEdge> homologyEdges) {
		log.info("Clustering orthologs and paralogs");
		List<List<HomologyUnit>> orthologyUnitClusters = new ArrayList<>();
		
		List<HomologyUnit> reference = new ArrayList<>();
		HashMap<String, Integer> indexOf = new HashMap<>();
		for (AnnotatedReferenceGenome genome : genomes) {
			for (HomologyUnit unit: genome.getHomologyUnits()) {
				String newId = String.format("%d___%s", genome.getId(), unit.getId());
				if(indexOf.containsKey(newId)) {
					log.info(String.format("The id combination is already in use. This is bad, very very bad. == %s", newId));
				}
				indexOf.put(newId, reference.size());
				reference.add(unit);
			}
		}
		
		log.info(String.format("Reference count  == %d", reference.size()));
		log.info(String.format("Edges count == %d ", homologyEdges.size()));
		log.info(String.format("Avg vertex degree == %f", (double)homologyEdges.size()/(double)reference.size()));
		
		int maxNodeDegree = 0;
		SparseMatrix scoreMatrix = new SparseMatrix(reference.size(), reference.size());
		for (int i = 0; i < reference.size(); i++) {
			HomologyUnit unit = reference.get(i);
			Collection<HomologyEdge> edges = unit.getAllHomologyRelationships();
			if(edges.size() > maxNodeDegree) {
				maxNodeDegree = edges.size();
			}
			for(HomologyEdge edge : edges) {
				String id_composed = String.format("%d___%s", edge.getSubjectUnit().getGenomeId(), edge.getSubjectUnit().getId());
				scoreMatrix.set(i, indexOf.get(id_composed), edge.getScore());
			}
		}
		
		//AL AGREGAR EJES, PONER 20 PTS
		
		log.info(String.format("Max vertex degree == %d", maxNodeDegree));
		
		log.info("===== SCORES UNNORMALIZED =====");
		for (String l : scoreMatrix.getMatrixAsString()) {
			log.info(l);
		}
		
		SparseMatrix similarityMatrix = normalizeMatrix(scoreMatrix);
		
		log.info("Similarity matrix is ready");
		log.info("===== Similarity =====");
		for (String l : similarityMatrix.getMatrixAsString()) {
			log.info(l);
		}
		log.info("===== Similarity =====");
		
		SparseMatrix results = simulateFlow(similarityMatrix, 1000, 1);
		log.info("Results are ready");
		log.info("===== Results =====");
		for (String l : results.getMatrixAsString()) {
			log.info(l);
		}
		log.info("===== Results =====");
		
		return orthologyUnitClusters;
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
			List<HomologyUnit> currentCluster = new ArrayList<>();
			
			queue.add(filteredUnits.get(i));
			while(!queue.isEmpty()) {
				HomologyUnit currentUnit = queue.poll();
				if(currentCluster.contains(currentUnit)) {
					//Element already inside cluster.
					continue;
				} else {
					if(marked.get(currentUnit.getUniqueKey()) != null) {
						//Element inside different cluster, merge with current cluster.
						boolean merged = false; 
						for(int j = 0; j < partitions.size() && !merged; j++) {
							List<HomologyUnit> set = partitions.get(j);
							if(set.contains(currentUnit)) {
								merged = true;
								currentCluster.addAll(set);
								partitions.remove(j);
							}
						}
						
						if(!merged) log.warning(String.format("Did not find cluster to merge, but unit was marked. ID: %s", currentUnit.getUniqueKey()));
					} else {
						//Add element to cluster and add its edges to the queue
						marked.put(currentUnit.getUniqueKey(), true);
						currentCluster.add(currentUnit);
						Collection<HomologyEdge> edges = currentUnit.getAllHomologyRelationships();
						for(HomologyEdge edge : edges) {
							queue.add(edge.getSubjectUnit());
						}
					}
				}
			}
			
			if (currentCluster.size() > 0) partitions.add(currentCluster);
		}
			
		return partitions;
	}
	
	public List<List<HomologyUnit>> processPartition(List<HomologyUnit> partition) {
		List<List<HomologyUnit>> clusters = new ArrayList<>();
		
		if (partition.size() <= 10) {
			//Clique
		} else if (partition.size() <= 1000){
			//MCL
			
			HashMap<String, Integer> indexOf = new HashMap<>();
			for(int i = 0; i < partition.size(); i++) { 
				indexOf.put(partition.get(i).getUniqueKey(), i);
			}
			
			float[][] matrix = new float[partition.size()][partition.size()];
			for(int i = 0; i < partition.size(); i++) {
				HomologyUnit currentUnit = partition.get(i);
				for(HomologyEdge edge : currentUnit.getAllHomologyRelationships()) {
					matrix[i][indexOf.get(edge.getQueryUnit().getUniqueKey())] = (float) edge.getScore();
				}
			}
			
			MCLJob job = new MCLJob(matrix);
		}
		
		return clusters;
	}
	
	public SparseMatrix normalizeMatrix(SparseMatrix matrix) {
		SparseMatrix normalized = new SparseMatrix(matrix.length(), matrix.height());
		for(int i = 0; i < matrix.length(); i++) {
			List<ValuePair> tuples = matrix.getRowAsTuples(i);
			double sum = matrix.sumOfRow(i);
			if (sum > 0) {
				double selfLoop = sum/tuples.size();
				sum += (selfLoop);
				normalized.set(i, i, selfLoop/sum);
				for(ValuePair p : tuples) {
					normalized.set(i, p.index, p.value/sum);
				}
			} else {
				normalized.set(i, i, 1);
			}
		}
		return normalized;
	}
	
	public SparseMatrix simulateFlow(SparseMatrix similarityMatrix, int depth, int iterations) {
		SparseMatrix counts = new SparseMatrix(similarityMatrix.length(), similarityMatrix.height());
		
		for(int z = 0; z < iterations; z++) {
			int outOfOddsCount = 0;
			int currentNode = (int)(Math.random()*similarityMatrix.length());
			boolean isDeadNode = true;
			while(isDeadNode) {
				List<ValuePair> spread = similarityMatrix.getRowAsTuples(currentNode);
				isDeadNode = spread.size() <= 1;
				if (isDeadNode) {
					currentNode = (int)(Math.random()*similarityMatrix.length());
				}
			}
			
			
			for(int i = 0; i < depth; i++) {
				List<ValuePair> spread = similarityMatrix.getRowAsTuples(currentNode);
				
				double odds = Math.random(); 
				ValuePair next = null;
				for (int k = 0; k < spread.size() && next == null; k++) {
					ValuePair current = spread.get(k);
					odds -= current.value;
					if (odds <= 0) {
						next = current;
					}
				}
				
				if (next == null) {
					outOfOddsCount++;
					log.info(String.format("FAILED ODDS at %d. Spread was: [%s], odds value was %f", outOfOddsCount, Arrays.toString(spread.toArray()), odds));
					next = spread.get(spread.size()-1);
				}
				counts.set(currentNode, next.index, counts.get(currentNode, next.index) + 1);
				currentNode = next.index;
			}
		}
		
		return counts;
	}
}
