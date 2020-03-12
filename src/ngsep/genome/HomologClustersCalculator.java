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

import ngsep.genome.strucs.SparseMatrix;
import ngsep.genome.strucs.ValuePair;
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
		List<List<HomologyUnit>> orthologyUnitClusters=new ArrayList<>();
		
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
		
		SparseMatrix scoreMatrix = new SparseMatrix(reference.size(), reference.size());
		for (int i = 0; i < reference.size(); i++) {
			HomologyUnit unit = reference.get(i);
			Collection<HomologyEdge> edges = unit.getAllHomologyRelationships();
			for(HomologyEdge edge : edges) {
				String id_composed = String.format("%d___%s", edge.getSubjectUnit().getGenomeId(), edge.getSubjectUnit().getId());
				scoreMatrix.set(i, indexOf.get(id_composed), edge.getScore());
			}
		}
		SparseMatrix similarityMatrix = normalizeMatrix(scoreMatrix);
		
		log.info("Similarity matrix is ready");
		SparseMatrix results = simulateFlow(similarityMatrix, 1000, 10);
		log.info("Results are ready");
		log.info("===== Results =====");
		for (String l : results.getMatrixAsString()) {
			log.info(l);
		}
		log.info("===== Results =====");
		
		return orthologyUnitClusters;
	}
	
	public SparseMatrix normalizeMatrix(SparseMatrix matrix) {
		SparseMatrix normalized = new SparseMatrix(matrix.length(), matrix.height());
		for(int i = 0; i < matrix.length(); i++) {
			List<ValuePair> tuples = matrix.getRowAsTuples(i);
			double sum = matrix.sumOfRow(i);
			double selfLoop = sum/tuples.size();
			sum += (selfLoop);
			normalized.set(i, i, selfLoop/sum);
			for(ValuePair p : tuples) {
				normalized.set(i, p.index, p.value/sum);
			}
		}
		return normalized;
	}
	
	public SparseMatrix simulateFlow(SparseMatrix similarityMatrix, int depth, int iterations) {
		SparseMatrix counts = new SparseMatrix(similarityMatrix.length(), similarityMatrix.height());
		
		for(int z = 0; z < iterations; z++) {
			int outOfOddsCount = 0;
			int currentNode = (int) Math.random()*similarityMatrix.length();
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
					log.info(String.format("Failed to set odds at %d. Spread was: [%s], odds value was %f", outOfOddsCount, Arrays.toString(spread.toArray()), odds));
					next = spread.get(spread.size()-1);
				}
				counts.set(currentNode, next.index, counts.get(currentNode, next.index) + 1);
				currentNode = next.index;
			}
		}
		
		return counts;
	}
}
