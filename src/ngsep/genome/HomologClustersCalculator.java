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
				String newId = String.format("%s___%s", genome.getId(), unit.getId());
				if(indexOf.containsKey(newId)) log.info("The id combination is already in use. This is bad, very very bad.");
				indexOf.put(newId, reference.size());
				reference.add(unit);
			}
		}
		
		SparseMatrix similarityMatrix = new SparseMatrix(reference.size(), reference.size());
		for (int i = 0; i < reference.size(); i++) {
			HomologyUnit unit = reference.get(i);
			Collection<HomologyEdge> edges = unit.getAllHomologyRelationships();
			for(HomologyEdge edge : edges) {
				String id_composed = String.format("%s___%s", edge.getSubjectUnit().getGenomeId(), edge.getSubjectUnit().getId());
				similarityMatrix.set(i, indexOf.get(id_composed), edge.getScore());
			}
		}
		
		log.info("==== Similarity matrix ====");
		log.info(similarityMatrix.getMatrixAsString());
		
		return orthologyUnitClusters;
	}
}
