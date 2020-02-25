package ngsep.genome;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.logging.Logger;

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
		List<HomologyUnit> unitsWithOrthologs = new ArrayList<>();
		Map<String, Integer> unitsPositionMap = new HashMap<>();
		for(AnnotatedReferenceGenome genome:genomes) {
			for(HomologyUnit unit:genome.getHomologyUnits()) {
				if(unit.getTotalHomologs()>0) {
					unitsPositionMap.put(unit.getUniqueKey(), unitsWithOrthologs.size());
					unitsWithOrthologs.add(unit);
				}

			}
		}
		log.info("Total units with orthologs: "+unitsWithOrthologs.size());
		int [] group = new int [unitsWithOrthologs.size()];
		Arrays.fill(group, -1);

		Distribution distClusterSizes = new Distribution(0, 50, 1);
		//Build connected components
		for(int i=0;i<group.length;i++) {
			if(group[i]!=-1) continue;
			List<HomologyUnit> cluster = new ArrayList<>();
			int groupNumber = orthologyUnitClusters.size();
			Queue<HomologyUnit> agenda = new LinkedList<>();
			agenda.add(unitsWithOrthologs.get(i));
			while(agenda.size()>0) {
				HomologyUnit unit = agenda.remove();
				int pos = unitsPositionMap.get(unit.getUniqueKey());
				int unitGroup = group[pos];
				if(unitGroup==-1) {
					cluster.add(unit);
					group[pos]=groupNumber;
					//agenda.addAll(unit.getOrthologsAllGenomes());
				} else if(unitGroup!=groupNumber) log.warning("Possible connection between clusters "+unitGroup + " and "+groupNumber+" Unit: "+unit.getUniqueKey());

			}
			if(cluster.size()>1) {
				orthologyUnitClusters.add(cluster);
				distClusterSizes.processDatapoint(cluster.size());
			} else if (cluster.size()==0) {
				log.warning("Empty cluster from unit: "+unitsWithOrthologs.get(i).getUniqueKey()+" clusters: "+orthologyUnitClusters.size()+" groupNumber: "+groupNumber);
			}

		}
		log.info("Number of clusters: "+orthologyUnitClusters.size());
		//TODO: Report it better
		distClusterSizes.printDistributionInt(System.out);
		return orthologyUnitClusters;
	}

}
