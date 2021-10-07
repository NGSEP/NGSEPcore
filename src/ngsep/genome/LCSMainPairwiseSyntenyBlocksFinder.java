package ngsep.genome;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.logging.Logger;

import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class LCSMainPairwiseSyntenyBlocksFinder implements PairwiseSyntenyBlocksFinder {

	private Logger log = Logger.getAnonymousLogger();
	@Override
	public List<PairwiseSyntenyBlock> findSyntenyBlocks(AnnotatedReferenceGenome genome1, AnnotatedReferenceGenome genome2, List<HomologyCluster> clusters) {
		List<PairwiseSyntenyBlock> answer = new ArrayList<PairwiseSyntenyBlock>();
		Map<String,HomologyCluster> unitsToClusters = new HashMap<>();
		for(int i=0;i<clusters.size();i++) {
			HomologyCluster cluster = clusters.get(i);
			for(HomologyUnit unit:cluster.getHomologyUnitsCluster()) {
				unitsToClusters.put(unit.getUniqueKey(), cluster);
			}
		}
		Set<String> pairedSequenceNamesG2 = new HashSet<String>();
		QualifiedSequenceList sequencesG1 = genome1.getSequencesMetadata();
		for(QualifiedSequence chrG1:sequencesG1) {
			List<HomologyUnit> unitsChrG1 = genome1.getHomologyUnits(chrG1.getName());
			log.info("Total units G1 for "+chrG1.getName()+": "+unitsChrG1.size());
			//Whole chromosome synteny LCS algorithm
			String chrNameG2 = findBestChromosome(unitsChrG1, genome2.getId());
			if(chrNameG2!=null) {

				PairwiseSyntenyBlock block = alignOrthologyClusters(unitsToClusters,genome1.getId(), unitsChrG1,genome2.getId(),chrNameG2);
				List<HomologyUnit> unitsChrG2 = genome2.getHomologyUnits(chrNameG2);
				log.info("Sequence "+chrG1.getName()+" in first genome aligned to sequence "+chrNameG2+" in the second genome. Orthology units sequence genome 1 "+unitsChrG1.size()+". Orthology units sequence genome 2: "+unitsChrG2.size()+" LCS size: "+block.getHomologies().size());
				answer.add(block);
				pairedSequenceNamesG2.add(chrNameG2);
			} else {
				log.info("Mate sequence not found for "+chrG1.getName()+" Sequence orthology units: "+unitsChrG1.size());
			}
			
		}
		QualifiedSequenceList sequencesG2 = genome2.getSequencesMetadata();
		for(QualifiedSequence chrG2:sequencesG2) {
			List<HomologyUnit> unitsChrG2 = genome2.getHomologyUnits(chrG2.getName());
			if(pairedSequenceNamesG2.contains(chrG2.getName())) continue;
			log.info("Total units G2 for "+chrG2.getName()+": "+unitsChrG2.size());
			//Whole chromosome synteny LCS algorithm
			String chrNameG1 = findBestChromosome(unitsChrG2, genome1.getId());
			if(chrNameG1!=null) {
				PairwiseSyntenyBlock block = alignOrthologyClusters(unitsToClusters,genome2.getId(), unitsChrG2,genome1.getId(),chrNameG1);
				List<HomologyUnit> unitsChrG1 = genome2.getHomologyUnits(chrNameG1);
				log.info("Sequence "+chrG2.getName()+" in second genome aligned to sequence "+chrNameG1+" in the first genome. Orthology units sequence genome 1 "+unitsChrG1.size()+". Orthology units sequence genome 2: "+unitsChrG2.size()+" LCS size: "+block.getHomologies().size());
				answer.add(block);
			} else {
				log.info("Mate sequence not found for "+chrG2.getName()+" Sequence orthology units: "+unitsChrG2.size());
			}
		}
		return answer;
	}
	
	/**
	 * Selects the chromosome having the largest number of mates with the given units
	 * @param units to select the chromosome with the best fit
	 * @param genomeId Id of the genome to query orthologs
	 * @return String name of the most frequent chromosome in the mates of the given units
	 */
	private String findBestChromosome(List<HomologyUnit> units, int genomeId) {
		Map<String,Integer> chrMateCounts = new HashMap<>();

		//Go over the orthology units. Locate chromosome of each unique ortholog and update counts map

		for(HomologyUnit unit:units)
		{
			HomologyUnit mate = unit.getUniqueOrtholog(genomeId);
			if(mate == null) continue;
			String sequenceName = mate.getSequenceName();			
			if((chrMateCounts.containsKey(sequenceName)))
			{
				int value = chrMateCounts.get(sequenceName)+1;
				chrMateCounts.put(sequenceName, value);
			}
			else
			{
				chrMateCounts.put(sequenceName, 1);
			}

		}
		//Find the chromosome with the largest count
		int bestChromosomeCounts = 0;
		String bestChromosome = null;

		for(Map.Entry<String,Integer> entry : chrMateCounts.entrySet()) {
			String name = entry.getKey();
			int chromosomeProteins = entry.getValue();

			if(chromosomeProteins > bestChromosomeCounts)
			{
				bestChromosomeCounts = chromosomeProteins;
				bestChromosome = name;
			}
		}

		return bestChromosome;
	}
	
	/**
	 * Aligns the orthology clusters from two homologous chromosomes using LCS
	 * @param genome1Id Id of the first genome
	 * @param unitsChrG1 This list has only one chromosome and is sorted by position
	 * @param genome2Id Id of the second genome
	 * @param seqName2 Name of the sequence to align in the second genome
	 * @return List<OrthologyUnit> List of selected units of the first list making the LCS relative to the second list
	 */
	private PairwiseSyntenyBlock alignOrthologyClusters(Map<String,HomologyCluster> unitsToClusters, int genomeId1, List<HomologyUnit> unitsChrG1,  int genomeId2, String seqName2) {

		List<LocalHomologyCluster> unitsG1List = new ArrayList<>();
		List<LocalHomologyCluster> unitsG2List = new ArrayList<>();
		Map<Integer,Integer> clusterIdsToG1Indexes = new HashMap<>();
		Map<Integer,LocalHomologyCluster> clusterIdsToG2LocalUnits = new HashMap<>();
		for(int i=0; i<unitsChrG1.size(); i++)
		{
			HomologyUnit unitG1 = unitsChrG1.get(i);
			HomologyCluster cluster = unitsToClusters.get(unitG1.getUniqueKey());
			if(cluster == null) continue;
			if(clusterIdsToG1Indexes.containsKey(cluster.getClusterId())) continue;
			LocalHomologyCluster localCluster1 = cluster.getLocalCluster(unitG1);
			if(localCluster1==null) continue;
			
			List<LocalHomologyCluster> localClusters2 = cluster.getLocalClusters(genomeId2, seqName2);
			if(localClusters2==null || localClusters2.size()==0) continue;
			LocalHomologyCluster localCluster2 = localClusters2.get(0);
			int clusterId = cluster.getClusterId();
			clusterIdsToG1Indexes.put(clusterId,unitsG1List.size());
			clusterIdsToG2LocalUnits.put(clusterId,localCluster2);
			unitsG1List.add(localCluster1);			
			unitsG2List.add(localCluster2);
		}
		//Not needed because they are already sorted
		//Collections.sort(unitsG1List, GenomicRegionPositionComparator.getInstance());
		Collections.sort(unitsG2List, GenomicRegionPositionComparator.getInstance());		


		Set<Integer> lcsForward = findLCS(unitsG1List, unitsG2List, clusterIdsToG1Indexes);

		Collections.reverse(unitsG2List);

		Set<Integer> lcsReverse = findLCS(unitsG1List, unitsG2List, clusterIdsToG1Indexes);

		Set<Integer> lcs = lcsForward;
		boolean reverse = lcsReverse.size()>lcsForward.size();
		if(reverse) lcs = lcsReverse;
		
		//System.out.println("Positions for LCS: "+positions.length+" LCS: "+lcs.size());
		// Select the orthology units in G1 located at the indexes given by the output of LCS
		List<SyntenyVertex> block = new ArrayList<SyntenyVertex>();
		Set<String> unitsInLCS = new HashSet<String>();
		for(int i:lcs) {
			LocalHomologyCluster localCluster1 = unitsG1List.get(i);
			int clusterId = localCluster1.getParent().getClusterId();
			LocalHomologyCluster localCluster2 = clusterIdsToG2LocalUnits.get(clusterId);
			SyntenyVertex vertex = new SyntenyVertex(localCluster1, localCluster2);
			block.add(vertex);
			for(HomologyUnit unit:localCluster1.getHomologyUnitsCluster()) unitsInLCS.add(unit.getUniqueKey());
			for(HomologyUnit unit:localCluster2.getHomologyUnitsCluster()) unitsInLCS.add(unit.getUniqueKey());
		}
		block = completeLCS(block, unitsToClusters, genomeId1, unitsChrG1, genomeId2, seqName2, unitsInLCS, reverse);
		return new PairwiseSyntenyBlock(block);
	}
	
	

	private Set<Integer> findLCS(List<LocalHomologyCluster> unitsG1List, List<LocalHomologyCluster> unitsG2List, Map<Integer,Integer> clusterIdsToG1Indexes) {
		// Create int array with the positions in g2 for the units in g1. Input for LCS
		int []positions = new int[unitsG1List.size()];
		for(int j=0; j<unitsG2List.size(); j++) {
			LocalHomologyCluster unit=unitsG2List.get(j);
			int clusterId = unit.getParent().getClusterId();
			int i = clusterIdsToG1Indexes.get(clusterId);
			positions[i] = j;
		}
		// Run LCS

		Set<Integer> lcs = findLCS(positions);
		return lcs;
	}
	
	/**
	 * Calculates the longest common subsequence (LCS) of sorted entries in the given indexes array using dynamic programming 
	 * @param indexesMap Indexes to find the LCS
	 * @return SortedSet<Integer> Positions making the LCS
	 */
	public SortedSet<Integer> findLCS (int [] indexesMap) {
		SortedSet<Integer> answer = new TreeSet<>();
		int n = indexesMap.length;
		int [] [] m = new int [n][n+1];
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<n+1; j++)
			{
				if(i==0)
				{
					if(j==0 || indexesMap[i]>0) {
						m[i][j] = 0;
					}
					else {
						m[i][j] = 1;
					}

				}
				else if(j<=indexesMap[i])
				{
					m[i][j] = m[i-1][j];
				}
				else
				{
					m[i][j] = Math.max((m[i-1][j]),(m[i-1][indexesMap[i]]+1));
				}
				//System.out.print(" "+m[i][j]);
			}
			//System.out.println();
		}

		int i = m.length-1;
		int j = m[0].length-1;
		//System.out.println("LCS matrix size: "+i+"-"+j);
		while(i > 0 && j > 0)
		{
			//System.out.print("Position: "+i+"-"+j);
			int up = m[i-1][j];
			int diag = -1;
			if(j>indexesMap[i]) {
				diag = m[i-1][indexesMap[i]]+1;
			}
			if(diag >= up )
			{
				answer.add(i);
				j=indexesMap[i];
				i--;
			} else {
				i--;
			}
			//System.out.println(" Diag score: "+diag+" up score: "+up+" size answer: "+answer.size()+" next Position: "+i+"-"+j);
		}
		if(i==0 && indexesMap[i]==0)
		{
			answer.add(i);
		}
		return answer;
	}

	private List<SyntenyVertex> completeLCS(List<SyntenyVertex> vertices, Map<String,HomologyCluster> unitsToClusters, int genomeId1, List<HomologyUnit> unitsChrG1,  int genomeId2, String seqName2, Set<String> unitsInLCS, boolean reverse) {
		List<SyntenyVertex> answer = new ArrayList<SyntenyVertex>();
		int i = 0;
		int j=0;
		LocalHomologyCluster lastInSynteny2 = null;
		while(i<unitsChrG1.size() && j<vertices.size()) {
			//Recover data from next unit
			HomologyUnit unitG1 = unitsChrG1.get(i);
			SyntenyVertex vertex = vertices.get(j);
			LocalHomologyCluster nextInSynteny1 = vertex.getLocalRegion1();
			LocalHomologyCluster nextInSynteny2 = vertex.getLocalRegion2();
			if(nextInSynteny1.getFirst()<unitG1.getFirst()) {
				answer.add(vertex);
				lastInSynteny2 = nextInSynteny2;
				j++;
				continue;
			}
			if(unitsInLCS.contains(unitG1.getUniqueKey())) {
				i++;
				continue;
			}
			HomologyCluster cluster = unitsToClusters.get(unitG1.getUniqueKey());
			if(cluster==null) {
				i++;
				continue;
			}
			LocalHomologyCluster localCluster1 = cluster.getLocalCluster(unitG1);
			if(localCluster1==null) {
				i++;
				continue;
			}
			List<LocalHomologyCluster> localClusters2 = cluster.getLocalClusters(genomeId2, seqName2);
			if(localClusters2==null || localClusters2.size()==0) {
				i++;
				continue;
			}
			int leftLimit = (lastInSynteny2!=null)?lastInSynteny2.getLast():0;
			int rightLimit = nextInSynteny2.getFirst();
			if(reverse) {
				leftLimit = nextInSynteny2.getLast();
				rightLimit = (lastInSynteny2!=null)?lastInSynteny2.getFirst():Integer.MAX_VALUE;
			}
			for(LocalHomologyCluster c2:localClusters2) {
				if((c2.getFirst()>leftLimit) && c2.getLast()<rightLimit) {
					SyntenyVertex middleVertex = new SyntenyVertex(localCluster1, c2);
					answer.add(middleVertex);
					lastInSynteny2 = c2;
					break;
				}
			}
			i++;
			
		}
		//TODO: End of sequence
		for(;j<vertices.size();j++) answer.add(vertices.get(j));
		return answer;
	}
	

}
